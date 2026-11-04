"""Helper for 04_evaluate_and_plot.ipynb: score every available model identically, then
write the manuscript tables and figures.

Inputs : results/{EP}/inputs/ and results/{EP}/<model>/predictions.csv (common format)
Outputs: results/{EP}/evaluation/*.csv and plots/{EP}/Figure*.png|pdf

A model whose predictions.csv does not exist is skipped and reported; all figures and tables
adapt to the models present.
"""
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from lifelines import KaplanMeierFitter
from sksurv.nonparametric import CensoringDistributionEstimator

import landmark_common as lc
from landmark_common import (EVAL_GRID, HORIZONS, HISTORY_LENGTHS, IBS_NAME, MIN_CASES_AUC,
                             MIN_CONTROLS_AUC, MIN_G, PRIMARY_AUC_HORIZON, PRIMARY_HISTORY,
                             REFERENCE_HISTORY, RISK_HORIZON, CALIBRATION_HORIZON, R_MODELS,
                             SURVIVAL_COLUMNS, evaluate, survival_labels)

MODEL_LABELS = {"Cox_ridge": "Landmark Cox", "DeepSurv": "DeepSurv", "DynamicDeepHit": "Dynamic-DeepHit",
                "JointModel": "Joint model (JMbayes2)",
                "Null_KM": "No-covariate model (Kaplan–Meier)"}
NULL_MODEL, NULL_HISTORY = "Null_KM", "no_covariates"
COLORS = {"Cox_ridge": "#0072B2", "DeepSurv": "#D55E00", "DynamicDeepHit": "#CC79A7",
          "JointModel": "#E69F00"}
HISTORY_LABELS = {"baseline_visit": "First visit", "1_visit": "Latest 1", "2_visits": "Latest 2",
                  "3_visits": "Latest 3", "all_visits": "All prior", "basic_only": "Baseline features only",
                  "no_covariates": "No covariates"}
COHORT_LABELS = {"NACC_internal": "NACC (internal, out-of-fold)", "AIBL_external": "AIBL (external)"}
# Figure 2 axis: the published ADRS paper's terms in brackets.
FIGURE2_LABELS = {"baseline_visit": "First visit\n(init)", "1_visit": "Latest 1\n(last)", "2_visits": "Latest 2",
                  "3_visits": "Latest 3", "all_visits": "All prior\n(aggregated)"}
FEATURE_LABELS = {"age": "Age at first visit", "sex": "Sex", "edu": "Education", "APOE4": "APOE4", "CDR_SB": "CDR-SB", "GDS": "GDS",
                  "MMSE": "MMSE", "memory_recall1": "Logical memory I", "memory_recall2": "Logical memory II",
                  "Digit_Span_F": "Digit span forward", "Digit_Span_B": "Digit span backward",
                  "animal": "Animal fluency", "BNT": "Boston Naming Test", "BMI": "BMI",
                  "blood_pressure_diastolic": "Diastolic BP", "blood_pressure_systolic": "Systolic BP",
                  "Heart_rate": "Heart rate", "visit_timing": "Visit timing"}
METRIC_LABELS = {IBS_NAME: "IBS, 24–60 months (lower = better)",
                 "AUC_60": "AUC at 60 months (higher = better)"}
FIXED_RISK_THRESHOLD = None  # optional prespecified absolute-risk probability; default = training median
RISK_INDEX = int(np.flatnonzero(EVAL_GRID == RISK_HORIZON)[0])
GROUP_COLORS = {"Lower risk": "#009E73", "Higher risk": "#D55E00"}


# ---------------------------------------------------------------- loading
def load_inputs(endpoint, mode):
    folder = lc.inputs_dir(endpoint, mode)
    config = json.loads((folder / "config.json").read_text())
    cohort = pd.read_csv(folder / "landmark_cohort.csv", dtype={"ID": str})
    folds = pd.read_csv(folder / "outer_folds.csv", dtype={"ID": str})
    tables = {co: cohort[cohort.cohort == co].reset_index(drop=True) for co in ["NACC", "AIBL"]}
    fold_of = tables["NACC"].ID.map(dict(zip(folds.ID, folds.outer_test_fold)))
    if fold_of.isna().any():
        raise ValueError("Every NACC participant needs an outer fold.")
    splits = [(np.flatnonzero(fold_of.ne(k)), np.flatnonzero(fold_of.eq(k)))
              for k in sorted(fold_of.unique())]
    return config, tables, splits


def load_predictions(endpoint, mode, models, config):
    found, missing = {}, []
    for model in models:
        folder = lc.model_dir(endpoint, mode, model)
        if not (folder / "predictions.csv").exists():
            missing.append(model)
            continue
        info = json.loads((folder / "run_info.json").read_text())
        if info["shared_config"]["file_sha256"] != config["file_sha256"]:
            raise ValueError(f"{model}: predictions were made from different shared inputs; refit this model.")
        found[model] = pd.read_csv(folder / "predictions.csv", dtype={"ID": str})
    for model in missing:
        print(f"SKIPPED {model}: no predictions in {lc.model_dir(endpoint, mode, model)}")
    return found


def threshold_from_survival(training_survival):
    # No outcome-based cutpoint search, external percentiles, or external recalibration.
    if FIXED_RISK_THRESHOLD is not None:
        return float(FIXED_RISK_THRESHOLD), "prespecified_absolute_risk"
    return float(np.median(1 - training_survival[:, RISK_INDEX])), "training_participant_median_exploratory"


# ---------------------------------------------------------------- scoring
class Evaluation:
    """All metrics for one endpoint; attributes mirror the CSV files written."""

    def __init__(self, endpoint, mode, include_models, n_bootstrap=None):
        self.endpoint, self.mode = endpoint, mode
        self.n_bootstrap = n_bootstrap or lc.mode_settings(mode)["n_bootstrap"]
        self.config, self.tables, self.splits = load_inputs(endpoint, mode)
        self.predictions = load_predictions(endpoint, mode, include_models, self.config)
        self.models = [m for m in include_models if m in self.predictions]
        if "Cox_ridge" not in self.models:
            raise ValueError("Landmark Cox predictions are required (reference for model comparisons).")
        plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                             "savefig.facecolor": "white", "figure.facecolor": "white", "pdf.fonttype": 42})
        self.out = lc.ensure(lc.evaluation_dir(endpoint, mode))
        self.fig_dir = lc.ensure(lc.figure_dir(endpoint, mode))
        self.results = {"NACC_internal": {}, "AIBL_external": {}}
        self.summary_rows, self.curve_parts, self.fold_rows = [], [], []
        self.strat_parts, self.threshold_rows = [], []

    # One model/history: out-of-fold NACC predictions and locked AIBL predictions.
    def _assemble(self, model, history, pred):
        development, external = self.tables["NACC"], self.tables["AIBL"]

        def survival(fit, part, ids):
            s = pred[(pred.fit == fit) & (pred.set == part)].set_index("ID")
            ids = list(ids)
            if part == "train" and not has_train:
                return None
            if len(s) != len(ids) or set(s.index) != set(ids):
                raise ValueError(f"{model}/{history}: {fit}/{part} participants differ from the shared folds.")
            return s.loc[ids, SURVIVAL_COLUMNS].to_numpy(float)

        # Training-set predictions set the risk-group threshold. A model may omit them for some histories
        # (e.g. when they are costly); those histories are then scored without risk groups.
        has_train = (pred.set == "train").any()
        oof = np.full((len(development), len(EVAL_GRID)), np.nan)
        oof_threshold = np.full(len(development), np.nan)
        oof_fold = np.zeros(len(development), dtype=int)
        for fold, (tr, te) in enumerate(self.splits, 1):
            prediction = survival(f"outer_{fold}", "test", development.ID.iloc[te])
            oof[te], oof_fold[te] = prediction, fold
            if has_train:
                threshold, rule = threshold_from_survival(survival(f"outer_{fold}", "train", development.ID.iloc[tr]))
                oof_threshold[te] = threshold
                self.threshold_rows.append(dict(model=model, history=history, fit=f"outer_{fold}", threshold=threshold,
                                                rule=rule, training_n=len(tr), horizon=RISK_HORIZON))
            score, _ = evaluate(development.iloc[te], prediction)
            self.fold_rows.append({"model": model, "history": history, "fold": fold, **score})
        threshold, rule = np.nan, "no training-set predictions"
        if has_train:
            threshold, rule = threshold_from_survival(survival("full_NACC", "train", development.ID))
            self.threshold_rows.append(dict(model=model, history=history, fit="full_NACC", threshold=threshold,
                                            rule=rule, training_n=len(development), horizon=RISK_HORIZON))
        external_prediction = survival("full_NACC", "test", external.ID)
        for cohort, table, prediction, thresholds, folds in [
                ("NACC_internal", development, oof, oof_threshold, oof_fold),
                ("AIBL_external", external, external_prediction, np.full(len(external), threshold),
                 np.zeros(len(external), int))]:
            self.results[cohort][(history, model)] = prediction
            score, curve = evaluate(table, prediction)
            meta = dict(cohort=cohort, model=model, history=history)
            self.summary_rows.append({**meta, **score})
            self.curve_parts.append(curve.assign(**meta))
            if not has_train:
                continue
            export = table[["ID", "time", "event"]].copy()
            for h in HORIZONS:
                export[f"risk_{int(h)}m"] = 1 - prediction[:, np.flatnonzero(EVAL_GRID == h)[0]]
            export["risk_threshold"] = thresholds
            export["risk_group"] = np.where(1 - prediction[:, RISK_INDEX] >= thresholds, "Higher risk", "Lower risk")
            export["outer_fold"] = folds
            self.strat_parts.append(export.assign(**meta, threshold_rule=rule))

    def _null_predictions(self):
        """Standard reference: the same Kaplan–Meier curve (fitted on the training set) for everyone."""
        development, external = self.tables["NACC"], self.tables["AIBL"]
        parts = []

        def km(train):
            return KaplanMeierFitter().fit(train.time, train.event).survival_function_at_times(EVAL_GRID).to_numpy()

        fits = [(f"outer_{k}", development.iloc[tr], development.iloc[te], "NACC")
                for k, (tr, te) in enumerate(self.splits, 1)] + [("full_NACC", development, external, "AIBL")]
        for fit, train, test, test_cohort in fits:
            curve = km(train)
            for part, table, cohort in [("train", train, "NACC"), ("test", test, test_cohort)]:
                parts.append(lc.prediction_frame(table.ID, cohort, fit, part, NULL_MODEL, NULL_HISTORY,
                                                 np.tile(curve, (len(table), 1))))
        return pd.concat(parts, ignore_index=True)

    def score(self):
        for model in self.models:
            pred = self.predictions[model]
            for history in pred.history.unique():
                print(f"{self.endpoint}: scoring {model} / {history}", flush=True)
                self._assemble(model, history, pred[pred.history == history])
        self._assemble(NULL_MODEL, NULL_HISTORY, self._null_predictions())
        self.summary = pd.DataFrame(self.summary_rows)
        self.curves = pd.concat(self.curve_parts, ignore_index=True)
        self.stratification = pd.concat(self.strat_parts, ignore_index=True)
        self.summary.to_csv(self.out / "metrics_point_estimates.csv", index=False)
        self.curves.to_csv(self.out / "metric_curves.csv", index=False)
        pd.DataFrame(self.fold_rows).to_csv(self.out / "internal_fold_metrics.csv", index=False)
        self.stratification.to_csv(self.out / "risk_group_assignments.csv", index=False)
        pd.DataFrame(self.threshold_rows).to_csv(self.out / "risk_thresholds.csv", index=False)
        return self

    # Paired participant bootstrap: identical resamples for every model/history.
    def bootstrap(self):
        metrics = [IBS_NAME] + [f"AUC_{int(h)}" for h in HORIZONS] + [f"Brier_{int(h)}" for h in HORIZONS]
        reference_key = (PRIMARY_HISTORY, "Cox_ridge")
        marginal, history_rows, model_rows = [], [], []
        for cohort, source in [("NACC_internal", "NACC"), ("AIBL_external", "AIBL")]:
            table, predictions = self.tables[source], self.results[cohort]
            samples = {key: {m: [] for m in metrics} for key in predictions}
            history_keys = [k for k in predictions if k[0] in HISTORY_LENGTHS and k[0] != REFERENCE_HISTORY
                            and (REFERENCE_HISTORY, k[1]) in predictions]
            model_keys = [(PRIMARY_HISTORY, m) for m in self.models
                          if m in lc.PYTHON_MODELS and m != "Cox_ridge" and (PRIMARY_HISTORY, m) in predictions]
            model_keys += [(self._representative(m), m) for m in self.models if m in R_MODELS]
            rng = np.random.default_rng(lc.SEED + 200)
            for b in range(self.n_bootstrap):
                idx = rng.integers(0, len(table), len(table))
                for key, prediction in predictions.items():
                    try:
                        score, _ = evaluate(table.iloc[idx].reset_index(drop=True), prediction[idx],
                                            include_concordance=False)
                    except ValueError:
                        score = {m: np.nan for m in metrics}
                    for m in metrics:
                        samples[key][m].append(score[m])
                if (b + 1) % 50 == 0:
                    print(f"{self.endpoint} {cohort}: bootstrap {b + 1}/{self.n_bootstrap}", flush=True)
            points = self.summary[self.summary.cohort == cohort].set_index(["history", "model"])
            for key in predictions:
                for m in metrics:
                    marginal.append(_interval(np.array(samples[key][m]), self.n_bootstrap, cohort=cohort,
                                              history=key[0], model=key[1], metric=m, estimate=points.loc[key, m]))
            for keys, ref, rows, label in [(history_keys, None, history_rows, "history minus latest 1 visit (same model)"),
                                           (model_keys, reference_key, model_rows, "model minus Landmark Cox (latest 3)")]:
                for key in keys:
                    base = ref or (REFERENCE_HISTORY, key[1])
                    for m in metrics:
                        diff = np.array(samples[key][m]) - np.array(samples[base][m])
                        rows.append(_interval(diff, self.n_bootstrap, cohort=cohort, history=key[0], model=key[1],
                                              metric=m, contrast=label,
                                              estimate=points.loc[key, m] - points.loc[base, m]))
                    # Relative change in IBS (%), bootstrapped as a ratio within each resample.
                    k_ibs, b_ibs = np.array(samples[key][IBS_NAME]), np.array(samples[base][IBS_NAME])
                    rows.append(_interval(100 * (k_ibs - b_ibs) / b_ibs, self.n_bootstrap, cohort=cohort,
                                          history=key[0], model=key[1], metric="IBS_change_pct", contrast=label,
                                          estimate=100 * (points.loc[key, IBS_NAME] - points.loc[base, IBS_NAME])
                                          / points.loc[base, IBS_NAME]))
        self.marginal = pd.DataFrame(marginal)
        self.history_diff = pd.DataFrame(history_rows)
        self.model_diff = pd.DataFrame(model_rows)
        self.marginal.to_csv(self.out / "metrics_with_95CI.csv", index=False)
        self.history_diff.to_csv(self.out / "paired_history_differences.csv", index=False)
        self.model_diff.to_csv(self.out / "paired_model_differences.csv", index=False)
        self.primary = self.history_diff[(self.history_diff.history == PRIMARY_HISTORY)
                                         & self.history_diff.metric.isin([IBS_NAME, f"AUC_{int(PRIMARY_AUC_HORIZON)}"])]
        self.primary.to_csv(self.out / "primary_3_vs_1_visits.csv", index=False)
        self.history_diff[self.history_diff.history == "baseline_visit"].to_csv(
            self.out / "recency_first_vs_latest_visit.csv", index=False)
        return self

    def tables_for_manuscript(self):
        cols = [IBS_NAME] + [f"AUC_{int(h)}" for h in HORIZONS] + [f"Brier_{int(RISK_HORIZON)}"]
        history_order = list(HISTORY_LENGTHS) + ["basic_only", NULL_HISTORY]

        def shaped(frame, index, cols=cols):
            frame = frame[frame.metric.isin(cols)].copy()
            frame["value"] = [(f"{e:+.1f}% ({lo:+.1f}% to {hi:+.1f}%)" if m == "IBS_change_pct" else
                               f"{e:.3f} ({lo:.3f} to {hi:.3f})") if np.isfinite(lo) else f"{e:.3f} (unsupported)"
                              for m, e, lo, hi in zip(frame.metric, frame.estimate, frame.lower, frame.upper)]
            # The no-covariate model gives everyone the same risk: AUC is not defined.
            frame.loc[frame.model.eq(NULL_MODEL) & frame.metric.str.startswith("AUC"), "value"] = "not applicable"
            frame["model"] = pd.Categorical(frame.model, lc.ALL_MODELS + [NULL_MODEL], ordered=True)
            frame["history"] = pd.Categorical(frame.history, history_order, ordered=True)
            return frame.pivot_table(index=index, columns="metric", values="value", aggfunc="first",
                                     observed=True)[cols]

        self.table2 = shaped(self.marginal, ["cohort", "model", "history"])
        counts = {f"{src}_{kind}": (len(t), int((t.event.eq(1) & t.time.le(RISK_HORIZON)).sum()))
                  for src, t in self.tables.items() for kind in ["internal" if src == "NACC" else "external"]}
        cohorts = self.table2.index.get_level_values("cohort")
        self.table2.insert(0, f"events_by_{int(RISK_HORIZON)}m", [counts[c][1] for c in cohorts])
        self.table2.insert(0, "n", [counts[c][0] for c in cohorts])
        self.table3 = shaped(pd.concat([self.history_diff, self.model_diff]), ["cohort", "contrast", "model", "history"],
                             cols=[IBS_NAME, "IBS_change_pct"] + cols[1:])
        self.table2.to_csv(self.out / "Table2_performance_by_model_and_history.csv")
        self.table3.to_csv(self.out / "Table3_paired_differences.csv")
        return self

    def calibration_and_groups(self):
        rows = []
        for cohort, source in [("NACC_internal", "NACC"), ("AIBL_external", "AIBL")]:
            table = self.tables[source]
            for (history, model), prediction in self.results[cohort].items():
                for horizon in HORIZONS:
                    risk = 1 - prediction[:, np.flatnonzero(EVAL_GRID == horizon)[0]]
                    groups = (np.zeros(len(risk), dtype=int) if np.ptp(risk) < 1e-12 else
                              pd.qcut(pd.Series(risk), 3, labels=False, duplicates="drop").fillna(0).to_numpy(int))
                    for group in np.unique(groups):
                        mask = groups == group
                        sub = table.iloc[np.flatnonzero(mask)]
                        cases = int((sub.event.eq(1) & sub.time.le(horizon)).sum())
                        controls = int(sub.time.gt(horizon).sum())
                        record = dict(cohort=cohort, model=model, history=history, horizon=horizon, group=int(group) + 1,
                                      n=len(sub), cases=cases, controls=controls, predicted=float(risk[mask].mean()),
                                      observed=np.nan, lower=np.nan, upper=np.nan,
                                      status="insufficient events or follow-up")
                        if cases >= MIN_CASES_AUC and controls >= MIN_CONTROLS_AUC:
                            km = KaplanMeierFitter().fit(sub.time, sub.event)
                            censor = CensoringDistributionEstimator().fit(survival_labels(sub))
                            if float(censor.predict_proba(np.array([horizon]))[0]) >= MIN_G:
                                ci = km.confidence_interval_.loc[lambda d: d.index <= horizon].iloc[-1]
                                record.update(observed=1 - float(km.predict(horizon)), lower=1 - float(ci.iloc[1]),
                                              upper=1 - float(ci.iloc[0]), status="ok")
                        rows.append(record)
        self.calibration = pd.DataFrame(rows)
        self.calibration.to_csv(self.out / "calibration_groups.csv", index=False)
        group_rows, at_risk_rows = [], []
        for (cohort, model, history, group), s in self.stratification.groupby(["cohort", "model", "history", "risk_group"]):
            km = KaplanMeierFitter().fit(s.time, s.event)
            controls = int(s.time.gt(RISK_HORIZON).sum())
            cases = int((s.event.eq(1) & s.time.le(RISK_HORIZON)).sum())
            supported = controls >= MIN_CONTROLS_AUC and cases >= MIN_CASES_AUC
            if supported:
                supported = float(CensoringDistributionEstimator().fit(survival_labels(s))
                                  .predict_proba(np.array([RISK_HORIZON]))[0]) >= MIN_G
            ci = km.confidence_interval_.loc[lambda d: d.index <= RISK_HORIZON].iloc[-1]
            group_rows.append(dict(cohort=cohort, model=model, history=history, risk_group=group, n=len(s),
                                   events_by_horizon=cases, observed_beyond_horizon=controls,
                                   mean_predicted_risk=s[f"risk_{int(RISK_HORIZON)}m"].mean(),
                                   observed_risk=1 - float(km.predict(RISK_HORIZON)) if supported else np.nan,
                                   lower=1 - float(ci.iloc[1]) if supported else np.nan,
                                   upper=1 - float(ci.iloc[0]) if supported else np.nan,
                                   status="ok" if supported else "insufficient events or follow-up"))
            for t in range(0, int(RISK_HORIZON) + 1, 12):
                at_risk_rows.append(dict(cohort=cohort, model=model, history=history, risk_group=group,
                                         month=t, at_risk=int(s.time.ge(t).sum())))
        pd.DataFrame(group_rows).to_csv(self.out / "risk_group_outcomes.csv", index=False)
        pd.DataFrame(at_risk_rows).to_csv(self.out / "risk_group_numbers_at_risk.csv", index=False)
        return self

    # ------------------------------------------------------------ figures
    def _representative(self, model):
        """History used for model comparisons: latest 3 visits, or the model's only history if it has one."""
        histories = list(self.predictions[model].history.unique())
        return PRIMARY_HISTORY if PRIMARY_HISTORY in histories or len(histories) != 1 else histories[0]

    def _comparison_keys(self):
        """One representative history per model (latest 3 visits) for the model-comparison figures."""
        return [(self._representative(m), m) for m in self.models]

    def _save(self, fig, name):
        if self.mode == "smoke":
            fig.suptitle("EXECUTION CHECK ONLY: not results", color="#aa2222", y=1.02)
        fig.savefig(self.fig_dir / f"{name}.png", dpi=400, bbox_inches="tight")
        fig.savefig(self.fig_dir / f"{name}.pdf", bbox_inches="tight")
        plt.show()
        plt.close(fig)

    def figures(self):
        plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                             "savefig.facecolor": "white", "figure.facecolor": "white", "pdf.fonttype": 42})
        ep = self.endpoint
        metrics = [IBS_NAME, f"AUC_{int(PRIMARY_AUC_HORIZON)}"]
        histories = list(HISTORY_LENGTHS)

        # Figure 2: performance by visit history, every model on the same axis. Reference lines are
        # standard only: AUC = 0.5 (chance) and the no-covariate Kaplan–Meier model for the IBS.
        fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), layout="constrained")
        for r, cohort in enumerate(COHORT_LABELS):
            for c, metric in enumerate(metrics):
                ax = axes[r, c]
                for j, model in enumerate(self.models):
                    s = self.marginal[(self.marginal.cohort == cohort) & (self.marginal.model == model)
                                      & (self.marginal.metric == metric)].set_index("history").reindex(histories)
                    x = np.arange(len(histories)) + (j - (len(self.models) - 1) / 2) * .08
                    marker = "D" if model in R_MODELS else "o"
                    ax.plot(x[1:], s.estimate.iloc[1:], color=COLORS[model], marker=marker, label=MODEL_LABELS[model])
                    ax.plot(x[0], s.estimate.iloc[0], color=COLORS[model], marker="s", ls="none")
                    ax.vlines(x, s.lower, s.upper, color=COLORS[model], alpha=.6)
                if metric == IBS_NAME:
                    null = self.marginal[(self.marginal.cohort == cohort) & (self.marginal.model == NULL_MODEL)
                                         & (self.marginal.metric == metric)].estimate.iloc[0]
                    ax.axhline(null, color="#888888", ls="--", lw=1, label=MODEL_LABELS[NULL_MODEL])
                else:
                    ax.axhline(.5, color="#888888", ls="--", lw=1, label="Chance (AUC = 0.5)")
                    ax.set_ylim(bottom=.45, top=1)
                ax.axvline(.5, color="#bbbbbb", ls=":", lw=1)
                ax.set_xticks(range(len(histories)), [FIGURE2_LABELS[h] for h in histories])
                ax.set(title=f"{COHORT_LABELS[cohort]}", ylabel=METRIC_LABELS[metric],
                       xlabel="Visits available at month 48")
                ax.grid(axis="y", alpha=.15)
        handles, labels = [], []
        for ax in axes[0]:
            for h, l in zip(*ax.get_legend_handles_labels()):
                if l not in labels:
                    handles.append(h)
                    labels.append(l)
        fig.legend(handles, labels, loc="outside lower center", ncol=min(len(labels), 4), frameon=False, fontsize=8)
        fig.suptitle(f"{ep}: performance by visit history (95% bootstrap CI; square = first visit)")
        self._save(fig, "Figure2_performance_by_history")

        # Supplementary Figures S1 and S2: forest plots of paired differences.
        for name, frame, keys, xlabel, better in [
                ("FigureS1_three_vs_one_visits", self.history_diff,
                 [(PRIMARY_HISTORY, m) for m in self.models if (PRIMARY_HISTORY, m) in self.results["NACC_internal"]
                  and (REFERENCE_HISTORY, m) in self.results["NACC_internal"]],
                 "Latest 3 minus latest 1 visit", "latest 3"),
                ("FigureS2_model_comparison", self.model_diff,
                 [k for k in self._comparison_keys() if k != (PRIMARY_HISTORY, "Cox_ridge")],
                 "Model minus Landmark Cox (latest 3)", "the model")]:
            if not keys:
                continue
            height = 1.6 + .45 * len(keys) * 2
            fig, axes = plt.subplots(2, 2, figsize=(10, height), layout="constrained")
            for r, cohort in enumerate(COHORT_LABELS):
                for c, metric in enumerate(metrics):
                    ax = axes[r, c]
                    s = frame[(frame.cohort == cohort) & (frame.metric == metric)]
                    for j, (history, model) in enumerate(keys):
                        v = s[(s.history == history) & (s.model == model)]
                        if v.empty:
                            continue
                        v = v.iloc[0]
                        ax.plot(v.estimate, j, "o", color=COLORS[model])
                        ax.hlines(j, v.lower, v.upper, color=COLORS[model], lw=2)
                    ax.axvline(0, color="#555555", ls="--", lw=1)
                    ax.set_yticks(range(len(keys)), [MODEL_LABELS[m] for _, m in keys])
                    favour = f"{'below' if metric == IBS_NAME else 'above'} 0 favours {better}"
                    ax.set(title=COHORT_LABELS[cohort], xlabel=f"{xlabel}: Δ {metric.split('_')[0]} ({favour})")
                    ax.margins(y=.3)
                    ax.grid(axis="x", alpha=.15)
            reference = ("; Landmark Cox (latest 3) is the reference line at 0" if name.startswith("FigureS2") else "")
            fig.suptitle(f"{ep}: {xlabel.lower()} (95% paired bootstrap CI{reference})")
            self._save(fig, name)

        # Supplementary Figure S3: AUC over 24–60 months.
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained", sharey=True)
        for ax, cohort in zip(axes, COHORT_LABELS):
            for history, model in self._comparison_keys():
                c = self.curves[(self.curves.cohort == cohort) & (self.curves.history == history)
                                & (self.curves.model == model)]
                ax.plot(c.horizon, c.AUC, color=COLORS[model], lw=1.8,
                        label=f"{MODEL_LABELS[model]} ({HISTORY_LABELS[history]})")
            ax.set(ylim=(.5, 1), xlabel="Months after the month-48 landmark", title=COHORT_LABELS[cohort])
            ax.set_xticks(np.arange(24, 61, 6))
            ax.grid(alpha=.15)
        axes[0].set_ylabel("Time-dependent AUC")
        axes[0].legend(frameon=False, fontsize=8)
        fig.suptitle(f"{ep}: discrimination over time")
        self._save(fig, "FigureS3_auc_over_time")

        keys = self._comparison_keys()
        risk_label = "net MCI risk" if ep == "MCI" else "dementia/AD risk"

        # Figure 3: Kaplan–Meier by risk group; Figure S4: calibration at 60 months;
        # Figure S5: predicted-risk distribution. Rows = cohorts, columns = models.
        def grid():
            return plt.subplots(2, len(keys), figsize=(3.2 * len(keys), 6.4), layout="constrained", squeeze=False)

        fig6, ax6 = grid()
        # Kaplan–Meier grid without constrained layout: room below each panel for numbers at risk.
        fig7, ax7 = plt.subplots(2, len(keys), figsize=(3.4 * len(keys), 8.6), squeeze=False)
        fig7.subplots_adjust(hspace=.95, wspace=.5, bottom=.16, top=.9, left=.11, right=.98)
        at_risk_ticks = np.arange(0, RISK_HORIZON + 1, 12)
        figs1, axs1 = grid()
        for r, cohort in enumerate(COHORT_LABELS):
            for c, (history, model) in enumerate(keys):
                title = f"{MODEL_LABELS[model]}\n{COHORT_LABELS[cohort]}"
                s = self.calibration[(self.calibration.cohort == cohort) & (self.calibration.history == history)
                                     & (self.calibration.model == model) & (self.calibration.horizon == CALIBRATION_HORIZON)]
                good = s[s.status == "ok"]
                ax = ax6[r, c]
                ax.plot([0, 1], [0, 1], "--", color="#999999", lw=1)
                ax.scatter(good.predicted, good.observed, color=COLORS[model], zorder=3)
                ax.vlines(good.predicted, good.lower, good.upper, color=COLORS[model])
                if len(good) < len(s):
                    ax.text(.03, .96, f"{len(s) - len(good)} of {len(s)} groups unsupported",
                            transform=ax.transAxes, va="top", fontsize=7)
                ax.set(xlim=(0, 1), ylim=(0, 1), title=title, xlabel=f"Predicted 60-month {risk_label}",
                       ylabel="Observed (Kaplan–Meier)")
                ax.xaxis.set_major_formatter(PercentFormatter(1))
                ax.yaxis.set_major_formatter(PercentFormatter(1))
                st = self.stratification[(self.stratification.cohort == cohort) & (self.stratification.history == history)
                                         & (self.stratification.model == model)]
                ax = ax7[r, c]
                for group, color in GROUP_COLORS.items():
                    sub = st[st.risk_group == group]
                    if len(sub):
                        KaplanMeierFitter(label=f"{group} (n={len(sub)})").fit(sub.time, sub.event).plot_survival_function(
                            ax=ax, ci_show=True, color=color, show_censors=False)
                ax.set(xlim=(0, RISK_HORIZON), ylim=(0, 1), title=title, xlabel="Months after month 48",
                       ylabel="Event-free survival")
                ax.set_xticks(np.arange(0, RISK_HORIZON + 1, 12))
                ax.yaxis.set_major_formatter(PercentFormatter(1))
                ax.legend(frameon=False, fontsize=7, loc="lower left")
                ax.text(-.3, -.3, "Number at risk", transform=ax.transAxes, ha="left", fontsize=7, weight="bold")
                for i, (group, color) in enumerate(GROUP_COLORS.items()):
                    sub = st[st.risk_group == group]
                    y = -.39 - .09 * i
                    ax.text(-.3, y, group.split()[0], transform=ax.transAxes, ha="left", fontsize=7, color=color)
                    for t in at_risk_ticks:
                        ax.text(t, y, f"{int(sub.time.ge(t).sum())}", transform=ax.get_xaxis_transform(),
                                ha="center", fontsize=7, color=color)
                ax = axs1[r, c]
                for group, color in GROUP_COLORS.items():
                    v = st.loc[st.risk_group == group, f"risk_{int(RISK_HORIZON)}m"]
                    if len(v):
                        ax.hist(v, bins=np.linspace(0, 1, 26), alpha=.65, color=color, label=f"{group} (n={len(v)})")
                ax.set(xlim=(0, 1), title=title, xlabel=f"Predicted 60-month {risk_label}", ylabel="Participants")
                ax.xaxis.set_major_formatter(PercentFormatter(1))
                ax.legend(frameon=False, fontsize=7)
        fig6.suptitle(f"{ep}: calibration at 60 months after the landmark (risk tertiles, 95% CI)")
        fig7.suptitle(f"{ep}: event-free survival by predicted-risk group (training-median threshold)")
        figs1.suptitle(f"{ep}: distribution of predicted 60-month risk")
        self._save(fig7, "Figure3_risk_groups")
        self._save(fig6, "FigureS4_calibration")
        self._save(figs1, "FigureS5_risk_distribution")
        self._feature_importance_figure()
        (self.out / "figure_notes.txt").write_text(
            "Figure1: participant flow from the source data to the month-48 landmark cohort.\n"
            "Figure2: points = estimate; bars = 95% participant bootstrap CI conditional on fitted predictions; "
            "square = first visit; dashed lines = chance (AUC 0.5) and the no-covariate "
            "Kaplan–Meier model (IBS). Brackets: terms of the published ADRS paper.\n"
            "Figure3: thresholds are training-set medians of predicted 60-month risk, locked for AIBL; exploratory. "
            "Numbers at risk below each panel (also in risk_group_numbers_at_risk.csv).\n"
            "FigureS1/S2: paired bootstrap differences with identical resamples across models.\n"
            "FigureS4: grouped Kaplan–Meier observed risk with Greenwood limits; unsupported groups omitted.\n"
            "FigureS6: permutation importance on held-out NACC folds (latest 3 visits); mean and SD over folds and repeats.\n")
        return self

    def _feature_importance_figure(self):
        parts = []
        for model in [m for m in self.models if m in lc.PYTHON_MODELS]:
            path = lc.model_dir(self.endpoint, self.mode, model) / "feature_importance.csv"
            if path.exists():
                parts.append(pd.read_csv(path))
        if not parts:
            print("No feature_importance.csv found: rerun 02_python_models_*.ipynb for Figure S6.")
            return
        d = pd.concat(parts, ignore_index=True)
        summary = (d.groupby(["model", "feature"])[["auc60_drop", "ibs_increase"]]
                   .agg(["mean", "std"]).round(5))
        summary.columns = ["_".join(c) for c in summary.columns]
        summary.reset_index().to_csv(self.out / "feature_importance_summary.csv", index=False)
        models = list(dict.fromkeys(d.model))
        order = (summary.loc["Cox_ridge"].auc60_drop_mean.sort_values().index if "Cox_ridge" in models
                 else summary.loc[models[0]].auc60_drop_mean.sort_values().index)
        fig, axes = plt.subplots(1, len(models), figsize=(3.6 * len(models), 5.2), layout="constrained",
                                 sharey=True, squeeze=False)
        for ax, model in zip(axes[0], models):
            s = summary.loc[model].reindex(order)
            ax.barh(range(len(order)), s.auc60_drop_mean, xerr=s.auc60_drop_std, color=COLORS[model],
                    alpha=.85, error_kw=dict(lw=.8))
            ax.axvline(0, color="#555555", lw=.8)
            ax.set(title=MODEL_LABELS[model], xlabel="Drop in AUC at 60 months when shuffled")
            ax.grid(axis="x", alpha=.15)
        axes[0, 0].set_yticks(range(len(order)), [FEATURE_LABELS.get(f, f) for f in order])
        fig.suptitle(f"{self.endpoint}: permutation feature importance (latest 3 visits; held-out NACC folds; mean ± SD)")
        self._save(fig, "FigureS6_feature_importance")

    def _flow_figure(self, flow):
        prevalent = ("Dementia at or before the first usable visit" if self.endpoint == "AD"
                     else "MCI or dementia at or before the first usable visit")
        fig, axes = plt.subplots(1, 2, figsize=(14, 6.4))
        for ax, co in zip(axes, ["NACC", "AIBL"]):
            f = flow[co]
            ax.set_axis_off()
            ax.set(xlim=(0, 1), ylim=(0, 1))
            main = dict(boxstyle="round,pad=0.5", fc="#eef3f8", ec="#0072B2")
            side = dict(boxstyle="round,pad=0.4", fc="#f7f7f7", ec="#999999")
            boxes = [(.93, f"{co}: participants in the source data\nn = {f['source_participants']:,}"),
                     (.63, f"Usable, event-free first visit\nn = {f['usable_baseline']:,}"),
                     (.33, f"Included at month 48 (landmark)\nn = {f['included']:,}"),
                     (.08, f"Events within 60 months after month 48: {f['events_by_60m']:,}\n"
                           f"Events during follow-up: {f['events_total']:,}")]
            for y, text in boxes:
                ax.text(.3, y, text, ha="center", va="center", fontsize=9, bbox=main)
            for (y0, _), (y1, _) in zip(boxes[:-1], boxes[1:]):
                ax.annotate("", xy=(.3, y1 + .06), xytext=(.3, y0 - .06),
                            arrowprops=dict(arrowstyle="->", color="#0072B2"))
            exclusions = [(.78, f"Excluded\n• No usable predictor visit: {f['no_usable_visit']:,}\n"
                                f"• {prevalent}: {f['prevalent']:,}"),
                          (.48, f"Excluded\n• Event or end of follow-up by month 48: {f['event_or_censored_by_48']:,}\n"
                                f"• Fewer than 3 visits by month 48: {f['fewer_than_3_visits']:,}")]
            for y, text in exclusions:
                ax.annotate("", xy=(.56, y), xytext=(.3, y), arrowprops=dict(arrowstyle="->", color="#999999"))
                ax.text(.57, y, text, ha="left", va="center", fontsize=8, bbox=side)
            ax.set_title("NACC (development)" if co == "NACC" else "AIBL (external validation)", fontsize=11)
        fig.suptitle(f"{self.endpoint}: participant flow")
        self._save(fig, "Figure1_participant_flow")

    def cohort_tables(self):
        """Table 1 (landmark cohort characteristics) and the participant flow, per cohort."""
        folder = lc.inputs_dir(self.endpoint, self.mode)
        visits = pd.read_csv(folder / "visits_pre48.csv", dtype={"ID": str})
        flow = pd.read_csv(folder / "cohort_accounting.csv")
        source_dir = lc.longitudinal_dir(self.endpoint)
        diagnosis_labels = {"HC": "Cognitively unimpaired", "MCI": "MCI"}

        def mean_sd(x):
            return f"{x.mean():.1f} ({x.std():.1f})"

        def median_iqr(x):
            q = x.quantile([.5, .25, .75])
            return f"{q.iloc[0]:.1f} ({q.iloc[1]:.1f} to {q.iloc[2]:.1f})"

        def n_pct(mask, denominator):
            return f"{int(mask.sum())} ({100 * mask.sum() / denominator:.1f}%)"

        columns = {}
        for co, t in self.tables.items():
            v = visits[visits.cohort == co].groupby("ID").visit_month
            n_visits, span = v.size(), v.max() - v.min()
            follow = KaplanMeierFitter().fit(t.time, 1 - t.event).median_survival_time_
            apoe = t.APOE4.dropna()
            rows = {"Participants, n": f"{len(t)}"}
            for h in HORIZONS:
                rows[f"Events within {int(h)} months after month 48, n (%)"] = n_pct(t.event.eq(1) & t.time.le(h), len(t))
            rows[f"Events during follow-up (to {int(lc.ADMIN_HORIZON)} months), n (%)"] = n_pct(t.event.eq(1), len(t))
            rows["Median follow-up after month 48, months (reverse Kaplan–Meier)"] = f"{follow:.1f}"
            rows["Age at first visit, years, mean (SD)"] = mean_sd(t.age)
            rows["Female, n (%)"] = n_pct(t.sex.eq(2), len(t))
            rows["Education, years, mean (SD)"] = mean_sd(t.edu)
            rows["APOE4 carrier, n (%)"] = n_pct(apoe.eq(1), len(apoe)) + (f"; missing {t.APOE4.isna().sum()}" if t.APOE4.isna().any() else "")
            rows["Visits by month 48, median (IQR)"] = median_iqr(n_visits)
            for k, mask in [("3", n_visits.eq(3)), ("4", n_visits.eq(4)), ("5 or more", n_visits.ge(5))]:
                rows[f"  {k} visits by month 48, n (%)"] = n_pct(mask, len(n_visits))
            long = pd.read_csv(source_dir / f"longitudinal_{co}_{self.endpoint}.csv", dtype={"ID": str},
                               usecols=["ID", "Classification", "start"])
            long = long[long.ID.isin(set(t.ID)) & long.start.le(lc.LANDMARK_MONTH)].sort_values(["ID", "start"])
            for label, rows_dx in [("first visit", long.groupby("ID").head(1)), ("latest visit by month 48", long.groupby("ID").tail(1))]:
                for code, name in diagnosis_labels.items():
                    rows[f"Diagnosis at {label}: {name}, n (%)"] = n_pct(rows_dx.Classification.eq(code), len(rows_dx))
            rows["Months from first to latest visit, median (IQR)"] = median_iqr(span)
            rows["Months from latest visit to month 48, median (IQR)"] = median_iqr(t.months_since_last_measurement)
            rows["Latest CDR-SB, mean (SD)"] = mean_sd(t.CDR_SB)
            rows["Latest CDR-SB > 0, n (%)"] = n_pct(t.CDR_SB.gt(0), t.CDR_SB.notna().sum())
            rows["Latest MMSE, median (IQR)"] = median_iqr(t.MMSE)
            columns[f"{co} ({'development' if co == 'NACC' else 'external validation'})"] = rows
        self.table1 = pd.DataFrame(columns)
        self.table1.index.name = "Characteristic"
        self.table1.to_csv(self.out / "Table1_cohort_characteristics.csv")
        flow = flow.set_index("cohort")
        counts = {}
        for co in ["NACC", "AIBL"]:
            audit = pd.read_csv(source_dir / f"outcome_audit_{co}_{self.endpoint}.csv", dtype={"ID": str})
            usable = int(audit.status.eq("included").sum())
            assert usable == flow.loc[co, "source_participants"], "Flow counts disagree with the landmark export."
            counts[co] = {"source_participants": len(audit),
                          "no_usable_visit": int(audit.status.eq("no_usable_predictor_visit").sum()),
                          "prevalent": int(audit.status.eq("endpoint_or_competing_stop_by_baseline").sum()),
                          "usable_baseline": usable,
                          "event_or_censored_by_48": int(flow.loc[co, "event_or_censored_by_landmark"]),
                          "fewer_than_3_visits": int(flow.loc[co, "insufficient_history"]),
                          "included": int(flow.loc[co, "eligible_participants"]),
                          "events_by_24m": int(flow.loc[co, "events_by_24m"]),
                          "events_by_36m": int(flow.loc[co, "events_by_36m"]),
                          "events_by_60m": int(flow.loc[co, "events_by_60m"]),
                          "events_total": int(flow.loc[co, "events_by_admin_horizon"])}
        self.flow = pd.DataFrame(counts)
        self.flow.index.name = "Participants"
        self.flow.to_csv(self.out / "TableS1_participant_flow.csv")
        self._flow_figure(self.flow)
        return self

    def write_manifest(self):
        manifest = dict(endpoint=self.endpoint, run_mode=self.mode, models=self.models,
                        created=datetime.now(timezone.utc).isoformat(), n_bootstrap=self.n_bootstrap,
                        shared_config_file_sha256=self.config["file_sha256"],
                        primary_metric=IBS_NAME, comparison_reference=f"Cox_ridge / {PRIMARY_HISTORY}",
                        risk_threshold_rule="training-set median of predicted 60-month risk; locked for AIBL")
        (self.out / "evaluation_manifest.json").write_text(json.dumps(manifest, indent=2))
        return self


def _interval(values, n_requested, **meta):
    valid = values[np.isfinite(values)]
    lo, hi = (np.quantile(valid, [.025, .975]) if len(valid) >= max(5, int(np.ceil(.8 * n_requested)))
              else (np.nan, np.nan))
    return dict(**meta, lower=lo, upper=hi, valid_replicates=len(valid), requested_replicates=n_requested)


def run(endpoint, mode, include_models, n_bootstrap=None):
    return (Evaluation(endpoint, mode, include_models, n_bootstrap)
            .cohort_tables().score().bootstrap().tables_for_manuscript().calibration_and_groups().figures()
            .write_manifest())
