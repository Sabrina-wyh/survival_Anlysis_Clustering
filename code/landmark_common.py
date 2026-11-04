"""Shared configuration and helpers for the month-48 landmark analysis.

Used by 01_prepare_inputs.ipynb, 02_python_models_{AD,MCI}.ipynb and 04_evaluate_and_plot.ipynb.

Layout (standard mode; smoke mode writes the same tree under results/smoke/ and plots/smoke/):
    results/preprocessing/{AD,MCI}/              endpoint-safe longitudinal data (Step 1)
    results/{AD,MCI}/inputs/                     shared cohort, folds, visits, grid (Step 2)
    results/{AD,MCI}/<model>/predictions.csv     every model, common format (Steps 3-4)
    results/{AD,MCI}/evaluation/                 metrics, tables, calibration, risk groups (Step 5)
    plots/{AD,MCI}/Figure<k>_<name>.png|pdf      figures
"""
import os
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sksurv.util import Surv
from sksurv.nonparametric import CensoringDistributionEstimator
from sksurv.metrics import concordance_index_ipcw

PROJECT_DIR = Path(os.environ.get("SURVIVAL_PROJECT_DIR", "/Users/yihan/Documents/Github/survival_Anlysis_Clustering"))
CODE_DIR = PROJECT_DIR / "code"
ENDPOINTS = ["AD", "MCI"]
ENDPOINT_LABELS = {"AD": "recorded dementia/AD endpoint", "MCI": "net first-MCI risk"}

# ---------------------------------------------------------------- study design
LANDMARK_MONTH = 48.0
MIN_PRIOR_VISITS = 3
HISTORY_LENGTHS = {"baseline_visit": 1, "1_visit": 1, "2_visits": 2, "3_visits": 3, "all_visits": None}
PRIMARY_HISTORY = "3_visits"
REFERENCE_HISTORY = "1_visit"
ADMIN_HORIZON = 120.0
HORIZONS = np.array([24., 36., 60.])
EVAL_GRID = np.arange(24., 61., 1.)
SURVIVAL_COLUMNS = [f"S_{int(t)}" for t in EVAL_GRID]
IBS_NAME = f"IBS_{int(EVAL_GRID[0])}_{int(EVAL_GRID[-1])}"
C_NAME = f"Uno_C_{int(EVAL_GRID[-1])}"
CALIBRATION_HORIZON = float(EVAL_GRID[-1])
PRIMARY_AUC_HORIZON = 60.0
RISK_HORIZON = 60.0
MIN_G = 0.05
MIN_CASES_AUC = 3
MIN_CONTROLS_AUC = 10
BASIC = ["age", "sex", "edu", "APOE4"]
COGNITIVE = ["CDR_SB", "GDS", "MMSE", "memory_recall1", "memory_recall2",
             "Digit_Span_F", "Digit_Span_B", "animal", "BNT"]
DYNAMIC = COGNITIVE + ["BMI", "blood_pressure_diastolic", "blood_pressure_systolic", "Heart_rate"]
CATEGORICAL = ["sex", "APOE4"]
PREDICTORS = ["age", "sex", "edu", "APOE4"] + DYNAMIC
SEED = int(os.environ.get("SURVIVAL_SEED", "42"))
assert np.isin(HORIZONS, EVAL_GRID).all() and EVAL_GRID[-1] < ADMIN_HORIZON
assert MIN_PRIOR_VISITS >= max(n for n in HISTORY_LENGTHS.values() if n is not None)

# Model names, history labels and where each model is fitted.
PYTHON_MODELS = ["Cox_ridge", "DeepSurv", "DynamicDeepHit"]
# R models: none at present. The joint model (JMbayes2) was attempted but did not converge (the
# logical-memory association parameters were not identifiable); its code is in code/not used/.
# A model listed here is scored like the others if its predictions.csv exists.
R_MODELS = []
ALL_MODELS = PYTHON_MODELS + R_MODELS


def run_mode():
    mode = os.environ.get("SURVIVAL_RUN_MODE", "standard")
    assert mode in {"standard", "smoke"}, mode
    return mode


def mode_settings(mode):
    """Standard = the analysis. Smoke = execution check run only by checks/smoke_test.py
    (2 folds, very short training and bootstrap; outputs under results/smoke/, plots/smoke/)."""
    standard = mode == "standard"
    return dict(outer_folds=3 if standard else 2, max_epochs=150 if standard else 3,
                patience=20 if standard else 2, n_bootstrap=200 if standard else 8,
                ddh_burn_in=10 if standard else 1, importance_repeats=5 if standard else 1,
                full_tuning_grid=standard)


# ---------------------------------------------------------------- paths
def results_root(mode):
    return PROJECT_DIR / "results" if mode == "standard" else PROJECT_DIR / "results" / "smoke"


def plots_root(mode):
    return PROJECT_DIR / "plots" if mode == "standard" else PROJECT_DIR / "plots" / "smoke"


def endpoint_dir(endpoint, mode):
    return results_root(mode) / endpoint


def inputs_dir(endpoint, mode):
    return endpoint_dir(endpoint, mode) / "inputs"


def model_dir(endpoint, mode, model):
    return endpoint_dir(endpoint, mode) / model


def evaluation_dir(endpoint, mode):
    return endpoint_dir(endpoint, mode) / "evaluation"


def figure_dir(endpoint, mode):
    return plots_root(mode) / endpoint


def ensure(path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def longitudinal_dir(endpoint):
    """Step 1 output; built (or reused if unchanged) by prepare_history_information.py."""
    from prepare_history_information import ensure_information_inputs
    return ensure_information_inputs(PROJECT_DIR) / endpoint


# ---------------------------------------------------------------- cohort (verbatim from the analysis notebook)
def load_longitudinal(path):
    d = pd.read_csv(path, dtype={"ID": str})
    needed = {"ID", "start", "stop", "event", *PREDICTORS}
    missing = needed - set(d.columns)
    if missing:
        raise ValueError(f"Missing columns in {path.name}: {missing}")
    if d[["ID", "start", "stop", "event"]].isna().any().any():
        raise ValueError("ID and survival labels must be complete.")
    for c in ["start", "stop", "event"] + PREDICTORS:
        d[c] = pd.to_numeric(d[c], errors="raise")
    if not np.isfinite(d[["start", "stop", "event"]].to_numpy()).all():
        raise ValueError("Non-finite survival labels.")
    if not d.event.isin([0, 1]).all() or (d.stop <= d.start).any():
        raise ValueError("Require binary event and stop > start.")
    if d.duplicated(["ID", "start"]).any():
        raise ValueError("Duplicated participant/measurement times.")
    d = d.sort_values(["ID", "start"]).reset_index(drop=True)
    for _, g in d.groupby("ID", sort=False):
        if not np.isclose(g.start.iloc[0], 0):
            raise ValueError("Nonzero entry time: delayed entry needs a separate design.")
        if np.any(g.start.to_numpy()[1:] < g.stop.to_numpy()[:-1] - 1e-8):
            raise ValueError("Overlapping intervals detected.")
        if g.event.sum() > 1 or (g.event.eq(1).any() and g.event.iloc[-1] != 1):
            raise ValueError("Multiple events or records after an event.")
    return d



def build_history_tables(d, landmark, minimum_visits, history_lengths):
    records = {h: [] for h in history_lengths}
    audits, excluded = [], {"event_or_censored_by_landmark": 0, "insufficient_history": 0}
    n_source = d.ID.nunique()
    for pid, g in d.groupby("ID", sort=True):
        g = g.sort_values("start")
        event_rows = g[g.event.eq(1)]
        endpoint = float(event_rows.stop.iloc[0] if len(event_rows) else g.stop.max())
        if endpoint <= landmark:
            excluded["event_or_censored_by_landmark"] += 1
            continue
        past = g[g.start <= landmark]
        if len(past) < minimum_visits:
            excluded["insufficient_history"] += 1
            continue
        last = past.iloc[-1]
        shared = {"ID": str(pid), "time": min(endpoint - landmark, ADMIN_HORIZON),
                  "event": int(bool(len(event_rows)) and endpoint - landmark <= ADMIN_HORIZON)}
        shared.update({c: last[c] for c in DYNAMIC})
        shared.update({c: g.iloc[0][c] for c in BASIC})
        shared["months_since_last_measurement"] = landmark - float(last.start)
        for label, n_visits in history_lengths.items():
            selected = past.head(1) if label == "baseline_visit" else past if n_visits is None else past.tail(n_visits)
            last = selected.iloc[-1]
            record = shared.copy()
            record.update({c: last[c] for c in DYNAMIC})
            record["months_since_last_measurement"] = landmark - float(last.start)
            t = selected.start.to_numpy(float) / 12.
            values = selected[DYNAMIC].to_numpy(float)
            for j, c in enumerate(DYNAMIC):
                valid = np.isfinite(values[:, j])
                x, y = t[valid], values[valid, j]
                slope = np.nan
                if len(y) >= 2 and np.ptp(x) > 0:
                    xc = x - x.mean()
                    slope = float(xc @ (y-y.mean()) / (xc @ xc))
                record[c+"__slope_per_year"] = slope
                record[c+"__mean_deviation"] = float(y.mean()-last[c]) if len(y) and np.isfinite(last[c]) else np.nan
            records[label].append(record)
            audits.append({"ID": str(pid), "history": label, "n_available": len(past),
                "n_used": len(selected), "first_measurement": float(selected.start.iloc[0]),
                "last_measurement": float(selected.start.iloc[-1]),
                "history_span_months": float(selected.start.iloc[-1]-selected.start.iloc[0])})
    tables = {h: pd.DataFrame(rows).replace([np.inf, -np.inf], np.nan).reset_index(drop=True)
              for h, rows in records.items()}
    reference = tables[REFERENCE_HISTORY]
    if reference.empty:
        raise ValueError("No eligible participants for this landmark/history requirement.")
    shared_cols = ["ID", "time", "event"] + PREDICTORS + ["months_since_last_measurement"]
    for h, table in tables.items():
        if table.ID.duplicated().any() or (table.time <= 0).any():
            raise ValueError("Invalid participant-level outcomes.")
        match_cols = ["ID", "time", "event"] + BASIC if h == "baseline_visit" else shared_cols
        pd.testing.assert_frame_equal(reference[match_cols], table[match_cols])
    audit = pd.DataFrame(audits)
    if (audit.last_measurement > landmark).any():
        raise ValueError("Post-landmark measurements entered predictors.")
    flow = {"source_participants": n_source, **excluded, "eligible_participants": len(reference)}
    assert n_source == sum(excluded.values()) + len(reference)
    return tables, audit, flow




def load_cohort(endpoint, write_dir=None):
    """Landmark cohort and history tables for one endpoint, with outcome-audit checks.

    Returns raw longitudinal data, history tables {cohort: {history: table}}, accounting and audit.
    """
    data_dir = longitudinal_dir(endpoint)
    raw = {co: load_longitudinal(data_dir / f"longitudinal_{co}_{endpoint}.csv") for co in ["NACC", "AIBL"]}
    datasets, audit_parts, accounting_rows = {}, [], []
    for co, d in raw.items():
        if not d.APOE4.dropna().isin([0, 1]).all():
            raise ValueError(f"{co}: APOE4 must be binary carrier status; revise input mapping explicitly.")
        datasets[co], audit, flow = build_history_tables(d, LANDMARK_MONTH, MIN_PRIOR_VISITS, HISTORY_LENGTHS)
        for history, table in datasets[co].items():
            table["__cohort"], table["__history"] = co, history
        audit_parts.append(audit.assign(cohort=co))
        reference = datasets[co][REFERENCE_HISTORY]
        row = {"cohort": co, "landmark_month": LANDMARK_MONTH, **flow,
               "events_by_admin_horizon": int(reference.event.sum())}
        row.update({f"events_by_{int(h)}m": int((reference.event.eq(1) & reference.time.le(h)).sum()) for h in HORIZONS})
        accounting_rows.append(row)
    accounting = pd.DataFrame(accounting_rows)
    history_audit = pd.concat(audit_parts, ignore_index=True)
    for co in raw:
        outcome_audit = pd.read_csv(data_dir / f"outcome_audit_{co}_{endpoint}.csv", dtype={"ID": str})
        expected = outcome_audit[outcome_audit.landmark_eligible]
        actual = datasets[co][REFERENCE_HISTORY]
        assert set(expected.ID) == set(actual.ID), "Outcome audit and model cohort disagree"
        assert not expected.target_month.le(LANDMARK_MONTH).any(), "Known target event at/before landmark"
        if endpoint == "MCI":
            assert not expected.dementia_month.le(LANDMARK_MONTH).any(), "Dementia at/before MCI landmark"
        expected = expected.set_index("ID").loc[actual.ID]
        np.testing.assert_allclose(actual.time, np.minimum(expected.endpoint_month-LANDMARK_MONTH, ADMIN_HORIZON))
        np.testing.assert_array_equal(actual.event, expected.target_observed.eq(1) & expected.endpoint_month.le(LANDMARK_MONTH+ADMIN_HORIZON))
        if write_dir is not None:
            expected.to_csv(ensure(write_dir) / f"included_outcome_audit_{co}.csv")
    if write_dir is not None:
        accounting.to_csv(ensure(write_dir) / "cohort_accounting.csv", index=False)
        history_audit.to_csv(write_dir / "history_selection_audit.csv", index=False)
    return raw, datasets, accounting, history_audit


def make_outer_folds(reference, n_folds, seed=SEED):
    """Participant-level outer folds on NACC, stratified by event (identical everywhere)."""
    from sklearn.model_selection import StratifiedKFold
    counts = reference.event.value_counts()
    if len(counts) != 2 or counts.min() < n_folds * 3:
        raise ValueError("Insufficient events/non-events for configured outer folds.")
    return list(StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed).split(reference, reference.event))


def shared_outer_folds(endpoint, mode, reference, n_folds):
    """Outer folds for the model notebooks, guaranteed identical to those exported in step 1."""
    splits = make_outer_folds(reference, n_folds)
    check_folds_against_inputs(endpoint, mode, reference, splits)
    return splits


def check_folds_against_inputs(endpoint, mode, reference, outer_splits):
    """Every model must use the exported shared folds."""
    folder = inputs_dir(endpoint, mode)
    if not (folder / "config.json").exists():
        raise FileNotFoundError(f"No shared inputs in {folder}: run 01_prepare_inputs.ipynb first.")
    shared = pd.read_csv(folder / "outer_folds.csv", dtype={"ID": str})
    mine = {reference.iloc[i].ID: k for k, (_, te) in enumerate(outer_splits, 1) for i in te}
    if dict(zip(shared.ID, shared.outer_test_fold)) != mine:
        raise ValueError("Folds differ from the shared inputs: rerun 01_prepare_inputs.ipynb.")


# ---------------------------------------------------------------- metrics (verbatim from the analysis notebook)
def survival_labels(df):
    return Surv.from_arrays(df.event.to_numpy(bool), df.time.to_numpy(float))


def metric_curve(df, survival, grid=EVAL_GRID, auc_all_times=True):
    grid = np.asarray(grid, float)
    survival = np.asarray(survival, float)
    if survival.shape != (len(df), len(grid)):
        raise ValueError("Survival array must be participants × evaluation times.")
    if not np.isfinite(survival).all() or np.any((survival < -1e-6) | (survival > 1+1e-6)):
        raise ValueError("Invalid survival probabilities.")
    if np.any(np.diff(survival, axis=1) > 1e-5):
        raise ValueError("Predicted survival must be non-increasing.")
    y = survival_labels(df)
    censor = CensoringDistributionEstimator().fit(y)
    times, events = df.time.to_numpy(), df.event.to_numpy(bool)
    records = []
    for j, horizon in enumerate(grid):
        case = events & (times <= horizon)
        control = times > horizon
        row = {"horizon": horizon, "n": len(df), "cases": int(case.sum()),
               "controls": int(control.sum()), "g_censor": np.nan, "AUC": np.nan,
               "Brier": np.nan, "status": "ok"}
        try:
            if horizon >= times.max() or not control.any():
                raise ValueError("No follow-up beyond horizon")
            gh = float(censor.predict_proba(np.array([horizon]))[0])
            row["g_censor"] = gh
            if gh < MIN_G:
                raise ValueError("Insufficient censoring support")
            case_weights = censor.predict_ipcw(y[case]) if case.any() else np.empty(0)
            # The same lower bound applies to G at observed case times.
            if len(case_weights) and (case_weights > 1. / MIN_G).any():
                raise ValueError("Unstable case censoring weights")
            row["Brier"] = float((np.sum(survival[case, j]**2 * case_weights)
                                   + np.sum((1-survival[control, j])**2 / gh)) / len(df))
            if case.sum() >= MIN_CASES_AUC and control.sum() >= MIN_CONTROLS_AUC and (auc_all_times or horizon in HORIZONS):
                label = np.r_[np.ones(case.sum()), np.zeros(control.sum())]
                risk = np.r_[1-survival[case, j], 1-survival[control, j]]
                weights = np.r_[case_weights, np.full(control.sum(), 1./gh)]
                row["AUC"] = float(roc_auc_score(label, risk, sample_weight=weights))
            else:
                row["status"] = "AUC unavailable: too few cases/controls; Brier retained"
        except ValueError as exc:
            row["status"] = str(exc)
        records.append(row)
    return pd.DataFrame(records)


def evaluate(df, survival, include_concordance=True):
    curve = metric_curve(df, survival, auc_all_times=include_concordance)
    summary = {IBS_NAME: np.nan, C_NAME: np.nan}
    if curve.Brier.notna().all():
        summary[IBS_NAME] = float(np.trapezoid(curve.Brier, EVAL_GRID) / np.ptp(EVAL_GRID))
    for h in HORIZONS:
        row = curve.loc[curve.horizon.eq(h)].iloc[0]
        summary[f"AUC_{int(h)}"] = row.AUC
        summary[f"Brier_{int(h)}"] = row.Brier
    if include_concordance:
        try:
            last = curve.iloc[-1]
            if last.g_censor < MIN_G or pd.isna(last.Brier):
                raise ValueError("Unsupported concordance horizon")
            y = survival_labels(df)
            summary[C_NAME] = float(concordance_index_ipcw(y, y, 1-survival[:, -1], tau=float(EVAL_GRID[-1]))[0])
        except ValueError:
            pass
    return summary, curve

# ---------------------------------------------------------------- common prediction format
def prediction_frame(ids, cohort, fit, part, model, history, survival):
    survival = np.asarray(survival, float)
    frame = pd.DataFrame(survival, columns=SURVIVAL_COLUMNS)
    frame.insert(0, "history", history)
    frame.insert(0, "model", model)
    frame.insert(0, "set", part)
    frame.insert(0, "fit", fit)
    frame.insert(0, "cohort", cohort)
    frame.insert(0, "ID", list(map(str, ids)))
    return frame
