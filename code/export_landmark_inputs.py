"""Helper for 01_prepare_inputs.ipynb (Step 2): export the shared month-48 landmark inputs.

Every model (Landmark Cox, DeepSurv, Dynamic-DeepHit)
uses these participants, outcomes, outer folds and evaluation grid. Written to
results/{AD,MCI}/inputs/ (standard) or results/smoke/{AD,MCI}/inputs/ (smoke).
"""
import json

import numpy as np
import pandas as pd

import landmark_common as lc


def export(endpoint, mode):
    """Write the shared inputs for one endpoint and run mode; return a one-row-per-cohort summary."""
    out = lc.ensure(lc.inputs_dir(endpoint, mode))
    settings = lc.mode_settings(mode)
    raw, datasets, accounting, history_audit = lc.load_cohort(endpoint, write_dir=out)
    reference = datasets["NACC"][lc.REFERENCE_HISTORY]
    splits = lc.make_outer_folds(reference, settings["outer_folds"])
    folds = pd.DataFrame([{"ID": reference.iloc[i].ID, "outer_test_fold": f}
                          for f, (_, te) in enumerate(splits, 1) for i in te])

    landmark, admin = lc.LANDMARK_MONTH, lc.ADMIN_HORIZON
    cohort_parts, pre_parts, dev_parts = [], [], []
    for co in ["NACC", "AIBL"]:
        table = datasets[co][lc.REFERENCE_HISTORY]
        # Outcome after the landmark, baseline features and latest pre-landmark values.
        cohort_parts.append(table[["ID", "time", "event"] + lc.BASIC + lc.DYNAMIC
                                  + ["months_since_last_measurement"]].assign(cohort=co))
        rows = raw[co][raw[co].ID.isin(set(table.ID))].copy()
        base = table.set_index("ID")[lc.BASIC]
        pre = rows[rows.start <= landmark][["ID", "start"] + lc.DYNAMIC].rename(columns={"start": "visit_month"})
        pre_parts.append(pre.join(base, on="ID").assign(cohort=co))
        if co == "NACC":
            # Development history (used by the archived joint model in checks/joint_model/): full start-stop rows from baseline,
            # administratively censored at landmark + ADMIN_HORIZON (same as the outcome).
            end = landmark + admin
            dev = rows[rows.start < end][["ID", "start", "stop", "event"] + lc.DYNAMIC].copy()
            truncated = dev.stop > end
            dev.loc[truncated, "event"] = 0
            dev.loc[truncated, "stop"] = end
            dev = dev.join(base, on="ID")
            last = dev.groupby("ID").agg(stop=("stop", "max"), event=("event", "max"))
            check = table.set_index("ID")[["time", "event"]]
            assert np.allclose(last.loc[check.index, "stop"] - landmark, check.time), "Outcome time mismatch"
            assert (last.loc[check.index, "event"].to_numpy() == check.event.to_numpy()).all(), "Event mismatch"
            dev_parts.append(dev)

    files = {
        "landmark_cohort.csv": pd.concat(cohort_parts, ignore_index=True),
        "visits_pre48.csv": pd.concat(pre_parts, ignore_index=True).sort_values(["cohort", "ID", "visit_month"]),
        "development_history.csv": pd.concat(dev_parts, ignore_index=True).sort_values(["ID", "start"]),
        "outer_folds.csv": folds,
    }
    for name, frame in files.items():
        frame.to_csv(out / name, index=False)
    grid = [float(x) for x in lc.EVAL_GRID]
    config = dict(
        endpoint=endpoint, run_mode=mode, landmark_month=landmark, administrative_horizon_months=admin,
        evaluation_grid_months=grid, survival_columns=lc.SURVIVAL_COLUMNS,
        horizons_months=[float(h) for h in lc.HORIZONS], risk_horizon_months=lc.RISK_HORIZON,
        outer_folds=settings["outer_folds"], seed=lc.SEED, basic_predictors=lc.BASIC,
        dynamic_predictors=lc.DYNAMIC, minimum_prior_visits=lc.MIN_PRIOR_VISITS,
        fits={**{f"outer_{k}": "train: NACC not in fold k; test: NACC fold k" for k in range(1, settings["outer_folds"] + 1)},
              "full_NACC": "train: all NACC; test: AIBL"},
        prediction_format="ID, cohort, fit, set (train/test), model, history, S_<month> = P(event-free at landmark + month | event-free at landmark)",
        input_sha256={co: lc.sha256_file(lc.longitudinal_dir(endpoint) / f"longitudinal_{co}_{endpoint}.csv") for co in raw},
        file_sha256={name: lc.sha256_file(out / name) for name in files})
    (out / "config.json").write_text(json.dumps(config, indent=2))

    cohort = files["landmark_cohort.csv"]
    summary = cohort.groupby("cohort", sort=False).agg(participants=("ID", "size"), events=("event", "sum"))
    for h in lc.HORIZONS:
        summary[f"events_by_{int(h)}m"] = cohort[cohort.time.le(h) & cohort.event.eq(1)].groupby("cohort").size()
    visits = files["visits_pre48.csv"].groupby(["cohort", "ID"]).size().groupby("cohort")
    summary["median_visits_by_48"] = visits.median()
    return summary.reset_index().assign(endpoint=endpoint, run_mode=mode, outer_folds=settings["outer_folds"])
