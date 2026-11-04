"""Outcome-safe inputs for the month-48 information-history analyses.

Snapshot of the user's field harmonisation; outcome construction is independent
of predictor missingness. Original preprocessing notebooks/CSVs are not modified.
"""
from pathlib import Path
import os, json, hashlib
import numpy as np
import pandas as pd

BASIC = ['age', 'sex', 'edu', 'APOE4']
DYNAMIC = ['CDR_SB', 'GDS', 'MMSE', 'memory_recall1', 'memory_recall2',
           'Digit_Span_F', 'Digit_Span_B', 'animal', 'BNT', 'BMI',
           'blood_pressure_diastolic', 'blood_pressure_systolic', 'Heart_rate']
PREDICTORS = BASIC + DYNAMIC
SCHEMA_VERSION = 'diagnosis_first_information_v1'

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def harmonize_sources(project, cache_dir):
    RAW_DIR = Path(os.environ.get('SURVIVAL_RAW_DIR', str(project/'raw_data')))
    naming = pd.concat([pd.read_excel(RAW_DIR/'naming_convension.xlsx', sheet_name=s)
                        for s in ['keys', 'Demographic', 'neuropsychological']], ignore_index=True)
    def impute_demographic_data(df):
        df["Demographic.ApoE genotype"] = df.groupby("AIBL ID")["Demographic.ApoE genotype"].transform(lambda x: x.ffill().bfill())
    
        education_mapping = {
            '0-6': (4, 6),
            '7-8': (7, 8),
            '9-12': (9, 12),
            '13-15': (13, 15),
            '15+': (16, 22)
        }
    
        def impute_years_of_education_exact(row):
            if pd.isna(row['Demographic.Years of Education Exact']) and row['Demographic.Years of Education'] in education_mapping:
                low, high = education_mapping[row['Demographic.Years of Education']]
                return np.random.default_rng(42).integers(low, high + 1) 
            return row['Demographic.Years of Education Exact']
    
        # Apply imputation function to each row
        df['Demographic.Years of Education Exact'] = df.apply(impute_years_of_education_exact, axis=1)
    
        return df
    
    def generate_visit_column(df):
        df = df.sort_values(['ID', 'age'], kind='stable').copy()
        first_age = df.groupby('ID')['age'].transform('first')
        # Preserve the original age-derived integer-month definition.
        df['VISIT'] = ((df['age'] - first_age) * 12).astype(int)
        return df.sort_values(['ID', 'VISIT'], kind='stable').reset_index(drop=True)
    
    
    aibl_columns = set(naming['AIBL'].dropna()) - {'Age'}
    aibl_columns.update(['AIBL ID', 'Demographic.YearMonthOfBirth', 'Progress Summary.Date of NP assessment',
                         'Demographic.Years of Education', 'Demographic.Years of Education Exact'])
    AIBL_raw = pd.read_excel(RAW_DIR / 'AIBL_2023.xlsx', usecols=lambda c: c in aibl_columns)
    year_month = AIBL_raw['Demographic.YearMonthOfBirth'].astype('string').str.strip().str.replace(r'\.0$', '', regex=True)
    AIBL_raw['Demographic.YearMonthOfBirth'] = pd.to_datetime(year_month + '01', format='%Y%m%d', errors='coerce')
    AIBL_raw['Progress Summary.Date of NP assessment'] = pd.to_datetime(
        AIBL_raw['Progress Summary.Date of NP assessment'], dayfirst=True, errors='coerce')
    AIBL_raw['Age'] = ((AIBL_raw['Progress Summary.Date of NP assessment'] - AIBL_raw['Demographic.YearMonthOfBirth']).dt.days / 365).round(2)
    AIBL_raw = impute_demographic_data(AIBL_raw)
    
    
    AIBL_selected = AIBL_raw[naming['AIBL'].dropna().tolist()].copy()
    AIBL_selected.replace('NI', np.nan, inplace=True)
    AIBL_selected = AIBL_selected.dropna(subset=['Age']).reset_index(drop=True)
    AIBL_rename_dict = dict(zip(naming['AIBL'], naming['New_Name']))
    AIBL_selected = AIBL_selected.rename(columns=AIBL_rename_dict)
    AIBL_selected = generate_visit_column(AIBL_selected)
    
    def add_age_column_NACC(df):
        df = df.copy()
        # Parse calendar fields before cleaning assessment-specific missing codes.
        df['VISITDATE'] = pd.to_datetime(dict(year=df['VISITYR'], month=df['VISITMO'], day=df['VISITDAY']), errors='coerce')
        df['BIRTHDATE'] = pd.to_datetime(dict(year=df['BIRTHYR'], month=df['BIRTHMO'], day=1), errors='coerce')
        df['AGE'] = ((df['VISITDATE'] - df['BIRTHDATE']).dt.total_seconds() / (365.25 * 24 * 60 * 60)).round(2)
        df['AGE'] = df['AGE'].where(df['AGE'].between(18, 120))
        if 'NACCAGE' in df:
            fallback = pd.to_numeric(df['NACCAGE'], errors='coerce')
            df['AGE'] = df['AGE'].fillna(fallback.where(fallback.between(18, 120)))
        return df
    
    # NACC UDS3 RDD: EDUC p55, NACCAGE p56; physical pp125-130;
    # CDR, GDS and C1 sections; genetic RDD NACCAPOE p6.
    # Valid 8/9 scores and valid 88/99 vital signs MUST NOT be removed globally.
    NACC_VALID_RANGES = {
        'EDUC': (0,36), 'SEX': (1,2), 'NACCAPOE': (1,6), 'NACCUDSD': (1,4),
        'NACCBMI': (10,100), 'BPSYS': (70,230), 'BPDIAS': (30,140), 'HRATE': (33,160),
        'CDRSUM': (0,18), 'NACCGDS': (0,15), 'NACCMMSE': (0,30),
        'LOGIMEM': (0,25), 'MEMUNITS': (0,25), 'DIGIF': (0,12), 'DIGIB': (0,12),
        'ANIMALS': (0,77), 'BOSTON': (0,30),
    }
    
    def clean_nacc_selected(df, naming):
        df = df.copy()
        audit = []
        ranges = dict(NACC_VALID_RANGES)
        for row in naming.itertuples(index=False):
            if str(row.New_Name).startswith('NPIQ') and pd.notna(row.NACC):
                ranges[row.NACC] = (1,3) if str(row.New_Name).endswith('_Sev') else (0,1)
        for column, (low, high) in ranges.items():
            if column not in df:
                continue
            values = pd.to_numeric(df[column], errors='coerce')
            invalid = values.notna() & ~values.between(low, high)
            audit.append({'variable':column, 'valid_min':low, 'valid_max':high,
                          'nonmissing_input':int(values.notna().sum()), 'masked_out_of_range':int(invalid.sum())})
            df[column] = values.where(values.between(low, high))
        return df, pd.DataFrame(audit)
    
    
    nacc_columns = set(naming['NACC'].dropna()) - {'AGE'}
    nacc_columns.update(['NACCID','VISITYR','VISITMO','VISITDAY','BIRTHYR','BIRTHMO','NACCAGE'])
    NACC_raw = pd.read_csv(RAW_DIR / 'investigator_ftldlbd_nacc66.csv', usecols=lambda c: c in nacc_columns, low_memory=False)
    NACC_raw = add_age_column_NACC(NACC_raw)
    NACC_raw, missing_code_audit = clean_nacc_selected(NACC_raw, naming)
    missing_code_audit.to_csv(cache_dir / 'NACC_missing_code_audit.csv', index=False)
    NACC_selected = NACC_raw[naming['NACC'].dropna().tolist()].copy()
    NACC_selected = NACC_selected.dropna(subset=['AGE']).reset_index(drop=True)
    
    
    NACC_rename_dict = dict(zip(naming['NACC'], naming['New_Name']))
    NACC_selected = NACC_selected.rename(columns=NACC_rename_dict)
    NACC_selected = generate_visit_column(NACC_selected)
    NACC_selected.head(4)
    
    NACC_selected['ID'] = NACC_selected['ID'].str.replace('NACC', '').astype(int)
    NACC_selected['APOE'] = NACC_selected['APOE'].map({1:'E3/E3' , 2:'E3/E4', 3:'E3/E2', 4:'E4/E4', 5:'E4/E2', 6:'E2/E2'})
    NACC_selected['APOE4'] = NACC_selected['APOE'].apply(lambda x: 1 if isinstance(x, str) and 'E4' in x.split('/') else (0 if isinstance(x, str) else None))
    # NACC_selected['height'] = NACC_selected['height'] * 2.54 # inches to cm
    # NACC_selected['weight'] = NACC_selected['weight'] * 0.453592 # lbs to kg
    NACC_selected['Digit_Span_F'] = (NACC_selected['Digit_Span_F'] / 12) * 16 # convert the range from 0-12 to 0-16 be consistent (later in modelling will do normalisation)
    NACC_selected['Digit_Span_B'] = (NACC_selected['Digit_Span_B'] / 12) * 14 # convert the range from 0-12 to 0-14 be consistent (later in modelling will do normalisation)
    
    NACC_selected['Classification'] = NACC_selected['Classification'].map({1:'HC' , 3:'MCI', 4:'AD'})
    NACC_selected = NACC_selected[NACC_selected['Classification'].isin(['HC', 'MCI', 'AD'])]
    NACC_selected.head(4)
    
    AIBL_selected['sex'] = AIBL_selected['sex'].map({'Female': 2, 'Male': 1})
    AIBL_selected['APOE4'] = AIBL_selected['APOE'].apply(lambda x: 1 if isinstance(x, str) and 'E4' in x.split('/') else (0 if isinstance(x, str) else None))
    
    NPIQ_cols_convert_AIBL = ['NPIQ_Delusions', 'NPIQ_Hallucinations', 'NPIQ_Agitation', 'NPIQ_Depression', 'NPIQ_Anxiety', 'NPIQ_Elation', 'NPIQ_Apathy', 'NPIQ_Disinhibition','NPIQ_Irritability', 'NPIQ_Motor', 'NPIQ_Nighttime', 'NPIQ_Appetite']
    AIBL_selected[NPIQ_cols_convert_AIBL] = AIBL_selected[NPIQ_cols_convert_AIBL].replace({'Yes': 1, 'No': 0})
    NPIQ_cols_numeric_AIBL = ['NPIQ_Hallucinations_Sev', 'NPIQ_Depression_Sev', 'NPIQ_Anxiety_Sev', 'NPIQ_Elation_Sev', 'NPIQ_Apathy_Sev', 'NPIQ_Disinhibition_Sev', 'NPIQ_Irritability_Sev', 'NPIQ_Motor_Sev', 'NPIQ_Nighttime_Sev', 'NPIQ_Appetite_Sev']  
    AIBL_selected[NPIQ_cols_numeric_AIBL] = AIBL_selected[NPIQ_cols_numeric_AIBL].apply(pd.to_numeric, errors='coerce')
    
    AIBL_selected['Classification'] = AIBL_selected['Classification'].replace('MCIX', 'MCI')
    AIBL_selected = AIBL_selected[AIBL_selected['Classification'].isin(['HC', 'MCI', 'AD'])] # Keep only rows where the 'classification' column contains 'HC', 'MCI', or 'AD'
    AIBL_selected.head(4)
    return {'NACC': NACC_selected, 'AIBL': AIBL_selected}

def construct_endpoint(source, cohort, endpoint):
    """Use every harmonized diagnosis row; select measurement rows independently.

    MCI: first recorded MCI; dementia without prior MCI is an explicit competing
    stop (event=0). Resulting single-event predictions are NET MCI risk, not CIF.
    """
    if endpoint not in {'AD', 'MCI'}:
        raise ValueError('Endpoint must be AD or MCI')
    source = source.copy()
    source['ID'] = source.ID.astype(str)
    source = source.sort_values(['ID', 'VISIT', 'age'], kind='stable')
    intervals, audits = [], []
    for pid, full in source.groupby('ID', sort=False):
        first_ad = full.loc[full.Classification.eq('AD'), 'VISIT'].min()
        first_mci = full.loc[full.Classification.eq('MCI'), 'VISIT'].min()
        target = first_ad if endpoint == 'AD' else first_mci
        competing = first_ad if endpoint == 'MCI' else np.nan
        target_observed = pd.notna(target) and (pd.isna(competing) or target < competing)
        terminal = float(target if target_observed else
                         competing if pd.notna(competing) else full.VISIT.max())
        stop_type = 'target_event' if target_observed else 'dementia_competing_stop' if pd.notna(competing) else 'last_diagnosis'
        # Preserve the user's row-completeness policy, without using it for outcomes.
        candidates = full.loc[:, ~full.columns.str.startswith('NPIQ')].copy()
        threshold = 0.0 if cohort == 'NACC' else 0.3
        candidates = candidates.loc[candidates.isna().mean(axis=1) <= threshold]
        candidates = candidates.assign(_missing=candidates[PREDICTORS].isna().sum(axis=1))
        candidates = (candidates.sort_values(['VISIT', '_missing', 'age'], ascending=[True, True, False], kind='stable')
                      .drop_duplicates('VISIT').drop(columns='_missing'))
        baseline = float(candidates.VISIT.min()) if len(candidates) else np.nan
        a = dict(ID=pid, cohort=cohort, endpoint=endpoint, source_visits=len(full),
                 baseline_source_month=baseline, first_target_source_month=target,
                 first_dementia_source_month=first_ad, source_last_diagnosis_month=float(full.VISIT.max()),
                 terminal_source_month=terminal, terminal_type=stop_type, target_observed=int(target_observed),
                 status='included', n_pre48_visits=0, landmark_eligible=False)
        if not len(candidates):
            a['status'] = 'no_usable_predictor_visit'
        elif terminal <= baseline or (pd.notna(target) and target <= baseline):
            a['status'] = 'endpoint_or_competing_stop_by_baseline'
        elif endpoint == 'MCI' and candidates.iloc[0].Classification != 'HC':
            a['status'] = 'MCI_requires_HC_baseline'
        else:
            # Positive intervals only. A terminal/event row never supplies predictors.
            g = candidates[candidates.VISIT < terminal].copy()
            if g.empty:
                a['status'] = 'no_pre_endpoint_measurement'
            else:
                g['start'] = g.VISIT.to_numpy(float)-baseline
                g['stop'] = np.r_[g.VISIT.to_numpy(float)[1:], terminal]-baseline
                g['event'] = 0
                g.iloc[-1, g.columns.get_loc('event')] = int(target_observed)
                g['__competing_stop'] = 0
                if stop_type == 'dementia_competing_stop':
                    g.iloc[-1, g.columns.get_loc('__competing_stop')] = 1
                assert (g.stop > g.start).all() and g.start.iloc[0] == 0
                a['n_pre48_visits'] = int(g.start.le(48).sum())
                a['landmark_eligible'] = bool(terminal-baseline > 48 and a['n_pre48_visits'] >= 3)
                intervals.append(g[['ID', 'Classification']+PREDICTORS+['start', 'stop', 'event', '__competing_stop']])
        a['endpoint_month'] = terminal-baseline if pd.notna(baseline) else np.nan
        a['target_month'] = target-baseline if pd.notna(baseline) and pd.notna(target) else np.nan
        a['dementia_month'] = first_ad-baseline if pd.notna(baseline) and pd.notna(first_ad) else np.nan
        audits.append(a)
    longitudinal = pd.concat(intervals, ignore_index=True)
    audit = pd.DataFrame(audits)
    audit_eligible = audit[audit.landmark_eligible]
    assert not audit_eligible.target_month.le(48).any()
    if endpoint == 'MCI':
        assert not audit_eligible.dementia_month.le(48).any()
    return longitudinal, audit


EXPECTED_PREPROCESSING_CODE = {'preprocessing(AIBL+NACC) simplify.ipynb': 'acd3ba43042869712249077dc4a7a21487128e11b445e1250006ed5810843953', 'preprocessing(AIBL+NACC) simplify MCI.ipynb': 'fcce43a24a8a3df28bfeaa300c99cca1c74e2c3302da90516ac168c32c721202'}

def preprocessing_notebook(project, name):
    """The published preprocessing notebooks are kept for reference in code/ or code/not used/."""
    for folder in [project/'code', project/'code'/'not used']:
        if (folder/name).exists():
            return folder/name
    raise FileNotFoundError(f'{name} not found in code/ or code/not used/: it defines the harmonisation snapshot.')

def check_preprocessing_snapshot(project):
    for name, expected in EXPECTED_PREPROCESSING_CODE.items():
        nb = json.loads(preprocessing_notebook(project, name).read_text())
        codes = [''.join(c['source']) for c in nb['cells'] if c['cell_type']=='code']
        actual = hashlib.sha256('\n'.join(codes).encode()).hexdigest()
        if actual != expected:
            raise ValueError(f'{name} code has changed since the harmonisation snapshot. Review and update prepare_history_information.py before reusing these inputs; existing CSVs cannot silently override refined preprocessing.')

def ensure_information_inputs(project, output_root=None, force=False):
    """Rebuild only when raw inputs, harmonisation snapshot, or this code change.

    A shared source cache avoids rereading the wide AIBL spreadsheet for each endpoint.
    CSVs remain local; the original preprocessed_data directory is never overwritten.
    """
    project = Path(project).expanduser().resolve()
    check_preprocessing_snapshot(project)
    output_root = Path(output_root or os.environ.get('SURVIVAL_INFORMATION_DATA_DIR',
                       str(project/'results'/'preprocessing')))
    output_root.mkdir(parents=True, exist_ok=True)
    raw_dir = Path(os.environ.get('SURVIVAL_RAW_DIR', str(project/'raw_data')))
    paths = {n: raw_dir/n for n in ['naming_convension.xlsx', 'AIBL_2023.xlsx', 'investigator_ftldlbd_nacc66.csv']}
    fingerprints = {n: sha256(p) for n, p in paths.items()}
    fingerprints['prepare_history_information.py'] = sha256(__file__)
    provenance = {n: sha256(preprocessing_notebook(project, n)) for n in
                  ['preprocessing(AIBL+NACC) simplify.ipynb', 'preprocessing(AIBL+NACC) simplify MCI.ipynb']}
    required = [output_root/f'{ep}'/f'{prefix}_{co}_{ep}.csv' for ep in ['AD', 'MCI']
                for co in ['NACC', 'AIBL'] for prefix in ['longitudinal', 'outcome_audit']]
    manifest_path = output_root/'input_manifest.json'
    if not force and manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest.get('fingerprints') == fingerprints and manifest.get('preprocessing_notebooks') == provenance:
            if all(p.exists() and sha256(p) == manifest.get('output_sha256', {}).get(str(p.relative_to(output_root))) for p in required):
                return output_root
    print('Rebuilding outcome-safe inputs from raw source files; AIBL loading may take several minutes.', flush=True)
    sources = harmonize_sources(project, output_root)
    summary = []
    for ep in ['AD', 'MCI']:
        ep_dir = output_root/ep
        ep_dir.mkdir(exist_ok=True)
        for co, source in sources.items():
            longitudinal, audit = construct_endpoint(source, co, ep)
            longitudinal.to_csv(ep_dir/f'longitudinal_{co}_{ep}.csv', index=False)
            audit.to_csv(ep_dir/f'outcome_audit_{co}_{ep}.csv', index=False)
            eligible = audit[audit.landmark_eligible]
            summary.append(dict(cohort=co, endpoint=ep, source_participants=len(audit),
                                usable_baseline_participants=longitudinal.ID.nunique(),
                                landmark_eligible=len(eligible),
                                events_after_48_by_108=int((eligible.target_observed.eq(1)&eligible.endpoint_month.le(108)).sum()),
                                competing_stops_after_48=int(eligible.terminal_type.eq('dementia_competing_stop').sum())))
            print(summary[-1], flush=True)
    pd.DataFrame(summary).to_csv(output_root/'cohort_summary.csv', index=False)
    manifest = dict(schema_version=SCHEMA_VERSION, fingerprints=fingerprints, preprocessing_notebooks=provenance,
                    time_origin='first usable predictor visit; months retain user age-derived integer-month convention',
                    outcomes='all harmonized diagnosis rows, before predictor missingness and duplicate selection',
                    visit_selection='no backward endpoint-dependent spacing or future visit-count filter',
                    mci_estimand='net first-MCI survival; dementia without prior MCI censored at its diagnosis; not cumulative incidence',
                    output_sha256={str(p.relative_to(output_root)): sha256(p) for p in required})
    manifest_path.write_text(json.dumps(manifest, indent=2))
    return output_root


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', default='/Users/yihan/Documents/Github/survival_Anlysis_Clustering')
    parser.add_argument('--output')
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args()
    print(ensure_information_inputs(args.project, args.output, args.force))
