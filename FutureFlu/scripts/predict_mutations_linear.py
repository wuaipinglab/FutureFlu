# -*- coding: utf-8 -*-
"""
Linear mutation-prevalence prediction for FutureFlu HA1 seasons.

Fit historical HA1 site prevalences and score candidate mutations for
a target season / hemisphere / subtype.
"""

import pandas as pd
import numpy as np
from datetime import datetime
from itertools import combinations
import re
import os


def site_prevalence(seq, predict_season, semisphere, subtype):


    ha1_ranges = {
        'H3N2': (17, 345),
        'H1N1': (18, 344),
        'Victoria': (16, 362)
    }
    ha1_range = ha1_ranges[subtype]

    years = sorted([y for y in seq['season'].unique() if y > 2009 and y < predict_season])

    AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']


    seq = seq.copy()
    seq['submission_date'] = seq['submission_date'].fillna(seq['collection_date'])
    seq['collection_date'] = pd.to_datetime(seq['collection_date'])
    seq['submission_date'] = pd.to_datetime(seq['submission_date'])


    site_columns = seq.filter(regex='^X\d+').columns
    site_columns = [col for col in site_columns
                   if ha1_range[0] <= int(re.findall(r'\d+', col)[0]) <= ha1_range[1]]

    for col in site_columns:
        seq[col] = seq[col].str.upper()

    annual_data = {}
    prevalence_data = []

    def get_consensus_sequence(data, site_columns):

        consensus = {}
        for col in site_columns:
            valid_data = data[col][data[col] != 'X']
            if len(valid_data) > 0:
                consensus[col] = valid_data.mode().iloc[0]
        return consensus

    def calculate_significance(new_aa_current, new_aa_prev, total_current, total_prev):

        from scipy import stats
        import numpy as np


        observed = np.array([
            [new_aa_current, total_current - new_aa_current],  # current season
            [new_aa_prev, total_prev - new_aa_prev]            # previous season
        ])


        if np.any(observed == 0):
            return 0.0

        _, p_value = stats.chi2_contingency(observed)[:2]
        return p_value

    for year in years:

        if semisphere == "North":
            start = f"{year}-09-01"
            end = f"{year+1}-02-01"
            submission_end = f"{predict_season}-02-01"

        else:
            start = f"{year}-02-01"
            end = f"{year}-09-01"
            submission_end = f"{predict_season-1}-09-01"


        season_data = seq[
            (seq['collection_date'] >= start) &
            (seq['collection_date'] < end) &
            (seq['submission_date'] < submission_end)
        ]

        # ===
        year_freq_data = {}
        year_counts_data = {}
        for site_col in site_columns:
            valid_data = season_data[site_col][season_data[site_col] != 'X']
            freq = valid_data.value_counts(normalize=True)
            counts = valid_data.value_counts()
            for aa in AMINO_ACIDS:
                col_name = f"{site_col}{aa}"
                year_freq_data[col_name] = freq.get(aa, 0.0)
            year_counts_data[site_col] = counts

        prevalence_data.append(pd.Series(year_freq_data, name=year))

        # ===
        annual_data[year] = {
            'freqs': year_freq_data,
            'counts': year_counts_data,
            'dominant_mutations': []
        }

        # ===
        if year != years[0]:

            if semisphere == "North":
                prev_start = f"{year-1}-09-01"
                prev_end = f"{year}-02-01"
            else:
                prev_start = f"{year-1}-02-01"
                prev_end = f"{year-1}-09-01"

            prev_data = seq[
                (seq['collection_date'] >= prev_start) &
                (seq['collection_date'] < prev_end) &
                (seq['submission_date'] < submission_end)
            ]


            current_consensus = get_consensus_sequence(season_data, site_columns)
            prev_consensus = get_consensus_sequence(prev_data, site_columns)


            prev_counts_data = {}
            for site_col in site_columns:
                valid_data = prev_data[site_col][prev_data[site_col] != 'X']
                prev_counts_data[site_col] = valid_data.value_counts()

            identified_mutations = []
            for site_col in site_columns:
                if (site_col in current_consensus and site_col in prev_consensus and
                    current_consensus[site_col] != prev_consensus[site_col]):

                    old_aa = prev_consensus[site_col]
                    new_aa = current_consensus[site_col]


                    current_counts = year_counts_data[site_col]
                    prev_counts = prev_counts_data[site_col]

                    new_aa_current = current_counts.get(new_aa, 0)
                    new_aa_prev = prev_counts.get(new_aa, 0)
                    total_current = sum(current_counts)
                    total_prev = sum(prev_counts)


                    p_value = calculate_significance(new_aa_current, new_aa_prev,
                                                   total_current, total_prev)
                    # if subtype == 'H3N2' or subtype == 'Victoria':
                    if p_value < 0.05:
                        mut = f"{site_col[1:]}{new_aa}"
                        identified_mutations.append(mut)

            annual_data[year]['dominant_mutations'] = sorted(identified_mutations)

    # ===
    prevalence_result = pd.concat(prevalence_data, axis=1).T


    prevalence_result['dominant_mutation'] = [
        ', '.join(annual_data[y]['dominant_mutations']) for y in years
    ]


    mutation_columns = [col for col in prevalence_result.columns
                       if col not in ['dominant_mutation', 'dominant_clades']]
    column_order = mutation_columns + ['dominant_mutation']

    prevalence_result = prevalence_result[column_order]

    return prevalence_result.reset_index().rename(columns={'index': 'season'})

########################################
# 2.
########################################
def gmeasure(prev_data, theta_range):
    """
        prev_data (DataFrame): site_prevalence
        theta_range (list): theta

        DataFrame
    """
    # ======================================================

    # ======================================================

    valid_years = prev_data['season'].tolist()
    year_theta_gsum = {year: {theta: 0.0 for theta in theta_range} for year in valid_years}

    # ======================================================

    # ======================================================
    for theta in theta_range:
        for col in prev_data.columns:

            if col.startswith('X'):

                values = prev_data[col].fillna(0).values
                n_years = len(values)


                mut = np.zeros(n_years, dtype=int)
                start = 0


                for r in range(n_years):


                    if values[r] >= theta and np.any(values[start:r] < theta):
                        low_pos = np.where(values[:r] < theta)[0]
                        if low_pos.size > 0:
                            a = low_pos[-1]

                            mut[a+1:r+1] = 1

                            start = r + 1


                yearly_gsum = values * mut
                for idx, year in enumerate(valid_years):
                    year_theta_gsum[year][theta] += yearly_gsum[idx]

    # ======================================================

    # ======================================================

    gsum_df = pd.DataFrame.from_dict(year_theta_gsum, orient='index')
    gsum_df.columns = [f'theta={theta:.2f}' for theta in theta_range]
    gsum_df = gsum_df.reset_index().rename(columns={'index': 'season'})

    return gsum_df


def setup_logging(base_path):
    """
    Args:
        base_path
    Returns:
    """
    # os.path.join
    log_file = os.path.join(base_path, "best_result_para.txt")


    def log_print(*args, **kwargs):
        # log_print
        with open(log_file, 'a', encoding='utf-8') as f:
            # *args
            print(*args, **kwargs, file=f)


    with open(log_file, 'w', encoding='utf-8') as f:
        f.write(f"分析开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("-" * 50 + "\n")

    return log_print


########################################
# 3. Fitness
########################################

def fit_regression(input_gmeasure, input_epi_data, log_print):
    """
    Args:
        input_gmeasure (DataFrame): gmeasure
        input_epi_data (DataFrame)
        log_print (function)

    Returns:
        tuple: (
    """

    import statsmodels.api as sm
    from statsmodels.genmod.families import Binomial

    # ===
    gmeasure = input_gmeasure.set_index('season').copy()
    epi_data = input_epi_data.copy()


    epi_data['Season'] = epi_data['Season'].astype(int)
    gmeasure.index = gmeasure.index.astype(str)

    # ===
    valid_gmeasure_years = [y for y in gmeasure.index if y.isdigit()]
    gmeasure_years = pd.to_numeric(valid_gmeasure_years, errors='coerce')
    gmeasure_years = gmeasure_years[~np.isnan(gmeasure_years)].astype(int).tolist()


    common_years = sorted(set(epi_data['Season']).intersection(gmeasure_years))


    best_overall = {
        'theta': None,
        'r_squared': -float('inf'),
        'years_used': None,
        'start_year': None
    }

    # ===
    max_year = max(common_years)
    log_print("\n开始不同年份范围的拟合分析:")


    for years_back in range(3, len(common_years) + 1):
        start_year = max_year - years_back + 1
        current_years = [y for y in common_years if start_year <= y <= max_year]

        if len(current_years) < 3:
            continue

        log_print(f"\n使用{years_back}年数据进行拟合 ({start_year}-{max_year}):")


        current_epi = epi_data[epi_data['Season'].isin(current_years)]
        current_seq = gmeasure.loc[[str(y) for y in current_years]]

        # ===
        results = []
        for col in current_seq.columns:
            if not col.startswith('theta='):
                continue

            theta = float(col.split('=')[1])


            X = current_seq[col].values.reshape(-1, 1)


            y_raw = current_epi['Positivity_Rate'].values
            y = y_raw / 100.0


            y = np.clip(y, 0.0001, 0.9999)


            X = sm.add_constant(X)

            try:

                model = sm.GLM(y, X, family=Binomial(link=sm.families.links.logit()))
                result = model.fit(disp=0)


                null_model = sm.GLM(y, np.ones(len(y)).reshape(-1, 1), family=Binomial(link=sm.families.links.logit()))
                null_result = null_model.fit(disp=0)
                r_squared = 1 - (result.deviance / null_result.deviance)

                results.append({
                    'Theta': theta,
                    'R_squared': round(r_squared, 3)
                })

            except Exception as e:
                log_print(f"拟合theta={theta}时出错: {str(e)}")
                continue

        # ===
        if results:
            results_df = pd.DataFrame(results)
            results_df = results_df.sort_values(['R_squared', 'Theta'], ascending=[False, True])
            best_current = results_df.iloc[0]

            log_print(f"当前最佳 - Theta: {best_current['Theta']:.2f}, R²: {best_current['R_squared']:.3f}")


            if best_current['R_squared'] > best_overall['r_squared']:
                best_overall['theta'] = best_current['Theta']
                best_overall['r_squared'] = best_current['R_squared']
                best_overall['years_used'] = years_back
                best_overall['start_year'] = start_year

    # ===
    if best_overall['theta'] is not None:
        log_print("\n=== 最终最优拟合结果 ===")
        log_print(f"使用年份范围: {best_overall['start_year']}-{max_year} ({best_overall['years_used']}年)")
        log_print(f"最优Theta: {best_overall['theta']:.2f}")
        log_print(f"最优R²: {best_overall['r_squared']:.3f}")
        return best_overall['theta'], best_overall['r_squared'], best_overall['years_used']
    else:
        log_print("\n警告: 未找到有效的拟合结果")
        return None, None, None
########################################
# 4.
########################################

def predict_mutations_multi_model(predict_season, theta, prev_data, mutation_history, output_file=None):
    """
        predict_season (int)
        theta (float)
        prev_data (DataFrame): site_prevalence
        mutation_history (DataFrame): risk_mutations
        output_file (str, optional)

        DataFrame
    """
    import os
    import pandas as pd
    import numpy as np


    models = ['linear']


    all_predictions = pd.DataFrame(columns=['predict_season', 'risk_mutation', 'previous_prevalence', 'predicted_prevalence', 'delta', 'model'])

    # ===

    mutation_columns = [col for col in prev_data.columns if col.startswith('X')]
    historical_data = prev_data[prev_data.index < predict_season][mutation_columns]

    # ===

    mutation_dominant_years = {}
    if 'dominant_mutation' in prev_data.columns:
        for idx, row in prev_data.iterrows():
            year = row.name if hasattr(row, 'name') else prev_data.iloc[idx]['season']
            if year >= predict_season:
                continue

            dominant_muts = row['dominant_mutation']
            if isinstance(dominant_muts, str) and dominant_muts.strip():
                for mut in dominant_muts.split(','):
                    mut = mut.strip()
                    if mut:

                        if mut not in mutation_dominant_years or year > mutation_dominant_years[mut]:
                            mutation_dominant_years[mut] = year


    dominant_mutations = set(mutation_dominant_years.keys())

    # ===
    for model_type in models:
        print(f"使用{model_type}模型进行预测...")


        risk_muts = []


        pred_prev = np.zeros(len(historical_data.columns))
        deltas = np.zeros(len(historical_data.columns))

        # ===
        for col_idx, col in enumerate(historical_data.columns):

            freqs = historical_data[col].values


            if len(freqs) < 2:
                pred_prev[col_idx] = freqs[-1] if len(freqs) > 0 else 0
                deltas[col_idx] = 0
                continue

            try:

                last_value = freqs[-1]
                prev_value = freqs[-2]


                if model_type == 'linear':

                    delta = last_value - prev_value
                    predicted = last_value + delta
                    deltas[col_idx] = delta

                elif model_type == 'exponential':


                    if prev_value > 0.0001 and last_value > 0.0001:
                        growth_rate = np.log(last_value / prev_value)
                        predicted = last_value * np.exp(growth_rate)
                    else:

                        delta = last_value - prev_value
                        predicted = last_value + delta

                elif model_type == 'logistic':
                    # Logistic

                    if 0.0001 < prev_value < 0.9999 and 0.0001 < last_value < 0.9999:

                        f_i = np.log((last_value / (1 - last_value)) / (prev_value / (1 - prev_value)))

                        exp_term = np.exp(f_i)
                        predicted = (last_value * exp_term) / (last_value * exp_term + (1 - last_value))
                    else:

                        delta = last_value - prev_value
                        predicted = last_value + delta

                elif model_type == 'gompertz':

                    if len(freqs) >= 3 and all(f > 0.0001 for f in freqs[-3:]):
                        try:

                            x0, x1, x2 = freqs[-3], freqs[-2], freqs[-1]


                            r1 = np.log(x1/x0)
                            r2 = np.log(x2/x1)


                            if abs(r1) > 0.0001:
                                decay = r2/r1


                                r3 = r2 * decay


                                predicted = x2 * np.exp(r3)
                            else:

                                predicted = x2 * np.exp(r2)
                        except:

                            growth_rate = np.log(last_value / prev_value)
                            predicted = last_value * np.exp(growth_rate)
                    else:

                        if prev_value > 0.0001 and last_value > 0.0001:
                            growth_rate = np.log(last_value / prev_value)
                            predicted = last_value * np.exp(growth_rate)
                        else:
                            delta = last_value - prev_value
                            predicted = last_value + delta

                pred_prev[col_idx] = np.clip(predicted, 0, 1)
            except Exception as e:
                print(f"预测{col}时出错({model_type}模型): {e}")

                pred_prev[col_idx] = freqs[-1] if len(freqs) > 0 else 0
                deltas[col_idx] = 0

        # ===
        for col_idx, col in enumerate(historical_data.columns):
            pred = pred_prev[col_idx]


            formatted_mut = col[1:] if col.startswith('X') else col

            # ===

            if formatted_mut in mutation_dominant_years:

                start_year = mutation_dominant_years[formatted_mut]

                historical_years = prev_data[prev_data.index < predict_season].index.tolist()
                try:
                    start_idx = historical_years.index(start_year)
                    freqs = historical_data[col].values[start_idx:]
                except ValueError:

                    freqs = historical_data[col].values
            else:

                freqs = historical_data[col].values

            # ===

            condition1 = (pred >= theta) and np.any(freqs < theta)


            condition2 = False
            if theta/10 >= 0.01:
                condition2 = (pred >= theta/10) and np.any(freqs < theta/10)
            elif theta*10 < 1:
                condition2 = (pred >= theta*10) and np.any(freqs < theta*10)


            current_prev = historical_data[col].iloc[-1] if len(historical_data[col]) > 0 else 0
            condition3 = current_prev < 0.75

            if (condition1 or condition2) and condition3:
                risk_muts.append(col)

        # ===
        model_records = []
        for mut in risk_muts:

            formatted_mut = mut[1:] if mut.startswith('X') else mut


            prev_prev = historical_data[mut].iloc[-1] if len(historical_data[mut]) > 0 else 0


            mut_idx = list(historical_data.columns).index(mut)
            predicted_prev = pred_prev[mut_idx]
            delta_val = deltas[mut_idx]

            model_records.append({
                'predict_season': predict_season,
                'risk_mutation': formatted_mut,
                'previous_prevalence': round(prev_prev, 4),
                'predicted_prevalence': round(predicted_prev, 4),
                'delta': round(delta_val, 4) if model_type == 'linear' else None,
                'model': model_type
            })

        # ===
        if model_records:
            current_predictions = pd.DataFrame(model_records)
            all_predictions = pd.concat([all_predictions, current_predictions], ignore_index=True)

            print(f"{model_type}模型预测发现{len(model_records)}个风险突变")
        else:
            print(f"{model_type}模型未预测到任何风险突变")


    if not all_predictions.empty:
        all_predictions = all_predictions.sort_values(by=['risk_mutation', 'model'])


    if output_file is not None:

        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)


        all_predictions.to_csv(output_file, index=False)
        print(f"已将所有模型预测结果保存至: {output_file}")

    print("所有模型预测完成!")
    return all_predictions

def analyze_risk_mutations(sequences_df, mutations_df, predict_season, semisphere, model=None):
    """
    Args:
        sequences_df (DataFrame)
        mutations_df (DataFrame)
        predict_season (int)
        semisphere (str): 'North'
        model (str, optional)

    Returns:
        DataFrame
            - risk_mutation_group
            - count
            - model
    """

    if semisphere == "North":
        start = f"{predict_season-1}-09-01"
        end = f"{predict_season}-02-01"
        submission_end = f"{predict_season}-02-01"
    else:
        start = f"{predict_season-1}-02-01"
        end = f"{predict_season-1}-09-01"
        submission_end = f"{predict_season-1}-09-01"


    sequences_df['collection_date'] = pd.to_datetime(sequences_df['collection_date'])
    sequences_df['submission_date'] = pd.to_datetime(sequences_df['submission_date'])


    date_mask = (
        (sequences_df['collection_date'] >= start) &
        (sequences_df['collection_date'] < end) &
        (sequences_df['submission_date'] < submission_end)
    )
    filtered_sequences = sequences_df[date_mask]


    if model is not None:

        if 'model' in mutations_df.columns:
            mutations_df = mutations_df[mutations_df['model'] == model]
        else:
            print(f"警告: 突变数据中没有'model'列，将使用所有突变数据")


    if mutations_df.empty:
        return pd.DataFrame(columns=['risk_mutation_group', 'count', 'model'])


    def parse_mutation(mutation_str):
        position = int(''.join(filter(str.isdigit, mutation_str)))
        amino_acid = mutation_str.replace(str(position), '')
        return position, amino_acid


    result_dfs = []


    if model is not None:
        models_to_analyze = [model]
    elif 'model' in mutations_df.columns:
        models_to_analyze = mutations_df['model'].unique()
    else:

        models_to_analyze = ['undefined']


    for current_model in models_to_analyze:

        if current_model != 'undefined' and 'model' in mutations_df.columns:
            model_mutations = mutations_df[mutations_df['model'] == current_model]
        else:
            model_mutations = mutations_df.copy()


        risk_mutations = [parse_mutation(mut) for mut in model_mutations['risk_mutation']]

        if not risk_mutations:
            continue


        model_sequences = filtered_sequences.copy()


        def find_risk_mutations(row):
            mutations_found = []
            for pos, aa in risk_mutations:
                column_name = f'X{pos}'
                if column_name in row:
                    if aa == '-':
                        if pd.isna(row[column_name]) or row[column_name] == '-':
                            mutations_found.append(f"{pos}{aa}")
                    else:
                        if row[column_name] == aa:
                            mutations_found.append(f"{pos}{aa}")
            return ','.join(sorted(mutations_found)) if mutations_found else None


        model_sequences['risk_mutation_group'] = model_sequences.apply(find_risk_mutations, axis=1)


        model_sequences = model_sequences[model_sequences['risk_mutation_group'].notna()]


        if not model_sequences.empty:
            mutation_counts = model_sequences['risk_mutation_group'].value_counts().reset_index()
            mutation_counts.columns = ['risk_mutation_group', 'count']

            mutation_counts['model'] = current_model
            result_dfs.append(mutation_counts)


    if result_dfs:
        return pd.concat(result_dfs, ignore_index=True)
    else:

        return pd.DataFrame(columns=['risk_mutation_group', 'count', 'model'])
