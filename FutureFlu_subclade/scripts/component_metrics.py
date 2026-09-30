"""
Package-local component metric helpers for subclade risk scoring.

Self-contained helpers used by the 2025/2026 subclade pipeline.
"""
from __future__ import annotations

import pandas as pd
import numpy as np

def count_mutations(mutation_group):
    if pd.isna(mutation_group) or mutation_group == '':
        return 0
    return len(mutation_group.split(','))


def calculate_single_mutation_mi(mutation, seq_df, all_mutation_groups):
    site     = int(''.join(filter(str.isdigit, mutation)))
    col_name = f"X{site}"

    if col_name not in seq_df.columns:
        return 0

    has_mutation      = seq_df[seq_df[col_name] == mutation[-1]]
    total_occurrences = len(has_mutation)

    if total_occurrences == 0:
        return 0


    co_mutations = set()
    for group in all_mutation_groups:
        muts_in_group = [m.strip() for m in group.split(',')]
        if mutation in muts_in_group:
            for m in muts_in_group:
                if m != mutation:
                    co_mutations.add(m)

    # solo
    solo_seqs = has_mutation.copy()
    for co_mut in co_mutations:
        co_site = int(''.join(filter(str.isdigit, co_mut)))
        co_col  = f"X{co_site}"
        if co_col in seq_df.columns:
            solo_seqs = solo_seqs[solo_seqs[co_col] != co_mut[-1]]

    solo_occurrences = len(solo_seqs)
    return solo_occurrences / total_occurrences


def calculate_group_mutual_information(mutation_matrix):
    if mutation_matrix.shape[1] <= 1:
        return 0
    n = len(mutation_matrix)
    marginal_probs = mutation_matrix.mean()
    patterns = mutation_matrix.apply(lambda x: ''.join(x.astype(str)), axis=1)
    joint_counts = patterns.value_counts()
    MI = 0
    for pattern, count in joint_counts.items():
        p_joint = count / n
        if p_joint > 0:
            binary_pattern = [int(b) for b in pattern]
            p_indep = 1
            for i, mut in enumerate(mutation_matrix.columns):
                p_i = (marginal_probs[mut] if binary_pattern[i] == 1
                       else 1 - marginal_probs[mut])
                p_indep *= p_i
            if p_indep > 0:
                MI += p_joint * np.log2(p_joint / p_indep)
    return MI / mutation_matrix.shape[1]


def get_mutation_matrix_simple(seq_df, mutations):
    mutation_matrix = pd.DataFrame(index=seq_df.index)
    for mut in mutations:
        site = int(''.join(filter(str.isdigit, mut)))
        col_name = f"X{site}"
        if col_name in seq_df.columns:
            mutation_matrix[mut] = (seq_df[col_name] == mut[-1]).astype(int)
        else:
            mutation_matrix[mut] = 0
    return mutation_matrix


def get_matching_sequences(mutation_group, sequence_df, all_mutation_groups):
    if pd.isna(mutation_group):
        return pd.DataFrame()
    mutations = [m.strip() for m in mutation_group.split(',')]
    matching = sequence_df.copy()
    for mut in mutations:
        site = int(''.join(filter(str.isdigit, mut)))
        col_name = f"X{site}"
        if col_name in sequence_df.columns:
            matching = matching.loc[matching[col_name] == mut[-1]]
    for other_group in all_mutation_groups:
        other_muts = [m.strip() for m in other_group.split(',')]
        if set(other_muts) > set(mutations):
            for extra_mut in set(other_muts) - set(mutations):
                site = int(''.join(filter(str.isdigit, extra_mut)))
                col_name = f"X{site}"
                if col_name in sequence_df.columns:
                    matching = matching.loc[matching[col_name] != extra_mut[-1]]
    return matching


def get_clade_info_from_matching(matching):
    hidden = {"", "unknown", "unassigned", "nan", "none"}
    if matching is None or len(matching) == 0:
        return "unknown"
    col = "subclade" if "subclade" in matching.columns else "clade"
    if col not in matching.columns:
        return "unknown"
    clean = matching[col].fillna("").astype(str).str.strip()
    clean = clean[~clean.str.lower().isin(hidden)]
    if clean.empty:
        return "unknown"
    counts = clean.value_counts()
    total = int(counts.sum())
    return ", ".join(f"{name} ({count / total * 100:.1f}%)" for name, count in counts.items())


def calculate_antigenic_novelty_from_matching(matching, antigenic_novelty_dict):
    if len(matching) == 0:
        return 0
    return sum(antigenic_novelty_dict.get(acc, 0) for acc in matching['accession_number'])


def get_clade_info_with_percentages(mutation_group, sequence_df, all_mutation_groups):
    matching = get_matching_sequences(mutation_group, sequence_df, all_mutation_groups)
    return get_clade_info_from_matching(matching)


def calculate_total_escape_value(mutation_group, mutation_escape, site_escape,
                                  virus_type, sites_df, mutations_df):
    if pd.isna(mutation_group):
        return 0
    mutations = [m.strip() for m in mutation_group.split(',')]
    total = 0
    need_adj = (virus_type == 'Victoria' and
                len(sites_df) != 585 and len(mutations_df) != 585)
    for mut in mutations:
        site = int(''.join(filter(str.isdigit, mut)))
        if need_adj and site >= 177:
            adj_site = site - (585 - len(sites_df))
            adj_mut = f"{adj_site}{mut[-1]}"
            if adj_mut in mutation_escape:
                total += mutation_escape[adj_mut]
            elif str(adj_site) in site_escape:
                total += site_escape[str(adj_site)]
        else:
            if mut in mutation_escape:
                total += mutation_escape[mut]
            elif str(site) in site_escape:
                total += site_escape[str(site)]
    return total


def calculate_metric_value(mutation_group, mutation_metric, site_metric,
                            virus_type, sites_df, mutations_df):
    if pd.isna(mutation_group):
        return 0
    mutations = [m.strip() for m in mutation_group.split(',')]
    total = 0
    need_adj = (virus_type == 'Victoria' and
                len(sites_df) != 585 and len(mutations_df) != 585)
    for mut in mutations:
        site = int(''.join(filter(str.isdigit, mut)))
        if need_adj and site >= 177:
            adj_site = site - (585 - len(sites_df))
            adj_mut = f"{adj_site}{mut[-1]}"
            if adj_mut in mutation_metric:
                total += mutation_metric[adj_mut]
            elif str(adj_site) in site_metric:
                total += site_metric[str(adj_site)]
        else:
            if mut in mutation_metric:
                total += mutation_metric[mut]
            elif str(site) in site_metric:
                total += site_metric[str(site)]
    return total


def calculate_prevalence(mutation_group, mutation_prevalence):
    if pd.isna(mutation_group):
        return 0
    mutations = [m.strip() for m in mutation_group.split(',')]
    total = sum(mutation_prevalence.get(m, 0) for m in mutations)
    return total / len(mutations) if mutations else 0


def calculate_antigenic_novelty(mutation_group, sequence_df, all_mutation_groups,
                                 antigenic_novelty_dict):
    matching = get_matching_sequences(mutation_group, sequence_df, all_mutation_groups)
    return calculate_antigenic_novelty_from_matching(matching, antigenic_novelty_dict)


def get_dominant_clade(clade_str):

    if pd.isna(clade_str) or clade_str.strip() in ("Unknown", ""):
        return None
    first_entry = clade_str.split(',')[0].strip()
    clade_name  = first_entry.split('(')[0].strip()
    return clade_name if clade_name else None


def filter_random_single_mutations(results_df):
    if results_df.empty:
        return results_df

    results_df = results_df.copy()
    results_df['_dom_clade'] = results_df['clade'].apply(get_dominant_clade)


    multi_rows     = results_df[results_df['mutation_count'] > 1]
    multi_mc_pairs = {}  # {(mutation, dominant_clade): max_seq_count}
    for _, row_b in multi_rows.iterrows():
        dc_b       = row_b['_dom_clade']
        seq_count_b = row_b['mutation_group_seq_count']
        if dc_b is None:
            continue
        for m in [x.strip() for x in row_b['risk_mutation_group'].split(',')]:
            key = (m, dc_b)
            if key not in multi_mc_pairs or seq_count_b > multi_mc_pairs[key]:
                multi_mc_pairs[key] = seq_count_b


    def should_exclude(row):
        if row['mutation_count'] != 1:
            return False
        m_a         = row['risk_mutation_group'].strip()
        dc_a        = row['_dom_clade']
        seq_count_a = row['mutation_group_seq_count']
        if dc_a is None:
            return False
        key = (m_a, dc_a)
        if key not in multi_mc_pairs:
            return False

        max_seq_count_b = multi_mc_pairs[key]
        return max_seq_count_b >= seq_count_a * 0.1

    mask = results_df.apply(should_exclude, axis=1)
    excluded_count = mask.sum()
    if excluded_count > 0:
        print(f"  [filter] 排除了 {excluded_count} 个随机单突变组", flush=True)

    return results_df[~mask].drop(columns=['_dom_clade'])
