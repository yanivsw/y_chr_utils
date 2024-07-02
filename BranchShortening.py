import pandas as pd
import numpy as np
import utils.YChrDataset as YChrDataset

def calculate_ci(series, lower=0.025, upper=0.975):
    return series.quantile(lower), series.quantile(upper)

def get_branch_lengths(input_gt, ancient_non_african_ind, african_ind, non_african_ind, sites=None):
    gt_df = input_gt[['REF', 'ALT', 'chimp', ancient_non_african_ind, non_african_ind, african_ind]].dropna(subset=['chimp', ancient_non_african_ind, non_african_ind, african_ind])
    gt_df = gt_df[~gt_df['ALT'].str.contains(",")]

    if sites == 'transitions':
        gt_df = gt_df[( ( (gt_df['REF'] == 'C') & (gt_df['ALT'] == 'T') ) |
                        ( (gt_df['REF'] == 'G') & (gt_df['ALT'] == 'A') ) |
                        ( (gt_df['REF'] == 'T') & (gt_df['ALT'] == 'C') ) |
                        ( (gt_df['REF'] == 'A') & (gt_df['ALT'] == 'G') ) )]
    elif sites == 'transversions':
        gt_df = gt_df[~( ( (gt_df['REF'] == 'C') & (gt_df['ALT'] == 'T') ) |
                         ( (gt_df['REF'] == 'G') & (gt_df['ALT'] == 'A') ) |
                         ( (gt_df['REF'] == 'T') & (gt_df['ALT'] == 'C') ) |
                         ( (gt_df['REF'] == 'A') & (gt_df['ALT'] == 'G') ) )]

    # Derived non-African shared
    non_african_shared = sum( (gt_df['chimp'] == gt_df[african_ind]) &
                              (gt_df['chimp'] != gt_df[non_african_ind]) &
                              (gt_df[non_african_ind] == gt_df[ancient_non_african_ind]) )

    # Derived present-day non-African
    present_day_non_african = sum( (gt_df['chimp'] == gt_df[african_ind]) &
                                   (gt_df['chimp'] == gt_df[ancient_non_african_ind]) &
                                   (gt_df['chimp'] != gt_df[non_african_ind]) )

    # Derived ancient non-African
    ancient_non_african = sum( (gt_df['chimp'] == gt_df[african_ind]) &
                               (gt_df['chimp'] == gt_df[non_african_ind]) &
                               (gt_df['chimp'] != gt_df[ancient_non_african_ind]) )

    # Derived African
    present_day_african = sum( (gt_df['chimp'] == gt_df[non_african_ind]) &
                               (gt_df['chimp'] == gt_df[ancient_non_african_ind]) &
                               (gt_df['chimp'] != gt_df[african_ind]) )

    # Discordant branches (errors or recurrent mutations)
    # discordant1 = sum( (gt_df['chimp'] == gt_df[ancient_non_african]) &
    #                    (gt_df['chimp'] != gt_df[african]) &
    #                    (gt_df[african] == gt_df[non_african]) )

    # discordant2 = sum( (gt_df['chimp'] == gt_df[non_african]) &
    #                    (gt_df['chimp'] != gt_df[ancient_non_african]) &
    #                    (gt_df[african] == gt_df[ancient_non_african]) )

    total = len(gt_df)

    return {
        'non_african_shared'      : non_african_shared,
        'present_day_non_african' : present_day_non_african,
        'ancient_non_african'     : ancient_non_african,
        'present_day_african'     : present_day_african,
        'total'                   : total
    }

def calculate_tmrca_afr_and_mut_rate(branch_lengths, ancient_age):
    mut_rate = (branch_lengths['present_day_non_african'] - branch_lengths['ancient_non_african']) / ancient_age / branch_lengths['total']
    tmrca_mh_african_branch = branch_lengths['present_day_african'] / (branch_lengths['total'] * mut_rate)
    tmrca_mh_non_african_branch = (branch_lengths['non_african_shared'] + branch_lengths['present_day_non_african']) / (branch_lengths['total'] * mut_rate)
    return mut_rate, tmrca_mh_african_branch, tmrca_mh_non_african_branch

def simulate_tmrca_afr_and_mut_rate_ci(branch_lengths, ancient_age, rng, n=1000):
    simulated = pd.DataFrame({key: rng.poisson(val, 1000) for key, val in branch_lengths.items()})

    simulated['mut_rate'] = (simulated['present_day_non_african'] - simulated['ancient_non_african']) / ancient_age / branch_lengths['total']
    simulated['tmrca_mh_african_branch'] = simulated['present_day_african'] / (branch_lengths['total'] * simulated['mut_rate'])
    simulated['tmrca_mh_non_african_branch'] = (simulated['non_african_shared'] + simulated['present_day_non_african']) / (branch_lengths['total'] * simulated['mut_rate'])
    return simulated

def get_tmrca_afr_and_mut_rate(input_gt, ancient_non_african, ancient_age, african, non_african_samples, sites=None):
    rng = np.random.default_rng()

    df_ci = pd.DataFrame()

    mut_rate_res = {}
    tmrca_mh_african_branch_res = {}
    tmrca_mh_non_african_branch_res = {}
    n_positions = {}
    for non_african in non_african_samples:
        branch_lengths = get_branch_lengths(input_gt, ancient_non_african, african, non_african, sites)

        mut_rate, tmrca_mh_african_branch, tmrca_mh_non_african_branch = calculate_tmrca_afr_and_mut_rate(branch_lengths, ancient_age)

        temp_ci = simulate_tmrca_afr_and_mut_rate_ci(branch_lengths, ancient_age, rng)
        temp_ci['non_African'] = non_african

        df_ci = pd.concat([df_ci, temp_ci], ignore_index=True)

        n_positions[non_african] = branch_lengths['total']

        mut_rate_res[non_african] = (mut_rate, *calculate_ci(temp_ci['mut_rate']))
        tmrca_mh_african_branch_res[non_african] = (tmrca_mh_african_branch, *calculate_ci(temp_ci['tmrca_mh_african_branch']))
        tmrca_mh_non_african_branch_res[non_african] = (tmrca_mh_non_african_branch, *calculate_ci(temp_ci['tmrca_mh_non_african_branch']))

    mut_rate_res['total'] = ((np.mean([x[0] for x in mut_rate_res.values()])), *calculate_ci(df_ci['mut_rate']))
    tmrca_mh_african_branch_res['total'] = ((np.mean([x[0] for x in tmrca_mh_african_branch_res.values()])), *calculate_ci(df_ci['tmrca_mh_african_branch']))
    tmrca_mh_non_african_branch_res['total'] = ((np.mean([x[0] for x in tmrca_mh_non_african_branch_res.values()])), *calculate_ci(df_ci['tmrca_mh_non_african_branch']))

    res = {
            'ancient_non_african' : ancient_non_african,
            'african' : african,
            'mut_rate' : mut_rate_res,
            'n_positions' : n_positions,
            'tmrca_mh_african_branch' : tmrca_mh_african_branch_res,
            'tmrca_mh_non_african_branch' : tmrca_mh_non_african_branch_res
        }

    return res

def get_archaic_branch_lengths(input_gt, african_ind, non_african_ind, archaic_ind, sites=None, remove_hets=False):
    gt_df = input_gt[['REF', 'ALT', 'chimp', archaic_ind, non_african_ind, african_ind]].dropna(subset=['chimp', archaic_ind, non_african_ind, african_ind])
    gt_df = gt_df[~gt_df['ALT'].str.contains(",")]

    if remove_hets:
        gt_df = gt_df.drop(gt_df[gt_df['Allele1'] != gt_df['Allele2']].index)

    if sites == 'transitions':
        gt_df = gt_df[( ( (gt_df['REF'] == 'C') & (gt_df['ALT'] == 'T') ) |
                        ( (gt_df['REF'] == 'G') & (gt_df['ALT'] == 'A') ) |
                        ( (gt_df['REF'] == 'T') & (gt_df['ALT'] == 'C') ) |
                        ( (gt_df['REF'] == 'A') & (gt_df['ALT'] == 'G') ) )]
    elif sites == 'transversions':
        gt_df = gt_df[~( ( (gt_df['REF'] == 'C') & (gt_df['ALT'] == 'T') ) |
                         ( (gt_df['REF'] == 'G') & (gt_df['ALT'] == 'A') ) |
                         ( (gt_df['REF'] == 'T') & (gt_df['ALT'] == 'C') ) |
                         ( (gt_df['REF'] == 'A') & (gt_df['ALT'] == 'G') ) )]

    # Derived MH shared
    modern_human_shared = sum( (gt_df['chimp'] == gt_df[archaic_ind]) &
                               (gt_df['chimp'] != gt_df[non_african_ind]) &
                               (gt_df[non_african_ind] == gt_df[african_ind]) )

    # Derived non-African
    non_african = sum( (gt_df['chimp'] == gt_df[archaic_ind]) &
                       (gt_df['chimp'] == gt_df[african_ind]) &
                       (gt_df['chimp'] != gt_df[non_african_ind]) )

    # Derived African
    african = sum( (gt_df['chimp'] == gt_df[archaic_ind]) &
                   (gt_df['chimp'] == gt_df[non_african_ind]) &
                   (gt_df['chimp'] != gt_df[african_ind]) )

    # Derived archaic
    archaic = sum( (gt_df['chimp'] == gt_df[non_african_ind]) &
                   (gt_df['chimp'] == gt_df[african_ind]) &
                   (gt_df['chimp'] != gt_df[archaic_ind]) )

    # # Discordant branches (denote errors or recurrent mutations)
    # dicordant1 = sum( (gt_df['chimp'] == gt_df[african]) &
    #                   (gt_df['chimp'] != gt_df[archaic]) &
    #                   (gt_df[archaic] == gt_df[non_african]) )

    # dicordant2 = sum( (gt_df['chimp'] == gt_df[non_african]) &
    #                   (gt_df['chimp'] != gt_df[african]) &
    #                   (gt_df[archaic] == gt_df[african]) )

    total = len(gt_df)

    return {
        'modern_human_shared' : modern_human_shared,
        'non_african'         : non_african,
        'african'             : african,
        'archaic'             : archaic,
        'total'               : total,
    }

def calculate_archaic_tmrca_and_age(branch_lengths, res_tmrca_afr, non_african):
    tmrca_arch_non_african_branch = ( ( branch_lengths['modern_human_shared'] + branch_lengths['non_african'] ) / branch_lengths['non_african'] ) * res_tmrca_afr['tmrca_mh_non_african_branch'][non_african][0]
    tmrca_arch_african_branch     = ( ( branch_lengths['modern_human_shared'] + branch_lengths['african'] ) / branch_lengths['african'] ) * res_tmrca_afr['tmrca_mh_african_branch'][non_african][0]

    # tmrca_arch_non_african_branch = ( branch_lengths['modern_human_shared'] + branch_lengths['non_african'] ) / ( res_tmrca_afr['mut_rate'][non_african][0] * branch_lengths['total'] )
    # tmrca_arch_african_branch     = ( branch_lengths['modern_human_shared'] + branch_lengths['african'] ) / ( res_tmrca_afr['mut_rate'][non_african][0] * branch_lengths['total'] )

    age_arch_non_african_branch = tmrca_arch_non_african_branch * ( ( branch_lengths['modern_human_shared'] + branch_lengths['non_african'] - branch_lengths['archaic'] ) / ( branch_lengths['modern_human_shared'] + branch_lengths['non_african'] ) )
    age_arch_african_branch     = tmrca_arch_african_branch * ( ( branch_lengths['modern_human_shared'] + branch_lengths['african'] - branch_lengths['archaic'] ) / ( branch_lengths['modern_human_shared'] + branch_lengths['african'] ) )

    return (
        tmrca_arch_non_african_branch,
        tmrca_arch_african_branch,
        age_arch_non_african_branch,
        age_arch_african_branch,
    )

def simulate_archaic_tmrca_and_age_ci(branch_lengths, res_tmrca_afr, non_african, rng):
    simulated = pd.DataFrame({key: rng.poisson(val, 1000) for key, val in branch_lengths.items()})

    simulated['tmrca_arch_non_african_branch'] = ( ( simulated['modern_human_shared'] + simulated['non_african'] ) / simulated['non_african'] ) * res_tmrca_afr['tmrca_mh_non_african_branch'][non_african][0]
    simulated['tmrca_arch_african_branch']     = ( ( simulated['modern_human_shared'] + simulated['african'] ) / simulated['african'] ) * res_tmrca_afr['tmrca_mh_african_branch'][non_african][0]

    simulated['age_arch_non_african_branch'] = simulated['tmrca_arch_non_african_branch'] * ( ( simulated['modern_human_shared'] + simulated['non_african'] - simulated['archaic'] ) / ( simulated['modern_human_shared'] + simulated['non_african'] ) )
    simulated['age_arch_african_branch']     = simulated['tmrca_arch_african_branch'] * ( ( simulated['modern_human_shared'] + simulated['african'] - simulated['archaic'] ) / ( simulated['modern_human_shared'] + simulated['african'] ) )
    return simulated

def get_archaic_tmrca_and_age(gt_df, archaic, african, non_african_samples, res_tmrca_afr, sites=None):
    rng = np.random.default_rng()

    tmrca_arch_non_african_branch_res = {}
    tmrca_arch_african_branch_res = {}
    age_arch_non_african_branch_res = {}
    age_arch_african_branch_res = {}
    total_sites = {}

    df_ci = pd.DataFrame()
    for non_african in non_african_samples:
        branch_lengths = get_archaic_branch_lengths(gt_df, african, non_african, archaic, sites)

        (
            tmrca_arch_non_african_branch,
            tmrca_arch_african_branch,
            age_arch_non_african_branch,
            age_arch_african_branch,
        ) = calculate_archaic_tmrca_and_age(branch_lengths, res_tmrca_afr, non_african)

        temp_ci = simulate_archaic_tmrca_and_age_ci(branch_lengths, res_tmrca_afr, non_african, rng)

        tmrca_arch_non_african_branch_res[non_african] = (tmrca_arch_non_african_branch, *calculate_ci(temp_ci['tmrca_arch_non_african_branch']))
        tmrca_arch_african_branch_res[non_african]     = (tmrca_arch_african_branch, *calculate_ci(temp_ci['tmrca_arch_african_branch']))

        age_arch_non_african_branch_res[non_african] = (age_arch_non_african_branch, *calculate_ci(temp_ci['age_arch_non_african_branch']))
        age_arch_african_branch_res[non_african]     = (age_arch_african_branch, *calculate_ci(temp_ci['age_arch_african_branch']))

        df_ci = pd.concat([df_ci, temp_ci], ignore_index=True)

        total_sites[non_african] = branch_lengths['total']

    tmrca_arch_non_african_branch_res['total'] = ((np.mean([x[0] for x in tmrca_arch_non_african_branch_res.values()])), *calculate_ci(df_ci['tmrca_arch_non_african_branch']))
    tmrca_arch_african_branch_res['total']     = ((np.mean([x[0] for x in tmrca_arch_african_branch_res.values()])), *calculate_ci(df_ci['tmrca_arch_african_branch']))

    age_arch_non_african_branch_res['total'] = ((np.mean([x[0] for x in age_arch_non_african_branch_res.values()])), *calculate_ci(df_ci['age_arch_non_african_branch']))
    age_arch_african_branch_res['total']     = ((np.mean([x[0] for x in age_arch_african_branch_res.values()])), *calculate_ci(df_ci['age_arch_african_branch']))

    total_sites['total'] = np.mean([x for x in total_sites.values()])

    res_tmrca = {
        'african' : african,
        'archaic' : archaic,
        'total_sites' : total_sites,
        'tmrca_arch_non_african_branch' : tmrca_arch_non_african_branch_res,
        'tmrca_arch_african_branch' : tmrca_arch_african_branch_res,
        'age_arch_non_african_branch' : age_arch_non_african_branch_res,
        'age_arch_african_branch' : age_arch_african_branch_res
    }

    return res_tmrca

def get_mutations_vs_individual(vcf_data, comp_individual, samples):
    rng = np.random.default_rng()

    res = {}
    for individual in samples:
        if individual not in vcf_data.columns:
            continue

        if individual == comp_individual:
            res[individual] = {
                'haplogroup'       : YChrDataset.GetHaplogroup(individual),
                'mutation_diff'    : [0, (0, 0)],
                'normalised_diff'  : [0, (0, 0)],
                'total'            : 0,
                'derived_test_ind' : (0, (0, 0)),
                'derived_comp_ind' : (0, (0, 0)),
            }
            continue

        genotype_df = vcf_data[['POS', 'REF', 'ALT', 'chimp', individual, comp_individual]].dropna(subset=['chimp', individual, comp_individual])
        genotype_df = genotype_df[~genotype_df['ALT'].str.contains(",")]

        derived_test_ind = genotype_df[(genotype_df['chimp'] == genotype_df[comp_individual]) &
                                       (genotype_df['chimp'] != genotype_df[individual])]
        derived_comp_ind = genotype_df[(genotype_df['chimp'] != genotype_df[comp_individual]) &
                                       (genotype_df['chimp'] == genotype_df[individual])]
        total = len(genotype_df)

        # derived_test_ind.value_counts(['REF', 'ALT'])
        # derived_comp_ind.value_counts(['REF', 'ALT'])

        num_derived_test_ind = len(derived_test_ind)
        num_derived_comp_ind = len(derived_comp_ind)

        mutation_diff = num_derived_test_ind - num_derived_comp_ind

        temp_ci = pd.DataFrame()
        temp_ci['derived_test_ind'] = rng.poisson(num_derived_test_ind, 1000)
        temp_ci['derived_comp_ind'] = rng.poisson(num_derived_comp_ind, 1000)
        temp_ci['mutation_diff'] = (temp_ci['derived_test_ind'] - temp_ci['derived_comp_ind'])

        res[individual] = {
            'haplogroup'       : YChrDataset.GetHaplogroup(individual),
            'mutation_diff'    : [mutation_diff, calculate_ci(temp_ci['mutation_diff'])],
            'normalised_diff'  : [mutation_diff/total,  tuple(x / total for x in calculate_ci(temp_ci['mutation_diff']))],
            'total'            : total,
            'derived_test_ind' : (num_derived_test_ind, calculate_ci(temp_ci['derived_test_ind'])),
            'derived_comp_ind' : (num_derived_comp_ind, calculate_ci(temp_ci['derived_comp_ind'])),
        }

    return res

def get_missing_mutations(input_gt, reference_individual, comparison_individual, ref_age, comp_age, expected_mutation_rate):
    rng = np.random.default_rng()

    if reference_individual == comparison_individual:
        return

    input_gt = input_gt[~input_gt['ALT'].str.contains(",")]

    genotype_df = input_gt[['POS', 'REF', 'ALT', 'chimp', comparison_individual, reference_individual]].dropna(subset=['chimp', comparison_individual, reference_individual])

    derived_ref = len( genotype_df[(genotype_df['chimp'] == genotype_df[comparison_individual]) &
                                   (genotype_df['chimp'] != genotype_df[reference_individual])] )
    derived_comp = len( genotype_df[(genotype_df['chimp'] != genotype_df[comparison_individual]) &
                                    (genotype_df['chimp'] == genotype_df[reference_individual])] )
    total = len(genotype_df)

    measured_mutation_diff = derived_ref - derived_comp
    expected_mutation_diff = round((comp_age - ref_age) * expected_mutation_rate * total)

    missing_mutations      = measured_mutation_diff - expected_mutation_diff
    proportion_missing     = ( missing_mutations / derived_ref )
    adjusted_mutation_rate = (1 - proportion_missing) * expected_mutation_rate

    temp_ci = pd.DataFrame()
    temp_ci['derived_ref']  = rng.poisson(derived_ref, 100)
    temp_ci['derived_comp'] = rng.poisson(derived_comp, 100)

    temp_ci['measured_mutation_diff'] = temp_ci['derived_ref'] - temp_ci['derived_comp']
    temp_ci['expected_mutation_diff'] = round((comp_age - ref_age) * expected_mutation_rate * total)
    temp_ci['missing_mutations']      = temp_ci['measured_mutation_diff'] - temp_ci['expected_mutation_diff']
    temp_ci['proportion_missing']     = (round(temp_ci['missing_mutations'] / temp_ci['derived_ref'], 5))
    temp_ci['adjusted_mutation_rate'] = (1 - temp_ci['proportion_missing']) * expected_mutation_rate

    return {
        'reference_ind' : [reference_individual, derived_ref],
        'comp_ind'      : [comparison_individual, derived_comp],
        'total'         : total,

        'measured_mutation_diff' : measured_mutation_diff,
        'expected_mutation_diff' : expected_mutation_diff,

        'missing_mutations'      : (missing_mutations, calculate_ci(temp_ci['missing_mutations'])),
        'prop_missing'           : (proportion_missing, calculate_ci(temp_ci['proportion_missing'])),
        'adjusted_mutation_rate' : (adjusted_mutation_rate, calculate_ci(temp_ci['adjusted_mutation_rate'])),

        # 'expected_comp_mutations' : expected_comp_mutations,
    }