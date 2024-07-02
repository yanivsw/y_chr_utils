import pandas as pd
import io
import numpy as np
import itertools

def read_bed_file(bed_path):
    import gzip

    open_func = gzip.open if '.gz' in bed_path else open
    with open_func(bed_path, 'rt') as f:
        bed_file = pd.read_csv(f, sep='\t', header=None)
    return bed_file

def get_filter_coordinates(filter_input, filter_name, chr_size):
    if isinstance(filter_input, str):
        filter = read_bed_file(filter_input)
    else:
        filter = filter_input

    filter_df = pd.DataFrame({filter_name: np.zeros(chr_size, dtype=int)}, index=np.arange(1, chr_size + 1))
    for x1, x2 in zip(filter[1], filter[2]):
        filter_df.loc[x1 + 1:x2 + 1, filter_name] = 1
    return filter_df

def merge_bed_ranges(bed_input):
    if isinstance(bed_input, str):
        input_bed_df = read_bed_file(bed_input)
    else:
        input_bed_df = bed_input

    input_bed_df = input_bed_df[input_bed_df[0] == 'Y']
    input_bed_df = input_bed_df.sort_values(by=1)
    input_bed_df["group"] = (input_bed_df[1] > input_bed_df[2].shift().cummax()).cumsum()

    merged = input_bed_df.groupby(['group', 0]).agg({1: min, 2: max})
    merged = merged.reset_index()[[0, 1, 2]]
    return merged

def write_bed(input_bed_df, output_loc):
    input_bed_df[[0, 1, 2]].to_csv(output_loc, sep='\t', index=None, header=None)

def count_bed_ranges(bed_input):
    if isinstance(bed_input, str):
        bed_file_data = read_bed_file(bed_input)
    else:
        bed_file_data = bed_input

    bed_file_data['range'] = bed_file_data[2] - bed_file_data[1]

    return bed_file_data

def read_vcf(path):
    import gzip
    import re

    open_func = gzip.open if '.gz' in path else open
    with open_func(path, 'rt') as f:
        lines = re.sub(r'^##.*\n?', '', f.read(), flags=re.MULTILINE)

    temp_df = pd.read_csv(io.StringIO(lines),
                          sep='\t',
                          dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                                 'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str}).rename(columns={'#CHROM': 'CHROM'})

    return temp_df

def filter_vcf(vcf_dataframe, bed_file, inverse=False):
    if isinstance(bed_file, str):
        bed_file_data = read_bed_file(bed_file)
    else:
        bed_file_data = bed_file

    filter_ranges = [list(range(x1 + 1, x2 + 1)) for x1, x2 in zip(bed_file_data[1], bed_file_data[2])]
    filter_list = list(itertools.chain.from_iterable(filter_ranges))

    if not inverse:
        vcf_out = vcf_dataframe[vcf_dataframe['POS'].isin(filter_list)]
    else:
        vcf_out = vcf_dataframe[~vcf_dataframe['POS'].isin(filter_list)]

    return vcf_out.reset_index(drop=True)

def convert_vcf_to_dataframe(path, vcf_samples=None, convert_hets=True):
    if isinstance(path, str):
        if '.pkl' in path:
            vcf_data = pd.read_pickle(path)
        else:
            vcf_data = read_vcf(path)

    vcf_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    if vcf_samples is None:
        vcf_samples = [col for col in vcf_data.columns.values if col not in vcf_cols]

    format_strings = vcf_data['FORMAT'].unique()
    format_str = max(format_strings, key=len)
    format_template = list(format_str.split(':'))

    out_df = pd.DataFrame()
    out_df[['CHROM', 'POS', 'REF', 'ALT']] = vcf_data[['CHROM', 'POS', 'REF', 'ALT']]
    # out_df[['ID']] = vcf_data[['ID']]
    # out_df[['FILTER']] = vcf_data[['FILTER']]

    out_df.columns = pd.MultiIndex.from_arrays([out_df.columns, ['']*len(out_df.columns)])

    for vcf_sample in vcf_samples:
        print(vcf_sample)

        # filter out rows with missing or incomplete genotype information
        genotype_mask = ~vcf_data[vcf_sample].isin(['.', './.:.:.:.:.:.:.:.', '.:.:.:.:.:.:.:.'])

        # Split genotype information into separate columns based on ':' delimiter
        temp_df = pd.DataFrame(index=vcf_data.index, columns=[vcf_sample, *format_template, 'Allele1', 'Allele2', 'Heterozygote'])

        temp_df[vcf_sample][genotype_mask] = vcf_data[vcf_sample][genotype_mask]
        temp_df.loc[genotype_mask, format_template] = temp_df[vcf_sample][genotype_mask].str.split(':', expand=True).values

        if ('0/0' in temp_df['GT'].values):
            # Split the 'GT' field into 'Allele1' and 'Allele2' and convert to numeric
            temp_df.loc[genotype_mask, ['Allele1', 'Allele2']] = temp_df['GT'][genotype_mask].str.split('/', expand=True).values
            temp_df['Allele1'] = pd.to_numeric(temp_df['Allele1'], errors='coerce')
            temp_df['Allele2'] = pd.to_numeric(temp_df['Allele2'], errors='coerce')

            # Check if the genotype is heterozygous
            temp_df.loc[genotype_mask, 'GT'] = temp_df[['Allele1', 'Allele2']][genotype_mask].max(axis=1).astype(object)
            temp_df.loc[genotype_mask, 'Heterozygote'] = 0

            heterozygote_mask = (temp_df['Allele1'] != temp_df['Allele2']) & genotype_mask
            temp_df.loc[heterozygote_mask, 'Heterozygote'] = 1
        else:
            temp_df['GT'] = pd.to_numeric(temp_df['GT'], errors='coerce')
            temp_df['Heterozygote'] = 0

        if 'DP' in format_template:
            if len(temp_df['DP'].unique()) > 2:
                temp_df['DP'] = pd.to_numeric(temp_df['DP'], errors='coerce')

        if 'GQ' in format_template:
            if len(temp_df['GQ'].unique()) > 2:
                temp_df['GQ'] = pd.to_numeric(temp_df['GQ'], errors='coerce')

        # Retain only the specified columns
        retain_cols = ['GT', 'DP', 'GQ', 'Allele1', 'Allele2', 'Heterozygote']
        temp_df = temp_df[[col for col in retain_cols if col in temp_df.columns]]

        # Create a MultiIndex for columns using the sample name and the column names
        temp_df.columns = pd.MultiIndex.from_product([[vcf_sample], temp_df.columns])

        out_df = pd.concat([out_df, temp_df], axis=1)

    # Convert heterozygous genotypes to homozygous based on the minimum posterior probability
    if ('PP' in format_str) & convert_hets:
        snp_ad_gt = ['AA', 'CC', 'GG', 'TT', 'AC', 'AG', 'AT', 'CG', 'CT', 'Gt']
        pd.options.mode.chained_assignment = None

        for vcf_sample in vcf_samples:
            het_mask = out_df[(vcf_sample, 'Heterozygote')] == 1

            if het_mask.any():
                print(vcf_sample)

                # Extract hets
                hets = pd.DataFrame(index = out_df[het_mask].index, columns = ['REF', 'ALT', vcf_sample, *format_template, *snp_ad_gt])
                hets[vcf_sample] = vcf_data[vcf_sample][het_mask]
                hets[format_template] = hets[vcf_sample].str.split(':', expand=True).values
                hets[snp_ad_gt] = hets['PP'].str.split(',', expand=True).values
                hets[['REF', 'ALT']] = vcf_data.loc[het_mask, ['REF', 'ALT']].values

                # Determine the minimum PP value and extract the corresponding base
                hets['min_PP'] = hets[['AA', 'CC', 'GG', 'TT']].astype(int).idxmin(axis=1).str[0].values

                # Set homozygous genotype to zero if min_PP matches REF
                hets['homozygous_GT'] = np.nan
                hets['homozygous_GT'][hets['min_PP'] == hets['REF']] = 0

                # Set homozygous genotype if min_PP does not match REF and is not NaN
                mask = (hets['min_PP'] != hets['REF']) & hets['min_PP'].notna()
                if mask.sum() > 0:
                    temp = hets[mask]
                    # Find the position of min_PP in ALT
                    temp['pos'] = temp.apply(lambda x: x['ALT'].find(x['min_PP']), axis=1)
                    # Update homozygous genotype based on the position of min_PP
                    temp['homozygous_GT'][temp['pos'] == 0] = 1
                    temp['homozygous_GT'][temp['pos'] == 2] = 2
                    temp['homozygous_GT'][temp['pos'] == -1] = 2
                    hets['homozygous_GT'][mask] = temp['homozygous_GT']

                if len(vcf_samples) > 1:
                    # Extract genotypes for heterozygous samples, excluding the current sample
                    # and get unique, sorted, non-null genotypes
                    temp_gt_df = get_genotypes_from_vcf_dataframe(out_df[het_mask])[vcf_samples]
                    temp_gt_df = temp_gt_df.drop(vcf_sample, axis=1)
                    hets['other_GTs'] = temp_gt_df.apply(lambda x: list(set(x)), axis = 1).apply(lambda v: sorted(pd.Series(v).dropna().values))

                    # For rows where the current sample is homozygous reference (0)
                    # and the only other genotype is also reference (0),
                    # set the ALT field to '.'
                    hets['remove_alt'] = hets.apply(lambda x: (x['homozygous_GT'] == 0) &
                                                              (x['homozygous_GT'] in x['other_GTs']) &
                                                              (len(x['other_GTs']) == 1), axis=1)
                    hets['ALT'][hets['remove_alt']] = '.'

                    # For rows where the current sample is homozygous reference (0) or homozygous for the first allele (1),
                    # and the ALT field contains multiple alleles,
                    # and no other sample has the second allele (2),
                    # simplify the ALT field to only include the first allele
                    hets['fix_alt_1'] = hets.apply(lambda x: ((x['homozygous_GT'] == 0) | (x['homozygous_GT'] == 1)) &
                                                             (',' in x['ALT']) &
                                                             (2 not in x['other_GTs']), axis=1)
                    hets['ALT'][hets['fix_alt_1']] = hets['ALT'].str.replace("(,).*","", regex=True)

                    # For rows where the current sample is homozygous for the second allele (2),
                    # and no other sample has this allele,
                    # and the minimum posterior probability base is not in ALT
                    # append the min_PP to ALT
                    hets['fix_alt_2'] = hets.apply(lambda x: (x['homozygous_GT'] == 2) &
                                                             (2 not in x['other_GTs']) &
                                                             (x['min_PP'] not in x['ALT']), axis=1)
                    hets['ALT'][hets['fix_alt_2']] = hets['ALT'] + ',' + hets['min_PP']
                else:
                    # For a single sample, identify rows where the sample is homozygous reference (0) and set the ALT field to '.'
                    hets['remove_alt'] = hets.apply(lambda x: (x['homozygous_GT'] == 0), axis=1)
                    hets['ALT'][hets['remove_alt']] = '.'

                out_df[vcf_sample, 'GT'][het_mask] = hets['homozygous_GT']
                out_df['ALT', ''][het_mask] = hets['ALT']

        # Get rows with multiple ALT alleles
        mask = out_df['ALT'].str.contains(",")

        # Extract genotypes and get unique, sorted genotypes
        temp_gt_df = pd.DataFrame()
        temp_gt_df['GT'] = get_genotypes_from_vcf_dataframe(out_df[mask])[vcf_samples].apply(lambda x: list(set(x)), axis = 1).apply(lambda v: sorted(pd.Series(v).dropna().values))

        # Identify rows that contain only the second allele but not the first allele
        fix_alt_mask = temp_gt_df[temp_gt_df.apply(lambda x: (2 in x['GT']) & (1 not in x['GT']), axis=1)]

        # Update genotypes and ALT field
        temp = out_df.loc[fix_alt_mask.index.values]
        mask = temp.loc[:,pd.IndexSlice[:, 'GT']] == 2
        temp.loc[:,pd.IndexSlice[:, 'GT']] = temp.where(~mask, 1.0).xs('GT', level=1, axis=1, drop_level=False).astype(object)
        temp['ALT'] = temp['ALT'].str.replace("^.*?,","", regex=True)

        out_df.loc[fix_alt_mask.index.values] = temp

    return out_df

def merge_vcfs(vcf1, vcf2):
    vcf1 = pd.read_pickle(vcf1) if '.pkl' in vcf1 else convert_vcf_to_dataframe(vcf1)
    vcf2 = pd.read_pickle(vcf2) if '.pkl' in vcf2 else convert_vcf_to_dataframe(vcf2)

    # Merge VCF dataframes on CHROM, POS, and REF columns
    merged_vcf_data = pd.merge(vcf1, vcf2, how='outer', on=['CHROM', 'POS', 'REF']).sort_values(by='POS')

    # Extract sample columns
    vcf_samples = [col for col in merged_vcf_data.columns if col not in ['CHROM', 'POS', 'REF', 'QUAL', 'FILTER', 'INFO', 'ID', 'ALT', 'ALT_x', 'ALT_y', 'FORMAT', 'FORMAT_x', 'FORMAT_y']]

    out_df = pd.DataFrame()
    out_df[['CHROM', 'POS', 'REF', 'ALT_x', 'ALT_y']] = merged_vcf_data[['CHROM', 'POS', 'REF', 'ALT_x', 'ALT_y']]

    # Replace '.' and NaN in ALT_x and ALT_y with empty strings
    out_df['ALT_x'].replace({'.': '', np.nan: ''}, inplace=True)
    out_df['ALT_y'].replace({'.': '', np.nan: ''}, inplace=True)

    # Initialize 'ALT' with ALT_x values
    out_df['ALT'] = out_df['ALT_x']

    # For rows where ALT_x is empty, use ALT_y values
    out_df.loc[(out_df['ALT_x'] == ''), 'ALT'] = out_df['ALT_y']

    # For rows where ALT_x and ALT_y are different and non-empty, concatenate them
    diff_non_empty = (out_df['ALT_x'] != out_df['ALT_y']) & (out_df['ALT_x'] != '') & (out_df['ALT_y'] != '')
    out_df.loc[diff_non_empty, 'ALT'] = out_df[diff_non_empty].apply(
        lambda x: x['ALT_x'] if x['ALT_y'] in x['ALT_x'] else (x['ALT_y'] if x['ALT_x'] in x['ALT_y'] else x['ALT_x'] + ',' + x['ALT_y']), axis=1)

    # Drop the now redundant ALT_x and ALT_y columns
    out_df.drop(['ALT_x', 'ALT_y'], axis=1, inplace=True)

    # Flatten MultiIndex if either input dataframe has MultiIndex columns
    if vcf1.columns.nlevels > 1 or vcf2.columns.nlevels > 1:
        out_df.columns = out_df.columns.get_level_values(0)

    out_df[vcf_samples] = merged_vcf_data[vcf_samples]
    out_df.sort_values(by=['POS'], inplace=True)

    drop_cols = [col for col in out_df.columns.values if 'ID' in col[0] or 'FILTER' in col[0]]
    out_df.drop(drop_cols, axis=1, inplace=True)

    return out_df

def convert_snpad_dataframe_to_vcf(input_df, out_file):
    format_str = 'GT:DP:A:C:G:T:PP:GQ'
    snpad_vcf_header = '''##fileformat=VCFv4.1
##FORMAT=<ID=A,Number=2,Type=Integer,Description="Number of A bases on forward and reverse strand">
##FORMAT=<ID=C,Number=2,Type=Integer,Description="Number of C bases on forward and reverse strand">
##FORMAT=<ID=G,Number=2,Type=Integer,Description="Number of G bases on forward and reverse strand">
##FORMAT=<ID=T,Number=2,Type=Integer,Description="Number of T bases on forward and reverse strand">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality (difference between lowest and second lowest PP value)">
##FORMAT=<ID=PP,Number=10,Type=Integer,Description="Normalized, Phred-scaled posteror probability for genotypes AA, CC, GG, TT, AC, AG, AT, CG, CT, GT, in this order">
'''

    sample_name = list(set(list(input_df.columns.levels[0])) - set(['CHROM', 'POS', 'REF', 'ALT']))[0]

    # Prepare the output DataFrame with required VCF columns
    out = input_df[['CHROM', 'POS', 'REF', 'ALT']].copy()
    out = out.rename(columns={'CHROM': '#CHROM'})
    out['ID'], out['QUAL'], out['FILTER'], out['INFO'], out['FORMAT'] = '.', input_df[sample_name, 'GQ'], '.', '.', format_str
    out[sample_name] = input_df[sample_name, 'snpAD']

    # Convert the DataFrame to a VCF formatted string and prepend the header
    out_str = snpad_vcf_header + out.to_csv(sep='\t', index=False)

    # Write the VCF string to a file
    with open(out_file, 'w+') as fh:
        fh.write(out_str)

def convert_gt_dataframe_to_vcf(input_df, out_file):
    format_str = 'GT'
    snpad_vcf_header = '''##fileformat=VCFv4.1
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
'''
    sample_name = list(set(list(input_df.columns.levels[0])) - set(['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'ID']))[0]

    out = input_df[['CHROM', 'POS', 'REF', 'ALT']].copy()
    out['#CHROM'] = out.pop('CHROM')  # Rename 'CHROM' to '#CHROM'
    out[sample_name] = input_df[sample_name, 'GT'].fillna('.').replace({0: '0', 1: '1'})

    # Add required VCF columns with default or specified values
    out['ID'], out['QUAL'], out['FILTER'], out['INFO'], out['FORMAT'] = '.', '.', '.', '.', format_str

    # Convert the DataFrame to a VCF formatted string and prepend the header
    out_str = snpad_vcf_header + out.to_csv(sep='\t', index=False)

    # Write the VCF string to a file
    with open(out_file, 'w+') as fh:
        fh.write(out_str)

def convert_vcf_dataframe_to_fasta_file(input_vcf_df, fasta_file, variant_sites_only=False, order=None):
    vcf_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

    if order:
        cols = [col for col in input_vcf_df.columns.values if col in vcf_cols] + order
        input_vcf_df = input_vcf_df[cols]

    output_vcf_df = pd.DataFrame()
    output_vcf_df = input_vcf_df[~input_vcf_df['ALT'].str.contains(",")]

    if variant_sites_only:
        output_vcf_df = output_vcf_df[output_vcf_df['ALT'] != '.']

    fasta_df = pd.DataFrame()

    # Convert genotype information to REF or ALT alleles, or 'N' for missing data
    for col in output_vcf_df.columns.values:
        if col in vcf_cols:
            continue
        temp = output_vcf_df[col].copy(deep=True).astype(object)
        temp[temp.isna()] = np.nan
        temp[temp == 0] = output_vcf_df['REF']
        temp[temp == 1] = output_vcf_df['ALT']
        temp[~temp.notna()] = 'N'
        fasta_df[col] = temp

    fasta_df = fasta_df.T

    with open(fasta_file, "w") as output_file:
        for row in fasta_df.iterrows():
            output_file.write('>' + row[0] + '\n')
            output_file.write("".join(row[1]) + '\n')

def get_genotypes_from_vcf_dataframe(input_vcf):
    genotypes = input_vcf.xs('GT', level=1, axis=1)

    output_df = pd.DataFrame()
    output_df[['CHROM', 'POS', 'REF', 'ALT']] = input_vcf[['CHROM', 'POS', 'REF', 'ALT']]
    output_df[genotypes.columns.values] = genotypes

    return output_df

def get_heterozygotes_from_vcf_dataframe(input_vcf):
    hets = input_vcf.xs('Heterozygote', level=1, axis=1)
    output_df = pd.DataFrame()
    output_df[['CHROM', 'POS', 'REF', 'ALT']] = input_vcf[['CHROM', 'POS', 'REF', 'ALT']]
    output_df[hets.columns.values] = hets

    return output_df

def prep_fasta_for_beast(input_vcf, output_fasta, identity_df_loc=None, identity_cutoff=0):
    temp = input_vcf.drop(['CHROM', 'POS', 'REF', 'ALT', 'chimp'], axis=1)

    # Remove rows where all genotypes (except chimp) are missing
    temp_df = input_vcf[~temp.isnull().all(axis=1)]
    temp_df = temp_df[~temp_df['ALT'].str.contains(",")]

    # Calculate mean across rows
    mean_values = temp_df.drop(['CHROM', 'POS', 'REF', 'ALT'], axis=1).mean(axis=1)

    # Remove the alternate allele if the mean is 0 (i.e. no alternate allele present)
    temp_df['ALT'][(mean_values == 0) & (temp_df['ALT'] != '.')] = '.'
    # Remove rows where the mean is 1 (i.e. only alternate alleles present)
    temp_df = temp_df[mean_values != 1]

    if identity_df_loc:
        identity_df = pd.read_csv(identity_df_loc, sep='\t', header=None, names=['identity'])
        identity_df['POS'] = identity_df.index + 1

        hg19_merge_identity = pd.merge(temp_df, identity_df, on=['POS'], how='inner')

        hg19_merge_identity_filtered = hg19_merge_identity[hg19_merge_identity['identity'] >= identity_cutoff]
        hg19_merge_identity_filtered = hg19_merge_identity_filtered.drop(['identity'], axis=1)

        invariant_sites = hg19_merge_identity_filtered[hg19_merge_identity_filtered['ALT'] == '.']['REF'].value_counts()
        output_vcf = hg19_merge_identity_filtered[hg19_merge_identity_filtered['ALT'] != '.']
    else:
        invariant_sites = temp_df[temp_df['ALT'] == '.']['REF'].value_counts()
        output_vcf = temp_df[temp_df['ALT'] != '.']

    # Convert VCF to FASTA
    convert_vcf_dataframe_to_fasta_file(output_vcf, output_fasta, variant_sites_only=True)

    return invariant_sites
