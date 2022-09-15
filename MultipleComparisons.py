import numpy as np
import scipy.stats as stats
from PvalueAdjust import BH_qvalues

def OneFactor_multiComparison(df_in, factor_nm, value_nm, condition0, stat_file, pair_ornot="Non", stat_style="Overwrite"): # 20220523: pair_ornot; 20220810: stat_style
    Group_f = sorted(list(set(df_in[factor_nm])))
    grpNum = {} # 20220704
    for gf in Group_f:
        df_gf = df_in[df_in[factor_nm] == gf]
        grpNum[gf] = len(df_gf)

    OUTstat = '' # 20220810
    if stat_style == 'Append':
        OUTstat = open(stat_file, 'a')
    else:
        OUTstat = open(stat_file, 'w')
    Pval_dict = {}; result = ''
    while True:
        group_f0 = Group_f[0]
        if len(Group_f) == 1:
            break
        Group_f = Group_f[1:]
        OUTstat.write('=' * 50 + ' Multiple comparisons of ' + value_nm + ' between different ' + factor_nm + ', under the ' + condition0 + ' ' + '=' * 50 + '\n\n')
        df_in0 = df_in[df_in[factor_nm] == group_f0]
        if len(df_in0) >= 3: # 20220424
            result0 = stats.shapiro(df_in0[value_nm])
            OUTstat.write('\tThe normal-distribution test of values in ' + group_f0 + '\n')
            OUTstat.write('\tMean: ' + str(format(np.mean(df_in0[value_nm]), '.4f')) + '\tMedian: ' + str(format(np.median(df_in0[value_nm]), '.4f')) + \
                          '\tNumber: ' + str(len(df_in0[value_nm])) + '\n')
            OUTstat.write('\t' + str(result0) + '\n\n')
            for group_f1 in Group_f:
                df_in1 = df_in[df_in[factor_nm] == group_f1]
                if len(df_in1) >= 3: # 20220424
                    result1 = stats.shapiro(df_in1[value_nm])
                    OUTstat.write('\t' + '=' * 90 + '\n\n')
                    OUTstat.write('\t\tThe normal-distribution test of values in ' + group_f1 + '\n')
                    OUTstat.write('\t\tMean: ' + str(format(np.mean(df_in1[value_nm]), '.4f')) + '\tMedian: ' + str(format(np.median(df_in1[value_nm]), '.4f')) + \
                                  '\tNumber: ' + str(len(df_in1[value_nm])) + '\n')
                    OUTstat.write('\t\t' + str(result1) + '\n\n')
                    if pair_ornot == 'Non':
                        result = stats.mannwhitneyu(df_in0[value_nm], df_in1[value_nm], alternative = 'two-sided')
                    else:
                        result = stats.wilcoxon(df_in0[value_nm], df_in1[value_nm], alternative = 'two-sided')
                    OUTstat.write('\t\tThe differential analysis between ' + group_f0 + ' and ' + group_f1 + ':\n\n')
                    OUTstat.write('\t\t' + str(result) + '\n\n')
                    samp_pair = ':'.join(sorted([group_f0, group_f1]))
                    Pval_dict[samp_pair] = float(result[1])
    OUTstat.write('=' * 100 + '\n\n\n')

    if len(Pval_dict.keys()) >= 3: # two-sample: one pval; three-sample: three pval-s
        Pval_sort = dict(sorted(Pval_dict.items(), key=lambda x:x[1], reverse=False))
        Pval_arr = [str(x) for x in Pval_sort.values()]

        Qval_arr = BH_qvalues(Pval_arr, 'Yes')
        Qval_new = [float(x) for x in Qval_arr]
        samp_pair = Pval_sort.keys()
        Qval_dict = dict(zip(samp_pair, Qval_new))

        Group_f = sorted(list(set(df_in[factor_nm])))
        samp_arr1 = Group_f[:-1]; samp_arr2 = Group_f[1:]
        OUTstat.write('\tThe Benjamini-Hochberg corrected p-value table: ( ' + condition0 + ' )\n\n')
        OUTstat.write('\t_Sample_\t' + '\t'.join(samp_arr1) + '\n')
        for samp_2 in samp_arr2:
            tmp_arr = []
            tmp_arr.append(samp_2)
            for samp_1 in samp_arr1:
                if samp_1 == samp_2:
                    break
                else:
                    samp_p0 = ':'.join(sorted([samp_1, samp_2]))
                    if samp_p0 in Qval_dict: # 20220424
                        tmp_arr.append(format(Qval_dict[samp_p0], '.2e'))
                    else:
                        tmp_arr.append('----')
            line_str = '\t'.join(tmp_arr)
            OUTstat.write('\t' + line_str + '\n')
        OUTstat.write('\n')

    OUTstat.write('\n\t' + '-' * 20 + '\n') # 20220704
    OUTstat.write('\tGroup\tLength\n')
    for gf in sorted(grpNum.keys()):
        OUTstat.write('\t' + gf + '\t' + format(grpNum[gf], ',') + '\n')
    OUTstat.write('\t' + '-' * 20 + '\n\n')

    OUTstat.close()

def TwoFactor_multiComparison(df_in, factor_nm1, factor_nm2, value_nm, stat_file, pair_ornot="Non", stat_style="Overwrite"): # 20220810: stat_style
    Group_f1 = sorted(list(set(df_in[factor_nm1])))
    Group_f2 = sorted(list(set(df_in[factor_nm2])))
    for f1_idx in range(len(Group_f1)): # 20220828
        g_factor1 = Group_f1[f1_idx]
        if f1_idx > 0: # 20220828
            stat_style = 'Append'
        df_in_f1 = df_in[df_in[factor_nm1] == g_factor1]
        if pair_ornot == 'Non':
            OneFactor_multiComparison(df_in_f1, factor_nm2, value_nm, g_factor1, stat_file, 'Non', stat_style)
        else:
            OneFactor_multiComparison(df_in_f1, factor_nm2, value_nm, g_factor1, stat_file, pair_ornot, stat_style)

    for g_factor2 in Group_f2:
        df_in_f2 = df_in[df_in[factor_nm2] == g_factor2]
        if pair_ornot == 'Non':
            OneFactor_multiComparison(df_in_f2, factor_nm1, value_nm, g_factor2, stat_file, 'Non', 'Append')
        else:
            OneFactor_multiComparison(df_in_f2, factor_nm1, value_nm, g_factor2, stat_file, pair_ornot, 'Append')

