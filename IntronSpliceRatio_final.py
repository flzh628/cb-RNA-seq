
import re
import pysam
import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sns
from plotnine import *
from time import strftime,localtime # 20210630  localtime, or ctime: Beijing Time
import argparse

parser = argparse.ArgumentParser(description='To calculate the splicing ratio of introns from cb-RNA-seq data.',
                                epilog='Author: Fengli Zhao\tzhaofl@sustech.edu.cn      20220301 20211110 20210630 20210618 20210329 20201231 20201222 20201216')

parser.add_argument('--inputs', type=str, required=True, help='Input bam files. Separated by comma.')
parser.add_argument('--geneflt', type=str, help='FPKM matrix file of gene, used to filter the low expressed genes. [No]', default='No')
parser.add_argument('--intron', type=str, help='The intron position file. [TAIR10.LongRNA.intron]', default='TAIR10.LongRNA.intron')
parser.add_argument('--ex_chr', type=str, help='The excluded Chr(s). Else, separated by comma. [No]', default='No')
parser.add_argument('--output', type=str, required=True, help='The name of output file.')

args = parser.parse_args()

infile = args.inputs
geneflt = args.geneflt
intron = args.intron
ex_chr = args.ex_chr
output = args.output

###################################################################################################################

def readcount(bam,chrom,start0,end0):
    openbam = pysam.AlignmentFile(bam, mode='rb', threads=2)
    new_dp = openbam.count_coverage(chrom, start0, end0, quality_threshold = 20)
    read = 0; index = 0; ii = 0
    while ii < 25:
        new = int(new_dp[0][ii]) + int(new_dp[1][ii]) + int(new_dp[2][ii]) + int(new_dp[3][ii])
        if ii == 0:
            read = new
            index = new
        else:
            if new <= index:
                index = new
            else:
                increase = new - index
                read += increase
                index = new
        ii += 1
    openbam.close()
    return read

###################################################################################################################

time0 = strftime("%Y-%m-%d %H:%M:%S", localtime())
print('=' * 100)
print('\t[' + time0 + ']\tBeginning to run the IntronSpliceRatio_final.py ...\n')

SS5_dict = {}; SS3_dict = {}; ex_chr_arr = []; FPKM_dict = {}; sample_arr = []
Read_5_str = {}; Read_3_str = {}
bamRE = re.compile(r'(\S+).bam')

if ex_chr != 'No':
    ex_chr_arr = ex_chr.split(',')

if geneflt != 'No':
    openflt = open(geneflt, 'r')
    idx = 0
    while True:
        lineflt = openflt.readline()
        if not lineflt:
            break
        else:
            newarr = lineflt.strip().split('\t')
            if idx == 0:
                for bamname in newarr[1:]:
                    nmsamp = bamRE.search(bamname)
                    samp = nmsamp.group(1)
                    sample_arr.append(samp)

            else:
                idxflt = 1
                while idxflt < len(newarr):
                    FPKM_dict.setdefault(newarr[0], {})[sample_arr[idxflt - 1]] = newarr[idxflt]
                    idxflt += 1

        idx += 1

    openflt.close()

Array00 = []; Array11 = []; Array22 = []; intron_file_head = ''

bams = infile.split(',')
for bam in bams:
    time1 = strftime("%Y-%m-%d %H:%M:%S", localtime())
    print('\t\t[' + time1 + ']\tDealling with the ' + bam + ' ...')

    nmbam = bamRE.search(bam)
    name = nmbam.group(1)

    openin = open(intron, 'r')
    while True:
        line = openin.readline()
        if not line:
            break
        else:
            array = line.strip('\n').split('\t')
            if line[:3] == 'LOC':
                intron_file_head = '\t'.join(array)
            else:
                if array[2] in ex_chr_arr:
                    continue
                elif geneflt != 'No':
                    if array[0] in FPKM_dict.keys():
                        fpkm = FPKM_dict[array[0]]
                        if name in fpkm.keys():
                            rna_intron = '_'.join(array[1:])
                            if float(fpkm[name]) < 1:
                                SS5_dict.setdefault(rna_intron, {})[name] = 'unknown'
                                Read_5_str.setdefault(rna_intron, {})[name] = 'uncounted'
                                SS3_dict.setdefault(rna_intron, {})[name] = 'unknown'
                                Read_3_str.setdefault(rna_intron, {})[name] = 'uncounted'
                            else:
                                num5i = 0; num5e = 0; num3i = 0; num3e = 0
                                if array[3] == '+':
                                    start1 = int(array[4]); end1 = int(array[4]) + 25
                                    start2 = int(array[4]) - 1 - 25; end2 = int(array[4]) - 1
                                    start3 = int(array[5]) - 25; end3 = int(array[5])
                                    start4 = int(array[5]) + 1; end4 = int(array[5]) + 1 + 25
                                else:
                                    start1 = int(array[5]) - 25; end1 = int(array[5])
                                    start2 = int(array[5]) + 1; end2 = int(array[5]) + 1 + 25
                                    start3 = int(array[4]); end3 = int(array[4]) + 25
                                    start4 = int(array[4]) - 1 - 25; end4 = int(array[4]) - 1

                                num5i = readcount(bam, array[2], start1, end1)
                                num5e = readcount(bam, array[2], start2, end2)
                                num3i = readcount(bam, array[2], start3, end3)
                                num3e = readcount(bam, array[2], start4, end4)

                                Read_5_str.setdefault(rna_intron, {})[name] = str(num5e) + '_' + str(num5i)
                                Read_3_str.setdefault(rna_intron, {})[name] = str(num3e) + '_' + str(num3i)

                                if num5e == 0:
                                    SS5_dict.setdefault(rna_intron, {})[name] = 'unknown'
                                else:
                                    SS5_dict.setdefault(rna_intron, {})[name] = format(num5i/num5e, '.4f')
                                if num3e == 0:
                                    SS3_dict.setdefault(rna_intron, {})[name] = 'unknown'
                                else:
                                    SS3_dict.setdefault(rna_intron, {})[name] = format(num3i/num3e, '.4f')

                                if (num5e != 0 and num3e != 0):
                                    if ((num5i/num5e >= 0 and num5i/num5e <= 1) and (num3i/num3e >= 0 and num3i/num3e <= 1)):
                                        Array00.append(format(num5i/num5e, '.4f')); Array11.append('5SS'); Array22.append(name)
                                        Array00.append(format(num3i/num3e, '.4f')); Array11.append('3SS'); Array22.append(name)

                        else:
                            print('Please make sure the sample names same !\n')
                    
                else:
                    rna_intron = '_'.join(array[1:])
                    num5i = 0; num5e = 0; num3i = 0; num3e = 0
                    if array[3] == '+':
                        start1 = int(array[4]); end1 = int(array[4]) + 25
                        start2 = int(array[4]) - 1 - 25; end2 = int(array[4]) - 1
                        start3 = int(array[5]) - 25; end3 = int(array[5])
                        start4 = int(array[5]) + 1; end4 = int(array[5]) + 1 + 25
                    else:
                        start1 = int(array[5]) - 25; end1 = int(array[5])
                        start2 = int(array[5]) + 1; end2 = int(array[5]) + 1 + 25
                        start3 = int(array[4]); end3 = int(array[4]) + 25
                        start4 = int(array[4]) - 1 - 25; end4 = int(array[4]) - 1

                    num5i = readcount(bam, array[2], start1, end1)
                    num5e = readcount(bam, array[2], start2, end2)
                    num3i = readcount(bam, array[2], start3, end3)
                    num3e = readcount(bam, array[2], start4, end4)

                    Read_5_str.setdefault(rna_intron, {})[name] = str(num5e) + '_' + str(num5i)
                    Read_3_str.setdefault(rna_intron, {})[name] = str(num3e) + '_' + str(num3i)

                    if num5e == 0:
                        SS5_dict.setdefault(rna_intron, {})[name] = 'unknown'
                    else:
                        SS5_dict.setdefault(rna_intron, {})[name] = format(num5i/num5e, '.4f')
                    if num3e == 0:
                        SS3_dict.setdefault(rna_intron, {})[name] = 'unknown'
                    else:
                        SS3_dict.setdefault(rna_intron, {})[name] = format(num3i/num3e, '.4f')

                    if (num5e != 0 and num3e != 0):
                        if ((num5i/num5e >= 0 and num5i/num5e <= 1) and (num3i/num3e >= 0 and num3i/num3e <= 1)):
                            Array00.append(format(num5i/num5e, '.4f')); Array11.append('5SS'); Array22.append(name)
                            Array00.append(format(num3i/num3e, '.4f')); Array11.append('3SS'); Array22.append(name)

    openin.close()

time2 = strftime("%Y-%m-%d %H:%M:%S", localtime())
print('\n\t[' + time2 + ']\tFinished the computing tasks, and begin to export data ...')

openOut = open(output,'w')
openOut.write(intron_file_head)

for bam in sorted(bams):
    nmbam = bamRE.search(bam)
    name = nmbam.group(1)
    openOut.write('\t' + name + '_5SS\t' + name + '_3SS')

for bam in sorted(bams):
    nmbam = bamRE.search(bam)
    name = nmbam.group(1)
    openOut.write('\t' + name + '_Read_5ei_str\t' + name + '_Read_3ei_str')

openOut.write('\n')

SS5_Arr = {}; SS3_Arr = {}

for rna_intron,value5 in SS5_dict.items():
    intronRE = re.compile(r'((\w+).\d+)_(\w+)_([+-])_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)')
    nmnm = intronRE.search(rna_intron)
    rna = nmnm.group(1); gene = nmnm.group(2); chrid = nmnm.group(3); strand = nmnm.group(4); start = nmnm.group(5)
    end = nmnm.group(6); gene_len = nmnm.group(7); intron_num = nmnm.group(8); intron_id = nmnm.group(9)
    openOut.write(gene + '\t' + rna + '\t' + chrid + '\t' + strand + '\t' + start + '\t' + end + '\t' + gene_len + '\t' + intron_num + '\t' + intron_id)
    Array = []; ArrNew = [];
    for name in sorted(value5.keys()):
        Array.append(str(value5[name]))
        SS5_Arr.setdefault(name,[]).append(value5[name])

    if rna_intron in SS3_dict.keys():
        value3 = SS3_dict[rna_intron]
        for name in sorted(value3.keys()):
            Array.append(str(value3[name]))
            SS3_Arr.setdefault(name,[]).append(value3[name])

    readStr_arr = []; readStr_New = []
    if rna_intron in Read_5_str.keys():
        str5 = Read_5_str[rna_intron]
        for name in sorted(str5.keys()):
            readStr_arr.append(str5[name])

    if rna_intron in Read_3_str.keys():
        str3 = Read_3_str[rna_intron]
        for name in sorted(str3.keys()):
            readStr_arr.append(str3[name])

    mid = int(len(Array)/2)
    for i in range(mid):
        ArrNew.append(Array[i])
        ArrNew.append(Array[i + mid])
        readStr_New.append(readStr_arr[i])
        readStr_New.append(readStr_arr[i + mid])

    rario_str = '\t'.join(ArrNew)
    read_str = '\t'.join(readStr_New)
    openOut.write('\t' + rario_str + '\t' + read_str + '\n')

openOut.close()

SS5_Samp_arr = {}; SS3_Samp_arr = {}; iii = 0
while iii < len(Array00):
    if Array11[iii] == '5SS':
        SS5_Samp_arr.setdefault(Array22[iii],[]).append(float(Array00[iii]))
    else:
        SS3_Samp_arr.setdefault(Array22[iii],[]).append(float(Array00[iii]))
    iii += 1

print('\n\n\t\t' + '=' * 50)
print('\t\t5\'SS data:\tSample\tNumber\tMean\tMedian')
for name in SS5_Samp_arr.keys():
    arr_1 = SS5_Samp_arr[name]
    print('\t\t\t\t' + name + '\t' + str(len(arr_1)) + '\t' + str(format(np.mean(arr_1), ".4f")) + '\t' + str(np.median(arr_1)))

print('\n\t\t3\'SS data:\tSample\tNumber\tMean\tMedian')
for name in SS3_Samp_arr.keys():
    arr_2 = SS3_Samp_arr[name]
    print('\t\t\t\t' + name + '\t' + str(len(arr_2)) + '\t' + str(format(np.mean(arr_2), ".4f")) + '\t' + str(np.median(arr_2)))

print('\t\t' + '=' * 50 + '\n')


print('\t\t' + 'Wilcox test results:' + '\n')

print('\t\t\t' + '5\'SS data:')
samp_Array5 = sorted(list(SS5_Samp_arr.keys()))
samp_array51 = samp_Array5[:-1]; samp_array52 = samp_Array5[1:]
str1_str1 = '\t'.join(samp_array51)
print('\t\t\t\tSample\t' + str1_str1)
for samp52 in samp_array52:
    array_tmp = []
    array_tmp.append(samp52)
    for samp51 in samp_array51:
        if samp51 == samp52:
            break
        else:
            array_1 = SS5_Samp_arr[samp51]; array_2 = SS5_Samp_arr[samp52]
            t_val = 0; p_val = 1
            if len(array_1) == len(array_2):
                t_val, p_val = stats.wilcoxon(array_1, array_2)
            else:
                t_val, p_val = stats.mannwhitneyu(array_1, array_2)
            array_tmp.append(str(p_val))
    new_str = '\t'.join(array_tmp)
    print('\t\t\t\t' + new_str)

print('\n\t\t\t' + '3\'SS data:')
samp_Array3 = sorted(list(SS3_Samp_arr.keys()))
samp_array31 = samp_Array3[:-1]; samp_array32 = samp_Array3[1:]
str1_str1 = '\t'.join(samp_array31)
print('\t\t\t\t\t' + str1_str1)
for samp32 in samp_array32:
    array_tmp = []
    array_tmp.append(samp32)
    for samp31 in samp_array31:
        if samp31 == samp32:
            continue
        else:
            array_1 = SS3_Samp_arr[samp31]; array_2 = SS3_Samp_arr[samp32]
            t_val = 0; p_val = 1
            if len(array_1) == len(array_2):
                t_val, p_val = stats.wilcoxon(array_1, array_2)
            else:
                t_val, p_val = stats.mannwhitneyu(array_1, array_2)
            array_tmp.append(str(p_val))
    new_str = '\t'.join(array_tmp)
    print('\t\t\t\t' + new_str)

time3 = strftime("%Y-%m-%d %H:%M:%S", localtime())
print('\n\t[' + time3 + ']\tBegin to draw pictures ...\n')

df3 = pd.DataFrame({'SS_ratio':Array00, 'Group':Array11, "Sample":Array22})
df3['SS_ratio']=df3.SS_ratio.astype(float)
boxplot_fig = (ggplot(df3, aes(x='Group', y='SS_ratio', fill='Sample')) 
+ geom_boxplot(outlier_alpha=0, outlier_size=0, position=position_dodge(0.85)) 
+ labs(x="", y="Splicing ratio")
+ guides(fill=guide_legend(title=""))
+ theme_matplotlib()
+ theme(legend_position="top")
+ theme(axis_title=element_text(family="arial", size=12), axis_text=element_text(family="arial", size=10), \
        legend_title=element_text(family="arial", size=12), legend_text=element_text(family="arial", size=10)))

pdffig3 = output + '_' + '.SS_ratio.flt.boxplot.pdf'
boxplot_fig.save(pdffig3, format='pdf')

for bam in bams:
    time4 = strftime("%Y-%m-%d %H:%M:%S", localtime())
    print('\t\t[' + time4 + ']\tDraw pictures for ' + bam + ' ...')

    nmbam = bamRE.search(bam)
    name = nmbam.group(1)
    array01 = SS5_Arr[name]
    array02 = SS3_Arr[name]
    array11 = []; array22 = []; array33 = []; array44 = []; ii = 0
    while ii < len(array01):
        if (array01[ii] != 'unknown' and array02[ii] != 'unknown'):
            array11.append(array01[ii])
            array22.append(array02[ii])
            if ((float(array01[ii]) >=0 and float(array01[ii]) <= 1) and (float(array02[ii]) >=0 and float(array02[ii]) <= 1)):
                array33.append(array01[ii])
                array44.append(array02[ii])

        ii += 1

    df1 = pd.DataFrame({'x':array11, 'y':array22})
    df1['x']=df1.x.astype(float); df1['y']=df1.y.astype(float)
    sns.set_theme(font='arial')
    sns_reg1 = sns.jointplot(x='x', y='y', data=df1, color='#7CBC47', kind='reg', space=0, height=5, ratio=5,
                            scatter_kws=dict(color='#7CBC47', alpha=0.8, s=3, marker='+'),
                            line_kws=dict(color='#D31A8A',alpha=1, lw=2),
                            marginal_kws=dict(bins=20))
#    sns_reg1.plot_joint(sns.kdeplot, color="r", zorder=0, levels=6)  # 等高线
    sns_reg1.set_axis_labels(xlabel='5\'SS ratio', ylabel='3\'SS ratio')
    pdffig1 = output + '_' + name + '.SS_ratio.pdf'
    sns_reg1.savefig(pdffig1, format='pdf')

    df2 = pd.DataFrame({'x':array33, 'y':array44})
    df2['x']=df2.x.astype(float); df2['y']=df2.y.astype(float)
    sns_reg2 = sns.jointplot(x='x', y='y', data=df2, color='#7CBC47', kind='reg', space=0, height=5, ratio=5,
                            scatter_kws=dict(color='#7CBC47', alpha=0.8, s=3, marker='+'),
                            line_kws=dict(color='#D31A8A',alpha=1, lw=2),
                            marginal_kws=dict(bins=20))
#    sns_reg2.plot_joint(sns.kdeplot, color="r", zorder=0, levels=6)
    sns_reg2.set_axis_labels(xlabel='5\'SS ratio', ylabel='3\'SS ratio')
    pdffig2 = output + '_' + name + '.SS_ratio.flt.pdf'
    sns_reg2.savefig(pdffig2, format='pdf')

time5 = strftime("%Y-%m-%d %H:%M:%S", localtime())
print('\n\t[' + time5 + ']\tFinished ! Please check the data and pictures !')
print('=' * 100)
