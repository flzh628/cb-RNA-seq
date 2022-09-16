
import re
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='To convert intronSS-ratio to gene level. Used for TAIR10.LongRNA.intron in SS-ratio calculation.',
                                 epilog='Author: Fengli Zhao\tzhaofl@sustech.edu.cn      20220328 20220317')

parser.add_argument('--inSS', type=str, required=True, help='Input file. The result of intronPos_SpliceRatio_final.py ')
parser.add_argument('--prefix', type=str, required=True, help='The prefix of output names.')

args = parser.parse_args()

inSS = args.inSS
prefix = args.prefix

stRE = re.compile(r'(\S+)_(Read_\d)ei_str')

line0 = ''; SS5_ex = {}; SS5_in = {}; SS3_ex = {}; SS3_in = {}; Sample = []; Type = []; ii = 0
INss = open(inSS, 'r')
OUTput = open(prefix + '.gene.SS-ratio', 'w')
while True:
    line = INss.readline()
    if not line:
        break
    arr = line.strip().split()
    arr0 = arr[9:]
    arr1 = arr0[int(len(arr0)/2):]
    if 'uncounted' not in arr1:
        if ii == 0:
            line0 = '\t'.join(np.array(arr).take([0,1,2,3,6,7]).tolist()) # only for str in arr
            OUTput.write(line0)
            for str0 in arr1:
                stSE = stRE.search(str0)
                Sample.append(stSE.group(1)); Type.append(stSE.group(2))
                OUTput.write('\t' + stSE.group(1) + '_' + stSE.group(2).replace('Read_', 'SS'))
            OUTput.write('\t' + '\t'.join(arr1) + '\n')
        elif ii ==1:
            line0 = '\t'.join(np.array(arr).take([0,1,2,3,6,7]).tolist())
            for i in range(len(arr1)):
                num1, num2 = arr1[i].split('_')
                if Type[i] == 'Read_5':
                    SS5_ex.setdefault(Sample[i], []).append(int(num1))
                    SS5_in.setdefault(Sample[i], []).append(int(num2))
                else:
                    SS3_ex.setdefault(Sample[i], []).append(int(num1))
                    SS3_in.setdefault(Sample[i], []).append(int(num2))
        else:
            line1 = '\t'.join(np.array(arr).take([0,1,2,3,6,7]).tolist())
            if line1 == line0:
                for i in range(len(arr1)):
                    num1, num2 = arr1[i].split('_')
                    if Type[i] == 'Read_5':
                        SS5_ex.setdefault(Sample[i], []).append(int(num1))
                        SS5_in.setdefault(Sample[i], []).append(int(num2))
                    else:
                        SS3_ex.setdefault(Sample[i], []).append(int(num1))
                        SS3_in.setdefault(Sample[i], []).append(int(num2))
            else:
                OUTput.write(line0); numStr_arr = []
                for i in range(len(Sample)):
                    tot_ex = 0; tot_in = 0
                    if Type[i] == 'Read_5':
                        tot_ex = sum(SS5_ex[Sample[i]])
                        tot_in = sum(SS5_in[Sample[i]])
                    else:
                        tot_ex = sum(SS3_ex[Sample[i]])
                        tot_in = sum(SS3_in[Sample[i]])
                    numStr_arr.append(str(tot_ex) + '_' + str(tot_in))
                    if tot_ex == 0:
                        OUTput.write('\tunknown')
                    else:
                        ratio = tot_in / tot_ex
                        if ratio > 1:
                            ratio = 1
                        OUTput.write('\t' + format(ratio, '.4f'))
                OUTput.write('\t' + '\t'.join(numStr_arr) + '\n')
                line0 = '\t'.join(np.array(arr).take([0,1,2,3,6,7]).tolist())
                SS5_ex = {}; SS5_in = {}; SS3_ex = {}; SS3_in = {}
                for i in range(len(arr1)):
                    num1, num2 = arr1[i].split('_')
                    if Type[i] == 'Read_5':
                        SS5_ex.setdefault(Sample[i], []).append(int(num1))
                        SS5_in.setdefault(Sample[i], []).append(int(num2))
                    else:
                        SS3_ex.setdefault(Sample[i], []).append(int(num1))
                        SS3_in.setdefault(Sample[i], []).append(int(num2))
        ii += 1
OUTput.write(line0); numStr_arr = []
for i in range(len(Sample)):
    tot_ex = 0; tot_in = 0
    if Type[i] == 'Read_5':
        tot_ex = sum(SS5_ex[Sample[i]])
        tot_in = sum(SS5_in[Sample[i]])
    else:
        tot_ex = sum(SS3_ex[Sample[i]])
        tot_in = sum(SS3_in[Sample[i]])
    numStr_arr.append(str(tot_ex) + '_' + str(tot_in))
    if tot_ex == 0:
        OUTput.write('\tunknown')
    else:
        ratio = tot_in / tot_ex
        if ratio > 1:
            ratio = 1
        OUTput.write('\t' + format(ratio, '.4f'))
OUTput.write('\t' + '\t'.join(numStr_arr) + '\n')

INss.close()
OUTput.close()
