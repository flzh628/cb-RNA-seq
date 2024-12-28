
import os
from scipy import stats
from PvalueAdjust import BH_qvalues
import argparse

parser = argparse.ArgumentParser(description='To identify the differential SS-ratio introns.',
                                 epilog='Author: Fengli Zhao\tzhaofl@sustech.edu.cn      20241228 20220801 20220401')

parser.add_argument('--ss_res', type=str, required=True, help='The file of two-sample SS-ratio result.')
parser.add_argument('--tot_Col', type=int, required=True, help='The total read number of Col-BAM.')
parser.add_argument('--tot_mut', type=int, required=True, help='The total read number of mut-BAM.')
parser.add_argument('--names', type=str, required=True, help='The names of Col and mut, separated by commas.')

args = parser.parse_args()

ss_res = args.ss_res
tot_Col = args.tot_Col
tot_mut = args.tot_mut
names = args.names

factor = tot_Col / tot_mut
Col, mut = names.split(',')

P_vals = []
INss = open(ss_res, 'r')
OUTput = open(ss_res + '.diff', 'w')
while True:
    line = INss.readline()
    if not line:
        break
    arr = line.strip().split()
    if 'LOC' in arr[0]:
        OUTput.write(line.strip() + '\t' + Col + '_Read_ei_str\t' + mut + '_Read_ei_str_norm\tFold_change\tP_val\n')
    else:
        if arr[9] == arr[10] and arr[9] == 'unknown':
            continue
        elif arr[11] == arr[12] and arr[11] == 'unknown':
            continue
        else:
            OUTput.write(line.strip())
            num1_5e, num1_5i = arr[-4].split('_'); num1_3e, num1_3i = arr[-3].split('_')
            num2_5e, num2_5i = arr[-2].split('_'); num2_3e, num2_3i = arr[-1].split('_')
            num1_e = int(num1_5e) + int(num1_3e); num1_i = int(num1_5i) + int(num1_3i)
            num2_e = int((int(num2_5e) + int(num2_3e)) * factor + 0.5); num2_i = int((int(num2_5i) + int(num2_3i)) * factor + 0.5)
            str1 = str(num1_e) + '_' + str(num1_i); str2 = str(num2_e) + '_' + str(num2_i)
            p_val = stats.fisher_exact([[num2_i, num2_e], [num1_i, num1_e]], alternative = 'two-sided')[1] #20241228: greater -> two-sided
            OUTput.write('\t' + str1 + '\t' + str2 + '\t' + format((num2_i/num2_e + 0.0001)/(num1_i/num1_e + 0.0001), '.2f') + '\t' + str(p_val) + '\n')
            P_vals.append(str(p_val))

INss.close()
OUTput.close()

Q_vals = BH_qvalues(P_vals, 'Not')

INres = open(ss_res + '.diff', 'r')
OUTres = open(ss_res + '.diff.q', 'w')
ii = 0
while True:
    line = INres.readline()
    if not line:
        break
    if ii == 0:
        OUTres.write(line.strip() + '\tq_value\n')
    else:
        OUTres.write(line.strip() + '\t' + str(Q_vals[ii - 1]) + '\n')
    ii += 1
INres.close()
OUTres.close()

os.remove(ss_res + '.diff')
os.rename(ss_res + '.diff.q', ss_res + '.diff')
