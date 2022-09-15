import os

def main():
    print('Just for p-values adjust using BH method !')

def BH_qvalues(p_values, sort_ornot):
    if p_values == []:
        return []
    Max_rank = {}; total_0 = len(p_values); Q_val = {}; q_values = []
    if sort_ornot == 'Not':
        openTmp = open('p-tmp', 'w')
        p_val_str0 = '\n'.join(p_values)
        openTmp.write(p_val_str0 + '\n')
        openTmp.close()
        os.system('msort -k n1 p-tmp >p-tmp.msort')
        os.remove('p-tmp')
        openTmpN = open('p-tmp.msort', 'r')
        ii_0 = 0
        while True:
            lineTmpN = openTmpN.readline()
            if not lineTmpN:
                break
            ii_0 += 1
            Max_rank[float(lineTmpN.strip())] = ii_0
        openTmpN.close()
        os.remove('p-tmp.msort')
    elif sort_ornot == 'Yes':
        for ii_1 in range(total_0):
            Max_rank[float(p_values[ii_1])] = ii_1 + 1
    else:
        os._exit()

    p_uniq_sort = sorted(Max_rank.keys())
    for p_val_0 in p_uniq_sort:
        rank = Max_rank[p_val_0]
        Q_val[p_val_0] = total_0 * float(p_val_0) / rank
    for i_0 in range(len(p_uniq_sort) - 2, -1, -1):
        if Q_val[p_uniq_sort[i_0]] > Q_val[p_uniq_sort[i_0 + 1]]:
            Q_val[p_uniq_sort[i_0]] = Q_val[p_uniq_sort[i_0 + 1]]
    for mm_0 in range(total_0):
        if float(p_values[mm_0]) in Q_val.keys():
            q_values.append(float(Q_val[float(p_values[mm_0])]))

    return q_values

if __name__ == '__main__':
    main()
