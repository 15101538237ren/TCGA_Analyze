# -*- coding:utf-8 -*-
import os, re

#获取所有染色体的dna序列,存储在sequence_dict对应chr_i的字典中
def get_all_dna_sequences(dna_dir, file_pre, chr_list):
    sequence_dict = {}
    for chr_i in chr_list:
        chr_i = str(chr_i)
        dna_path = dna_dir + os.sep + file_pre + chr_i + ".fa"
        dna_file = open(dna_path, "r")
        line_seq = dna_file.readline()
        match = re.search(r'chromosome:([^:]+):(\d*):(\d*):(\d*)',line_seq)
        if match:
            sequence_length = int(match.group(4))
            print "chr %s len : %d" % (chr_i, sequence_length)

            seq_arr = []
            while line_seq:
                line_seq = dna_file.readline()
                seq_arr.append(line_seq[0 : -1])
            sequence_dict[chr_i] = "".join(seq_arr)
        print "finish chr %s" % chr_i
    return [sequence_dict]

# 查询chr_i (int)在区间(start, end)的序列, sequence_list为get_all_dna_sequences返回的[sequence_dict]
def query_a_sequence(sequence_list, chr_i, start, end):
    return sequence_list[0][str(chr_i)][start + 1: end + 1]

if __name__ == '__main__':
    pass