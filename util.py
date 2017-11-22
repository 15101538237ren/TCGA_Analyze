# -*- coding:utf-8 -*-
import os, math, re, json, time, pickle, random
import numpy as np
# import matplotlib
# matplotlib.use('Agg') # or matplotlib.use("Pdf")
import matplotlib.pyplot as plt

#preparation files
data_dir = "data"
pickle_dir = "pkl"
fig_dir = "figure"
methy_data_dir = "methy_data"
tsv_dir = "tsv_data"
mean_and_var_fig_dir = "mean_and_std_figure"
run_needed = "run_needed"
dirs = [data_dir, pickle_dir, fig_dir, tsv_dir, methy_data_dir, mean_and_var_fig_dir, run_needed]
for dir_name in dirs:
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

manifest_path = run_needed + os.sep + "5_cancer_manifest.tsv"
json_file_path = run_needed + os.sep + "5_cancer_meta.json"
tumor_suppressed_gene_file = run_needed + os.sep + "gene_with_protein_product.tsv"

cancer_names = ["BRCA", "COAD", "LIHC", "LUAD", "LUSC"] #,"BLCA" ,"ESCA","HNSC" ,"KIRC", "KIRP", "PAAD", "READ", "THCA", "STAD","LGG","OV","GBM","LAML", "PRAD","UCEC","SARC", "UVM","CESC", "DLBC"

tumor_stages = ["normal","i","ia","ib","ii","iia","iib","iic","iii","iiia","iiib","iiic","iv","iva","ivb","ivc","x","not reported"]
tumor_stages_rm_ivc = ["normal","i","ia","ib","ii","iia","iib","iic","iii","iiia","iiib","iiic","iv","iva","ivb","x","not reported"]
tumor_stage_convert = {"normal":"normal","i":"i","ia":"i","ib":"i","ii":"ii","iia":"ii","iib":"ii","iic":"ii","iii":"iii","iiia":"iii","iiib":"iii","iiic":"iii","iv":"iv","iva":"iv","ivb":"iv","ivc":"iv","x":"x","not reported":"not reported"}
merged_stage = ["normal","i","ii","iii","iv","x","not reported"]

tumor_stages_xaxis = {}
for idx, item in enumerate(tumor_stages):
    tumor_stages_xaxis[item] = idx + 1

tumor_stages_xaxis2 = {}
for idx, item in enumerate(merged_stage):
    tumor_stages_xaxis2[item] = idx + 1

def read_whole_genenames(file_path):
    genes = []
    alias_dict = {}
    now_file = open(file_path,'r')
    lines = now_file.readline().split("\r")
    for line in lines:
        contents = line.split("\t")
        gene_name = contents[0]
        alias_dict[gene_name] = gene_name
        alias = contents[1].split("|")
        previous_symbols = contents[2].split("|")
        if len(alias) and alias[0] != "":
            for alias_name in alias:
                alias_dict[alias_name] = gene_name
        if len(previous_symbols) and previous_symbols[0]!="":
            for previous_symbol in previous_symbols:
                alias_dict[previous_symbol] = gene_name
        genes.append(gene_name)
    now_file.close()
    return [genes, alias_dict]

[TSG, alias_dict] = read_whole_genenames(tumor_suppressed_gene_file)
print "genome gene counts %d" % len(TSG)


# 通过manifest文件中的对应关系,将下载的文件名filename和uuid对应起来,方便互相查询(uuid->filename, filename->uuid)
def connect_filename_to_uuid():
    uuid_to_filename = {}
    filename_to_uuid = {}
    now_file = open(manifest_path,'r')

    #pass the header
    now_file.readline()

    str_pattern = r'([^\t]+)\t([^\t]+)'
    cancer_pattern = r'jhu-usc.edu_([^\.]+)*'
    uuid_dict = {cancer_name:[] for cancer_name in cancer_names}
    file_name_dict={cancer_name:[] for cancer_name in cancer_names}

    line = now_file.readline()
    while line:
        match_p = re.search(str_pattern, line)
        if match_p:
            uuid = match_p.group(1)
            file_name = match_p.group(2)

            uuid_to_filename[uuid] = file_name
            filename_to_uuid[file_name] = uuid
            try:
                cancer_name = re.search(cancer_pattern, file_name).group(1)
                uuid_dict[cancer_name].append(uuid)
                file_name_dict[cancer_name].append(file_name)
            except AttributeError:
                print file_name
            # print "%s\t%s" % (uuid, file_name)

        line=now_file.readline()
    now_file.close()
    print "connect_filename_to_uuid called"
    return [uuid_to_filename, filename_to_uuid, uuid_dict, file_name_dict]

#returned global vars from connect_filename_to_uuid()
[uuid_to_filename, filename_to_uuid, uuid_dict, file_name_dict] = connect_filename_to_uuid()

#get all the uuids from the filenames
def get_exist_uuids_from_filenames(filenames):
    uuids = []

    for filename in filenames:
        uuids.append(filename_to_uuid[filename])

    print "get_exist_uuids_from_filenames called"
    return uuids

#通过json文件, 获取每个uuid对应的一级癌症阶段merged_tumor_stage和二级癌症阶段tumor_stage
def connect_uuid_to_cancer_stage(cancer_name, uuid_list, json_file_path):
    stage_to_uuids = {stage_name:[] for stage_name in tumor_stages}
    merged_stage_to_uuids={stage_name:[] for stage_name in merged_stage}
    uuid_to_stage = {}

    json_obj = json.load(open(json_file_path,'r'))
    cnt_cases = 0
    for uuid in uuid_list:
        filename = uuid_to_filename[uuid]
        sample_type = filename.split(".")[5].split("-")[3]
        tumor_type = int(sample_type[0 : -1])
        normal = 1 if tumor_type > 9 else 0

        stage_save = "not reported"
        if not normal:
            # if cancer, detail stage classification
            for obj in json_obj:
                if obj["file_id"] != uuid:
                    continue
                else:

                    if "cases" in obj.keys():
                        if len(obj["cases"]):
                            if "diagnoses" in obj["cases"][0].keys():
                                if len(obj["cases"][0]["diagnoses"]):
                                    stage = obj["cases"][0]["diagnoses"][0]["tumor_stage"]
                                    if stage != "not reported":
                                        stage_save = stage.split(" ")[1]
                    break
        else:
            stage_save = "normal"
        cnt_cases += 1
        uuid_to_stage[uuid] = stage_save
        stage_to_uuids[stage_save].append(uuid)
        merged_stage_to_uuids[tumor_stage_convert[stage_save]].append(uuid)
    print "cancer_name %s total cases %d" % (cancer_name, cnt_cases)
    return [uuid_to_stage, stage_to_uuids, merged_stage_to_uuids]

#最重要的一个函数, 通过遍历每个下载的tcga甲基化数据文件,将需要的统计量缓存到pkl文件中, load=True直接加载这些缓存好的文件,load=False, 从头计算并缓存, whole_genes=True代表计算全基因组的统计量, 如果计算部分基因集, 如抑癌基因集,则把它设为False
#目前缓存并输出的统计量: profile[gene][stage_idx] = [all cases' methylation values in stage_idx(int) and gene]
#                   : profile_uuid[stage_name] = [all the uuids in stage_name], 每个病人一个uuid, 每个癌症阶段对应的病人uuid不同, uuid顺序与上面每个gene对应癌症阶段的甲基化水平列表的顺序相同.
# 此处输出和缓存的数据的癌症阶段都是按照二级阶段(更细节的分期),如果需要merged_stage,请用下面的convert_origin_profile_into_merged_profile,输入为此函数输出
def gene_and_cancer_stage_profile_of_dna_methy(cancer_name, data_path, pickle_filepath, uuids, load=False, whole_genes= True):

    if not load:
        profile = {}
        profile_uuid = {}
        for tumor_stage in tumor_stages:
            profile_uuid[tumor_stage] = []

        for gene in TSG:
            profile[gene] = []
            for ts in tumor_stages:
                profile[gene].append([])

        [uuid_to_stage, _, _] = connect_uuid_to_cancer_stage(cancer_name,uuids, json_file_path)
        tot_timelapse = 0.0
        for uidx, uuid in enumerate(uuids):
            t0 = time.time()
            file_path = uuid_to_filename[uuid]

            now_file = open(data_path + file_path,'r')
            now_file.readline()
            line = now_file.readline()
            temp_gene_methy_dict = {}
            for gene_symbol in TSG:
                temp_gene_methy_dict[gene_symbol] = []
            while line:
                line_contents = line.split("\t")
                try:
                    gene_symbols = line_contents[5].split(";")
                    positions_to_tss = line_contents[8].split(";")
                    beta_val = -1.0 if line_contents[1] == "NA" else float(line_contents[1])
                    gene_types = line_contents[6].split(";")
                    for idx, gene_symbol in enumerate(gene_symbols):
                        if gene_symbol != "." and (-1500 <= int(positions_to_tss[idx]) <= 1000) and beta_val > 0.0:
                            if not whole_genes:
                                if (gene_symbol in TSG):
                                    temp_gene_methy_dict[gene_symbol].append(beta_val)
                                    #one gene only add once for each cpg
                                    break
                            else:
                                if (gene_types[idx]=="protein_coding"):
                                    try:
                                        temp_gene_methy_dict[alias_dict[gene_symbol]].append(beta_val)
                                    except KeyError, e1:
                                        pass
                                        # print "KeyError : %s" % str(e1)
                                    #one gene only add once for each cpg
                                    break
                except IndexError, e:
                    print "line_contents :",
                    print line_contents
                line=now_file.readline()
            now_file.close()
            for gene_symbol in TSG:
                mean_methy_level_of_this_case = float(np.array(temp_gene_methy_dict[gene_symbol]).mean())
                profile[gene_symbol][tumor_stages_xaxis[uuid_to_stage[uuid]] - 1].append(mean_methy_level_of_this_case)
            profile_uuid[uuid_to_stage[uuid]].append(uuid)

            t1 = time.time()
            ratio = float(uidx + 1.0) / float(len(uuids))
            percentage = ratio * 100.0
            time_this_turn = t1 - t0
            tot_timelapse = tot_timelapse + time_this_turn
            remain_time = (tot_timelapse / ratio) * (1.0 - ratio)
            print "%d/%d, %.2f %%, Total Time: %.2fs, Time Left: %.2fs" %(uidx + 1.0, len(uuids), percentage, tot_timelapse, remain_time)

        pickle_file = open(pickle_filepath, 'wb')
        pickle.dump(profile,pickle_file,-1)
        pickle.dump(profile_uuid,pickle_file,-1)
        pickle_file.close()
    else:
        pickle_file = open(pickle_filepath,"rb")
        profile = pickle.load(pickle_file)
        profile_uuid = pickle.load(pickle_file)
        pickle_file.close()
        print "load pickle file %s finished" % (pickle_filepath)
    return [profile, profile_uuid]

#将缓存的数据中的二级癌症阶段按照一级阶段进行合并,输出格式与gene_and_cancer_stage_profile_of_dna_methy相同
def convert_origin_profile_into_merged_profile(origin_profile_list):
    [origin_profile, origin_profile_uuid] = origin_profile_list
    new_profile = {gene:[] for gene in TSG}
    new_profile_uuid = {stage:[] for stage in merged_stage}
    tmp_tumor_stages = tumor_stages if "ivc" in origin_profile_uuid.keys() else tumor_stages_rm_ivc
    for idx, item1 in enumerate(tmp_tumor_stages):
        new_profile_uuid[tumor_stage_convert[item1]].extend(origin_profile_uuid[item1])
    for gene in TSG:
        for stage in merged_stage:
            new_profile[gene].append([])

        for idx, item1 in enumerate(tmp_tumor_stages):
            # print "%s\t%s\t%d\t%d" % (gene, item1, len(origin_profile[gene][idx]), len(origin_profile_uuid[item1]))
            if len(origin_profile[gene][idx]) == len(origin_profile_uuid[item1]):
                for idx2, item2 in enumerate(origin_profile[gene][idx]):
                    new_profile[gene][tumor_stages_xaxis2[tumor_stage_convert[item1]] - 1].append(item2)
    return [new_profile, new_profile_uuid]

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

def run_cpg_stat_pipline():
    dna_dir = data_dir + os.sep + "GRCh38"
    file_pre = "Homo_sapiens.GRCh38.dna.chromosome."
    chr_list = range(1, 23)
    sequence_rtn = get_all_dna_sequences(dna_dir, file_pre, chr_list)

if __name__ == '__main__':
    run_cpg_stat_pipline()