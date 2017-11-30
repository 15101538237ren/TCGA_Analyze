# -*- coding:utf-8 -*-
import os, math, re, json, time, pickle, random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import requests
#preparation files
data_dir = "data"
pickle_dir = "pkl"
fig_dir = "figure"
methy_data_dir = "methy_data"
tsv_dir = "tsv_data"
mean_and_var_fig_dir = "mean_and_std_figure"
run_needed = "run_needed"
snv_base_dir = "snv"
dirs = [data_dir, pickle_dir, fig_dir, tsv_dir, methy_data_dir, mean_and_var_fig_dir, run_needed]
for dir_name in dirs:
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

manifest_path = run_needed + os.sep + "24_cancer_manifest.tsv"
json_file_path = run_needed + os.sep + "24_cancer_meta.json"
tumor_suppressed_gene_file = run_needed + os.sep + "gene_with_protein_product.tsv"

cancer_names = ["BRCA", "COAD", "LIHC", "LUAD", "LUSC","BLCA" ,"ESCA","HNSC" ,"KIRC", "KIRP", "PAAD", "READ", "THCA", "STAD","LGG","OV","GBM","LAML", "PRAD","UCEC","SARC", "UVM","CESC", "DLBC"] #,

tumor_stages = ["normal","i","ia","ib","ii","iia","iib","iic","iii","iiia","iiib","iiic","iv","iva","ivb","ivc","x","not reported"]
tumor_stages_rm_ivc = ["normal","i","ia","ib","ii","iia","iib","iic","iii","iiia","iiib","iiic","iv","iva","ivb","x","not reported"]
tumor_stage_convert = {"normal":"normal","i":"i","ia":"i","ib":"i","ii":"ii","iia":"ii","iib":"ii","iic":"ii","iii":"iii","iiia":"iii","iiib":"iii","iiic":"iii","iv":"iv","iva":"iv","ivb":"iv","ivc":"iv","x":"x","not reported":"not reported"}
merged_stage = ["normal","i","ii","iii","iv","x","not reported"]
merged_stage_n = ["normal","i","ii","iii","iv","x"]

tumor_stages_xaxis = {}
for idx, item in enumerate(tumor_stages):
    tumor_stages_xaxis[item] = idx + 1

tumor_stages_xaxis2 = {}
for idx, item in enumerate(merged_stage):
    tumor_stages_xaxis2[item] = idx + 1

def read_genenames(file_path):
    genes = []
    now_file = open(file_path,'r')
    line = now_file.readline()
    while line:
        gene_name = line[0:-1]
        genes.append(gene_name)
        line = now_file.readline()
    now_file.close()
    return [genes, []]

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

#将原来小类stage数据合并到大类stage
def convert_origin_to_new_profile(origin_list_of_a_gene):
    new_list = [[] for item in merged_stage]
    for idx, item1 in enumerate(tumor_stages):
        for item2 in origin_list_of_a_gene[idx]:
            new_list[tumor_stages_xaxis2[tumor_stage_convert[item1]] - 1].append(item2)
    return new_list

def plot_for_each_gene(cancer_name, gene_name, x, y, box_data, c, xrange, xticks, out_fig_path):
    plt.clf()
    plt.cla()
    fig, ax = plt.subplots()
    plt.xticks(xrange, xticks)

    ax.scatter(x, y, color=c, s=2.0)
    ax.boxplot(box_data,sym='')
    ax.set_xlim([0, len(merged_stage) - 0.5])
    ax.set_ylim([0, 1.0])
    ax.set_title(gene_name + " methylation for different cancer stage")
    plt.savefig(out_fig_path)

def merged_stage_scatter_and_box_plot(cancer_name, profile_arr, overwritten=False):
    profile = profile_arr[0]
    new_profile = {}
    figure_dir = fig_dir + os.sep + cancer_name
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)
    for gene_idx, gene in enumerate(TSG):
        if gene in profile.keys():
            out_fig_path = fig_dir + os.sep + cancer_name + os.sep + gene.lower() + '.png'
            if not overwritten:
                if os.path.exists(out_fig_path):
                    continue
            print "now plot scatter of %s %s" %(cancer_name, gene)
            xs = range(1,len(merged_stage)+1)
            new_profile[gene] = convert_origin_to_new_profile(profile[gene])

            new_x_profile = []
            new_y_profile = []
            for idx, arr in enumerate(new_profile[gene]):
                for item1 in arr:
                    ro = random.random() * 0.4 - 0.2
                    new_x_profile.append(idx + 1 + ro)
                    new_y_profile.append(item1)
                plot_for_each_gene(cancer_name, gene, new_x_profile, new_y_profile, new_profile[gene], "blue", xs, merged_stage,out_fig_path)

#put array data to .csv file, one value per line
def save_data_to_file(arr, path, precision = 4):
    file_out = open(path, "w")
    for item in arr:
        file_out.write(str(round(item, precision)) + "\n")
    file_out.close()

#保存cancer_name癌症,out_stage_list中阶段的DNA甲基化数据
def save_gene_methy_data(cancer_name, profile_list, out_stage_list, out_xy=False, out_all_stage=False):
    if not os.path.exists(methy_data_dir):
        os.makedirs(methy_data_dir)
    profile = profile_list[0]
    for gene in TSG:
        if gene in profile.keys():
            gene_data = profile[gene]
            merged_data = []
            ltws_xy = []
            ltws_y = []
            for idx, stage in enumerate(merged_stage):
                if stage in out_stage_list:
                    if out_xy:
                        methy_cases_vals = gene_data[idx]
                        for item_y in methy_cases_vals:
                            ro = random.random()*0.3 - 0.15
                            x = idx + 1 + ro
                            ltw = str(round(float(x), 4)) + "," + str(round(float(item_y), 6)) + "\n"
                            ltws_xy.append(ltw)
                            tmp_stage = "n" if stage == "normal" else stage
                            ltw2 = str(round(float(item_y), 6)) + "," + tmp_stage.ljust(4) + "\n"
                            ltws_y.append(ltw2)
                    stage_data = gene_data[idx]
                    merged_data.extend(stage_data)
                    save_data_to_file(stage_data, methy_data_dir + os.sep + gene + "_" + merged_stage[idx] + "_" + cancer_name + ".dat")
            if out_xy:
                out_xy_path = methy_data_dir + os.sep + gene + "_xy_" + cancer_name + ".dat"
                out_y_label_path = methy_data_dir + os.sep + gene + "_y_label_" + cancer_name + ".dat"
                out_xy_file = open(out_xy_path, "w")
                out_y_label_file = open(out_y_label_path, "w")
                out_xy_file.write("\n".join(ltws_xy))
                out_y_label_file.write("\n".join(ltws_y))
                out_xy_file.close()
                out_y_label_file.close()
            if out_all_stage:
                save_data_to_file(merged_data,  methy_data_dir + os.sep + gene + "_" + "all_stage" + "_" + cancer_name + ".dat")
    print "save methy data successfully!"

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

#将某癌症数据写入到tsv文件中
def dump_data_into_tsv_according_to_cancer_type_and_stage(cancer_name, uuid_list, outdir, profile_list):
    [profile, profile_uuid] = profile_list

    for stage_idx, stage_name in enumerate(merged_stage):
        if stage_name != "not reported" and len(profile_uuid[stage_name]):
            outfile_path = outdir + os.sep + cancer_name + "_" + stage_name + ".tsv"
            outfile = open(outfile_path, "w")
            header = "gene\t" + "\t".join(profile_uuid[stage_name]) + "\n"
            outfile.write(header)
            for gene in TSG:
                if gene in profile.keys():
                    gene_valid = True
                    methy_vals = []
                    for item in profile[gene][stage_idx]:
                        item_str = str(round(item,4))
                        if item_str == "nan":
                            gene_valid = False
                            break
                        else:
                            methy_vals.append(item_str)
                    if gene_valid:
                        data_str = gene + "\t" + "\t".join(methy_vals) + "\n"
                        outfile.write(data_str)
            outfile.close()
    print "%s dump_data_into_tsv_according_to_cancer_type_and_stage" % cancer_name

def query_stage_of_an_submitter_id(submitter_ids, query_size):
    cases_endpt = 'https://api.gdc.cancer.gov/cases'
    filt = {"op":"and","content":[{"op":"in","content":{"field":"submitter_id","value":submitter_ids}}]}
    params = {'filters':json.dumps(filt), 'expand':'diagnoses','fields':['diagnoses.submitter_id','diagnoses.tumor_stage'],'pretty':True,'size':query_size}
    # requests URL-encodes automatically
    response = requests.get(cases_endpt, params = params)
    ret_obj = response.json()
    # ret_obj = json.dumps(response.json(), indent=2)
    # print ret_obj

    # print json.dumps(ret_obj["data"]["hits"],indent=2)
    submitter_id_to_stage = {}
    for tidx in range(len(submitter_ids)):
        try:
            if "hits" in ret_obj["data"].keys() and  "diagnoses" in ret_obj["data"]["hits"][tidx].keys() and len(ret_obj["data"]["hits"][tidx]["diagnoses"]):
                stage_name = ret_obj["data"]["hits"][tidx]["diagnoses"][0]["tumor_stage"]
                submitter_id = ret_obj["data"]["hits"][tidx]["diagnoses"][0]["submitter_id"]
                stage = stage_name if stage_name == "not reported" else stage_name.split(" ")[1]
                submitter_id_rev = submitter_id.split("_")[0]
                submitter_id_to_stage[submitter_id_rev] = stage
        except IndexError,e:
            print e
    print submitter_id_to_stage
    return [submitter_id_to_stage]
def dna_mutation_data_transform_pipline():
    outdir ="snv_dat"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    #create gene_index_file
    with open(os.path.join(outdir, "gene_idx.txt"),"w") as gene_idx_file:
        gene_idx_file.write("\n".join([str(gidx+1) + "\t" + gene for gidx, gene in enumerate(TSG)]))

    colum_idxs = [0, 8, 15, 32] #要提取maf文件的列编号
    data_pkl_dict = {}
    for cancer_name in cancer_names:
        cancer_dir = os.path.join(snv_base_dir, cancer_name)
        if os.path.exists(cancer_dir):
            data_dict = {stage:{} for stage in merged_stage}
            uuid_to_stage_dict = {}
            file_names = os.listdir(cancer_dir)
            for file_name in file_names:
                if not file_name.startswith("."):
                    file_path = os.path.join(cancer_dir,file_name)
                    snv_file = open(file_path, "r")
                    print "%s start %s" % (cancer_name, file_path)
                    # pass head 6 lines
                    for i in range(6):
                        snv_file.readline()
                    line = snv_file.readline()
                    while line:
                        line_contents = line.split("\t")
                        [gene_name, variation_classification, bar_code, uuid] = [line_contents[i] for i in colum_idxs]
                        submitter_id = "-".join(bar_code.split("-")[0:3])

                        if submitter_id not in uuid_to_stage_dict.keys():
                            uuid_to_stage_dict[submitter_id] = query_stage_of_an_submitter_id(submitter_id)
                        stage = tumor_stage_convert[uuid_to_stage_dict[submitter_id]]
                        if gene_name in TSG and variation_classification == "Missense_Mutation":
                            if submitter_id not in data_dict[stage].keys():
                                data_dict[stage][submitter_id] = {gene : 0 for gene in TSG}
                            data_dict[stage][submitter_id][gene_name] += 1
                        line = snv_file.readline()
                    print "end %s" % file_path
            data_pkl_dict[cancer_name]=[data_dict, uuid_to_stage_dict]
    for cancer_name in cancer_names:
        cancer_dir = os.path.join(snv_base_dir, cancer_name)
        if os.path.exists(cancer_dir):
            for cancer_stage in merged_stage_n:
                submitter_ids = data_dict[cancer_stage].keys()
                #create sample id file
                with open(os.path.join(outdir, "Mut_"+cancer_name+"_"+cancer_stage+"_sample_id.txt"),"w") as sample_id_file:
                    sample_id_file.write("\n".join([str(gidx+1) + "\t" + submitter_id for gidx, submitter_id in enumerate(submitter_ids)]))
                with open(os.path.join(outdir, "Mut_"+cancer_name+"_"+cancer_stage+".dat"),"w") as data_file:
                    header = "\t".join([str(item) for item in range(len(submitter_ids) + 1)]) + "\n"
                    data_file.write(header)
                    data_str = []
                    for gidx,gene  in enumerate(TSG):
                        arr = [gidx + 1]
                        arr.extend([data_dict[cancer_stage][submitter_id][gene] for submitter_id in submitter_ids])
                        data_str.append("\t".join([str(item) for item in arr]))
                    data_file.write("\n".join(data_str))
def get_maf_submitter_ids():
    submitter_dict = {}
    outdir ="snv_dat"
    for cancer_name in cancer_names:
        cancer_dir = os.path.join(snv_base_dir, cancer_name)
        if os.path.exists(cancer_dir):
            output_cancer_dir = os.path.join(outdir, cancer_name)
            if not os.path.exists(output_cancer_dir):
                os.makedirs(output_cancer_dir)
            outfile_path = os.path.join(output_cancer_dir, cancer_name + "_submitter_ids.txt")
            submitter_dict[cancer_name] = {}
            file_names = os.listdir(cancer_dir)
            for file_name in file_names:
                if not file_name.startswith("."):
                    file_path = os.path.join(cancer_dir,file_name)
                    snv_file = open(file_path, "r")
                    print "%s start %s" % (cancer_name, file_path)
                    # pass head 6 lines
                    for i in range(6):
                        snv_file.readline()
                    line = snv_file.readline()
                    while line:
                        line_contents = line.split("\t")
                        bar_code = line_contents[15]
                        submitter_id = "-".join(bar_code.split("-")[0:3])
                        if submitter_id not in submitter_dict[cancer_name].keys():
                            submitter_dict[cancer_name][submitter_id] = 1
                        line = snv_file.readline()
                    print "end %s" % file_path
            with open(outfile_path,"w") as outfile:
                outfile.write("\n".join([str(gidx+1) + "\t" + submitter_id for gidx, submitter_id in enumerate(submitter_dict[cancer_name].keys())]))
                print "write %s successful" % outfile_path
def get_submitter_id_stages():
    outdir ="snv_dat"
    query_size = 300
    for cancer_name in cancer_names:
        cancer_dir = os.path.join(snv_base_dir, cancer_name)
        if os.path.exists(cancer_dir):
            output_cancer_dir = os.path.join(outdir, cancer_name)
            input_path = os.path.join(output_cancer_dir, cancer_name + "_submitter_ids.txt")
            output_path = os.path.join(output_cancer_dir, cancer_name + "_stages.txt")
            submitter_id_stage_dict = {}
            submitter_ids = []
            with open(input_path, "r") as input_file:
                line = input_file.readline()
                while line:
                    submitter_id = line.split("\t")[1][0:-1]
                    submitter_ids.append(submitter_id)
                    line = input_file.readline()
            print cancer_name
            len_submitters = len(submitter_ids)
            print len_submitters
            for i in range(0, len_submitters, query_size):
                sub_submitter_ids = submitter_ids[i:i+query_size]
                stages_of_subids_dict = query_stage_of_an_submitter_id(sub_submitter_ids, query_size)
                for k,v in stages_of_subids_dict[0].items():
                    submitter_id_stage_dict[k] = v
            print submitter_id_stage_dict

            with open(output_path,"w") as outfile:
                for gidx, submitter_id in enumerate(submitter_ids):
                    if submitter_id not in submitter_id_stage_dict.keys():
                        submitter_id_stage_dict[submitter_id] = "not reported"
                outfile.write("\n".join([str(gidx+1)  + "\t" + submitter_id_stage_dict[submitter_id] for gidx, submitter_id in enumerate(submitter_ids)])) #+ "\t" + submitter_id
                print "write %s successful" % input_path
if __name__ == '__main__':
    get_maf_submitter_ids()
    get_submitter_id_stages()
    # dna_mutation_data_transform_pipline()