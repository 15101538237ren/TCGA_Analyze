# -*- coding:utf-8 -*-
import os, math, re, json, time, pickle, random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#preparation files
data_dir = "data"
fig_dir = "figure"
methy_data_dir = "methy_data"
mean_and_var_fig_dir = "mean_and_std_figure"

manifest_path = data_dir + os.sep + "manifest.txt"
json_file_path = data_dir + os.sep + "metadata.json"
tumor_suppressed_gene_file = data_dir + os.sep + "gene_with_protein_product.tsv"

cancer_names = ["BRCA","COAD","LIHC","LUAD","LUSC"] #
cancer_markers = {"BRCA":'rs', "COAD":'gp', "LIHC":'bo',"LUAD":'kx',"LUSC":'c*'}

tumor_stages = ["normal","i","ia","ib","ii","iia","iib","iic","iii","iiia","iiib","iiic","iv","iva","ivb","x","not reported"]
tumor_stage_convert = {"normal":"normal","i":"i","ia":"i","ib":"i","ii":"ii","iia":"ii","iib":"ii","iic":"ii","iii":"iii","iiia":"iii","iiib":"iii","iiic":"iii","iv":"iv","iva":"iv","ivb":"iv","x":"x","not reported":"not reported"}
merged_stage = ["normal","i","ii","iii","iv","x","not reported"]
normal_keyword = "normal"

tumor_stages_xaxis = {}
for idx, item in enumerate(tumor_stages):
    tumor_stages_xaxis[item] = idx + 1

tumor_stages_xaxis2 = {}
for idx, item in enumerate(merged_stage):
    tumor_stages_xaxis2[item] = idx + 1

#checked
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

#checked
def connect_filename_to_uuid():
    uuid_to_filename = {}
    filename_to_uuid = {}
    now_file = open(manifest_path,'r')

    #pass the header
    now_file.readline()

    str_pattern = r'([^\t]+)\t([^\t]+)'

    line = now_file.readline()
    while line:
        match_p = re.search(str_pattern, line)
        if match_p:
            uuid = match_p.group(1)
            file_name = match_p.group(2)
            uuid_to_filename[uuid] = file_name
            filename_to_uuid[file_name] = uuid
            # print "%s\t%s" % (uuid, file_name)

        line=now_file.readline()
    now_file.close()
    print "connect_filename_to_uuid called"
    return [uuid_to_filename, filename_to_uuid]

#global vars
[uuid_to_filename, filename_to_uuid] = connect_filename_to_uuid()
[TSG, alias_dict] = read_whole_genenames(tumor_suppressed_gene_file)
print "gene number %d" % len(TSG)
#checked
def get_exist_uuids_from_filenames(filenames):
    uuids = []

    for filename in filenames:
        uuids.append(filename_to_uuid[filename])

    print "get_exist_uuids_from_filenames called"
    return uuids
#checked
def connect_uuid_to_cancer_stage(uuids, json_file_path):
    stage_to_uuids = {}
    uuid_to_stage = {}

    json_obj = json.load(open(json_file_path,'r'))

    for uuid in uuids:
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
            stage_save = normal_keyword
        uuid_to_stage[uuid] = stage_save
        if stage_save not in stage_to_uuids.keys():
            stage_to_uuids[stage_save] = [uuid]
        else:
            stage_to_uuids[stage_save].append(uuid)
        # print stage_save
    return [uuid_to_stage, stage_to_uuids]

#checked
def plot_for_each_gene(cancer_name, gene_name, x, y, box_data, c, xrange, xticks):
    plt.clf()
    plt.cla()
    fig, ax = plt.subplots()
    plt.xticks(xrange, xticks)

    ax.scatter(x, y, color=c, s=2.0)
    ax.boxplot(box_data,sym='')
    ax.set_xlim([0, len(merged_stage) - 0.5])
    ax.set_ylim([0, 1.0])
    ax.set_title(gene_name + " methylation for different cancer stage")

    figure_dir = fig_dir + os.sep + cancer_name
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)
    plt.savefig(figure_dir +os.sep + gene_name.lower() + '.png')

#put array data to .csv file, one value per line
def save_data_to_file(arr, path, precision = 4):
    file_out = open(path, "w")
    for item in arr:
        file_out.write(str(round(item, precision)) + "\n")
    file_out.close()

#checked
def gene_and_cancer_stage_profile_of_dna_methy(cancer_name, data_path, pickle_filepath,uuids, load=False, whole_genes= True):

    if not load:
        profile = {}
        x_profile = {}
        y_profile = {}
        for gene in TSG:
            profile[gene] = []
            x_profile[gene] = [] #profile 中 stage 编号的list
            y_profile[gene] = [] #profile 中 DNA methylation level的list
            for tumor_stage in tumor_stages:
                profile[gene].append([])

        [uuid_to_stage, stage_to_uuids] = connect_uuid_to_cancer_stage(uuids, json_file_path)
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
            ro = random.random()*0.3 - 0.15
            for gene_symbol in TSG:
                mean_methy_level_of_this_case = float(np.array(temp_gene_methy_dict[gene_symbol]).mean())
                profile[gene_symbol][tumor_stages_xaxis[uuid_to_stage[uuid]] - 1].append(mean_methy_level_of_this_case)
                x_profile[gene_symbol].append(tumor_stages_xaxis[uuid_to_stage[uuid]] + ro)
                y_profile[gene_symbol].append(mean_methy_level_of_this_case)

            t1 = time.time()
            ratio = float(uidx + 1.0) / float(len(uuids))
            percentage = ratio * 100.0
            time_this_turn = t1 - t0
            tot_timelapse = tot_timelapse + time_this_turn
            remain_time = (tot_timelapse / ratio) * (1.0 - ratio)
            print "%d/%d, %.2f %%, Total Time: %.2fs, Time Left: %.2fs" %(uidx + 1.0, len(uuids), percentage, tot_timelapse, remain_time)

        pickle_file = open(pickle_filepath, 'wb')
        pickle.dump(profile,pickle_file,-1)
        pickle.dump(x_profile,pickle_file,-1)
        pickle.dump(y_profile,pickle_file,-1)
        pickle_file.close()
    else:
        pickle_file = open(pickle_filepath,"rb")
        profile = pickle.load(pickle_file)
        x_profile = pickle.load(pickle_file)
        y_profile = pickle.load(pickle_file)
        pickle_file.close()
        print "load pickle file %s finished" % (pickle_filepath)
    return [profile, x_profile, y_profile]

#将原来小类stage数据合并到大类stage
def convert_origin_to_new_profile(origin_list_of_a_gene):
    new_list = [[] for item in merged_stage]
    for idx, item1 in enumerate(tumor_stages):
        for item2 in origin_list_of_a_gene[idx]:
            new_list[tumor_stages_xaxis2[tumor_stage_convert[item1]] - 1].append(item2)
    return new_list

def convert_origin_profile_into_merged_profile(origin_profile_list):
    origin_profile = origin_profile_list[0]
    new_profile = {gene:[] for gene in TSG}

    for gene in TSG:
        for stage in merged_stage:
            new_profile[gene].append([])
        for idx, item1 in enumerate(tumor_stages):
            for item2 in origin_profile[gene][idx]:
                new_profile[gene][tumor_stages_xaxis2[tumor_stage_convert[item1]] - 1].append(item2)
    return [new_profile]
def merged_stage_scatter_and_box_plot(profile_arr):
    profile = profile_arr[0]
    new_profile = {}

    for gene_idx, gene in enumerate(TSG):
        xs = range(1,len(merged_stage)+1)
        new_profile[gene] = convert_origin_to_new_profile(profile[gene])

        new_x_profile = []
        new_y_profile = []
        for idx, arr in enumerate(new_profile[gene]):
            for item1 in arr:
                ro = random.random() * 0.4 - 0.2
                new_x_profile.append(idx + 1 + ro)
                new_y_profile.append(item1)
            plot_for_each_gene(cancer_name, gene, new_x_profile, new_y_profile, new_profile[gene], "blue", xs, merged_stage)

def init_mean_and_var_dict():
    means_dict = {}
    variance_dict = {}
    for gene in TSG:
        means_dict[gene] = {}
        variance_dict[gene] = {}
        for cancer_name in cancer_names:
            means_dict[gene][cancer_name] = []
            variance_dict[gene][cancer_name]=[]
    return [means_dict, variance_dict]

#initialize global vars
[mean_dict, var_dict] = init_mean_and_var_dict()


def plot_mean_and_var(all_cancer_profiles):
    mean_path = mean_and_var_fig_dir + os.sep + "mean" + os.sep
    var_path = mean_and_var_fig_dir + os.sep + "std" + os.sep
    paths = [mean_path, var_path]
    for item in paths:
        if not os.path.exists(item):
            os.makedirs(item)

    for idx, cancer_name in enumerate(cancer_names):
        profile = all_cancer_profiles[idx]
        for gene in TSG:
            for idx_of_stage, stage_name in enumerate(merged_stage):
                methy_of_this_gene = profile[gene][idx_of_stage]
                #print "cancer %s gene %s stage %s len %d" % (cancer_name, gene, stage_name, len(methy_of_this_gene))
                if len(methy_of_this_gene):
                    mean_dict[gene][cancer_name].append(np.array(methy_of_this_gene).mean())
                    var_dict[gene][cancer_name].append(np.array(methy_of_this_gene).std())
                else:
                    mean_dict[gene][cancer_name].append(0.0)
                    var_dict[gene][cancer_name].append(0.0)
    dicts = [mean_dict, var_dict]
    dict_names = ["mean", "std"]
    if not os.path.exists(methy_data_dir):
        os.makedirs(methy_data_dir)
    for gene in TSG:
        for idx, dict in enumerate(dicts):
            dict_gene_i = dicts[idx][gene]
            plt.clf()
            plt.cla()
            fig, ax = plt.subplots()
            plt.xticks(range(1, len(merged_stage) - 1), merged_stage[0:-2])

            ax.set_title(dict_names[idx] + " methy-level of " + gene)
            ax.set_xlim([0, len(merged_stage) - 1.5])
            if idx == 1:
                ax.set_ylim([0, 0.2])
            else:
                ax.set_ylim([0, 1.0])
            plots = []
            for cancer_name in cancer_names:
                out_data_path = methy_data_dir + os.sep + gene + "_" + dict_names[idx] + "_" + cancer_name + ".dat"
                vals = dict_gene_i[cancer_name]
                x = range(1, len(vals) - 1)
                y = vals[0 : -2]
                out_data_file = open(out_data_path, "w")
                for idx_x, x_val in enumerate(x):
                    ltw = str(round(float(x_val), 2)) + "," + str(round(float(y[idx_x]), 4)) + "\n"
                    out_data_file.write(ltw)
                out_data_file.close()

                plot, = plt.plot(x, y, cancer_markers[cancer_name], ls="-")
                plots.append(plot)
            ax.legend(plots, cancer_names, loc='best')
            plt.savefig(paths[idx] + gene.lower() + '.png')
def save_gene_methy_data(cancer_name, profile_list):
    if not os.path.exists(methy_data_dir):
        os.makedirs(methy_data_dir)
    for gene in TSG:
        gene_data = profile_list[0][gene]
        merged_data = []
        out_xy_path = methy_data_dir + os.sep + gene + "_xy_" + cancer_name + ".dat"
        out_y_label_path = methy_data_dir + os.sep + gene + "_y_label_" + cancer_name + ".dat"
        out_xy_file = open(out_xy_path, "w")
        out_y_label_file = open(out_y_label_path, "w")
        for idx, stage in enumerate(merged_stage):
            if stage != "not reported":
                methy_cases_vals = gene_data[idx]
                for item_y in methy_cases_vals:
                    ro = random.random()*0.3 - 0.15
                    x = idx + 1 + ro
                    ltw = str(round(float(x), 4)) + "," + str(round(float(item_y), 6)) + "\n"
                    out_xy_file.write(ltw)
                    tmp_stage = "n" if stage == "normal" else stage
                    ltw2 = str(round(float(item_y), 6)) + "," + tmp_stage.ljust(4) + "\n"
                    out_y_label_file.write(ltw2)
                stage_data = gene_data[idx]
                merged_data.extend(stage_data)
                save_data_to_file(stage_data, methy_data_dir + os.sep + gene + "_" + merged_stage[idx] + "_" + cancer_name + ".dat")
        out_xy_file.close()
        out_y_label_file.close()
        save_data_to_file(merged_data,  methy_data_dir + os.sep + gene + "_" + "all_stage" + "_" + cancer_name + ".dat")
    print "save methy data successfully!"
if __name__ == '__main__':
    all_cancer_profiles = []
    for cancer_name in cancer_names:
        print "now analysing %s" % cancer_name
        data_path = data_dir + os.sep+ cancer_name + os.sep
        pickle_filepath = data_dir + os.sep + cancer_name + ".pkl"
        filenames = os.listdir(data_path)
        uuids = get_exist_uuids_from_filenames(filenames)
        temp_profile_list = gene_and_cancer_stage_profile_of_dna_methy(cancer_name,data_path, pickle_filepath,uuids, load=False, whole_genes= True)
        # new_profile_list = convert_origin_profile_into_merged_profile(temp_profile_list)
        # save_gene_methy_data(cancer_name, new_profile_list)
        #all_cancer_profiles.append(new_profile_list[0])
    #plot_mean_and_var(all_cancer_profiles)