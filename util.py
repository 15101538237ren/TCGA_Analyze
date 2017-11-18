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
json_file_path = run_needed + os.sep + "metadata.json"
tumor_suppressed_gene_file = run_needed + os.sep + "gene_with_protein_product.tsv"

cancer_names = ["BRCA", "COAD", "LIHC", "LUAD", "LUSC"]#, "BLCA" ,"ESCA" ,"HNSC" ,"KIRC" ,"KIRP" ,"PAAD" ,"PRAD" ,"READ" ,"THCA" ,"UCEC"] #
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
    cancer_pattern = r'jhu-usc.edu_([^\.]+)*'
    uuid_dict = {cancer_name:[] for cancer_name in cancer_names}
    file_name_dict={cancer_name:[] for cancer_name in cancer_names}

    line = now_file.readline()
    while line:
        match_p = re.search(str_pattern, line)
        if match_p:
            uuid = match_p.group(1)
            file_name = match_p.group(2)
            cancer_name = re.search(cancer_pattern, file_name).group(1)
            uuid_to_filename[uuid] = file_name
            filename_to_uuid[file_name] = uuid
            uuid_dict[cancer_name].append(uuid)
            file_name_dict[cancer_name].append(file_name)
            # print "%s\t%s" % (uuid, file_name)

        line=now_file.readline()
    now_file.close()
    print "connect_filename_to_uuid called"
    return [uuid_to_filename, filename_to_uuid, uuid_dict, file_name_dict]

#global vars
[uuid_to_filename, filename_to_uuid, uuid_dict, file_name_dict] = connect_filename_to_uuid()
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
            stage_save = normal_keyword
        cnt_cases += 1
        uuid_to_stage[uuid] = stage_save
        stage_to_uuids[stage_save].append(uuid)
        merged_stage_to_uuids[tumor_stage_convert[stage_save]].append(uuid)
    print "cancer_name %s total cases %d" % (cancer_name, cnt_cases)
    return [uuid_to_stage, stage_to_uuids, merged_stage_to_uuids]

#put array data to .csv file, one value per line
def save_data_to_file(arr, path, precision = 4):
    file_out = open(path, "w")
    for item in arr:
        file_out.write(str(round(item, precision)) + "\n")
    file_out.close()

#checked
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

#将原来小类stage数据合并到大类stage
def convert_origin_to_new_profile(origin_list_of_a_gene):
    new_list = [[] for item in merged_stage]
    for idx, item1 in enumerate(tumor_stages):
        for item2 in origin_list_of_a_gene[idx]:
            new_list[tumor_stages_xaxis2[tumor_stage_convert[item1]] - 1].append(item2)
    return new_list

def convert_origin_profile_into_merged_profile(origin_profile_list):
    [origin_profile, origin_profile_uuid] = origin_profile_list
    new_profile = {gene:[] for gene in TSG}
    new_profile_uuid = {stage:[] for stage in merged_stage}

    for gene in TSG:
        for stage in merged_stage:
            new_profile[gene].append([])

        for idx, item1 in enumerate(tumor_stages):
            # print "%s\t%s\t%d\t%d" % (gene, item1, len(origin_profile[gene][idx]), len(origin_profile_uuid[item1]))
            if len(origin_profile[gene][idx]) == len(origin_profile_uuid[item1]):
                for idx2, item2 in enumerate(origin_profile[gene][idx]):
                    new_profile[gene][tumor_stages_xaxis2[tumor_stage_convert[item1]] - 1].append(item2)
                    new_profile_uuid[tumor_stage_convert[item1]].append(origin_profile_uuid[item1][idx2])
    return [new_profile, new_profile_uuid]
#checked
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
def plot_scatter_for_one_gene(cancer_name, gene_name, profile_arr, out_dir="."):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    profile = profile_arr[0]
    out_fig_path = out_dir + os.sep + cancer_name + "_" + gene_name.lower() + '.png'
    xs = range(1,len(merged_stage)+1)
    new_profile = convert_origin_to_new_profile(profile[gene_name])
    new_x_profile = []
    new_y_profile = []
    for idx, arr in enumerate(new_profile):
        for item1 in arr:
            ro = random.random() * 0.4 - 0.2
            new_x_profile.append(idx + 1 + ro)
            new_y_profile.append(item1)
        plot_for_each_gene(cancer_name, gene_name, new_x_profile, new_y_profile, new_profile, "blue", xs, merged_stage,out_fig_path)
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
            if gene in profile.keys():
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
            if gene in dicts[idx].keys():
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
def save_cancer_std_and_mean_of_all_genes(cancer_name, cancer_profile_arr, stage_names, out_dir="."):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    cancer_profile = cancer_profile_arr[0]
    for gene in TSG:
        if gene in cancer_profile.keys():
            for idx_of_stage, stage_name in enumerate(stage_names):
                methy_of_this_gene = cancer_profile[gene][idx_of_stage]
                mean_dict[gene][cancer_name].append(np.array(methy_of_this_gene).mean())
                var_dict[gene][cancer_name].append(np.array(methy_of_this_gene).std())
    dict_names = ["mean", "std"]
    dicts = [mean_dict, var_dict]
    header = ["gene"]

    for stage_name in stage_names:
        header.append(stage_name)
    for idx, dict_val in enumerate(dicts):
        out_data_path = out_dir + os.sep + cancer_name + "_" + dict_names[idx] + ".dat"
        out_data_file = open(out_data_path, "w")
        out_data_file.write(",".join(header) + "\n")
        for gene in TSG:
            if gene in dict_val.keys():
                data = []
                flag = True
                for item in dict_val[gene][cancer_name]:
                    if item != item:
                        flag = False
                        break
                    else:
                        data.append(str(round(float(item), 4)))
                if flag:
                    ltw = gene + "," + ",".join(data) + "\n"
                    out_data_file.write(ltw)
        out_data_file.close()
    print "write std and means into %s successful" % out_dir
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
                    methy_cases_vals = gene_data[idx]
                    for item_y in methy_cases_vals:
                        ro = random.random()*0.3 - 0.15
                        x = idx + 1 + ro
                        if out_xy:
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

def get_stage_idx_from_stage_name(stage_name):
    for idx, stage in enumerate(merged_stage):
        if stage == stage_name:
            return idx
    return -1
#比较对比组和癌症组的基因数量变化
def cmp_gene_variations_in_mean_and_std(cancer_names, stats_names, stat_epsilons, control_stage, exp_stage, input_dir,out_dir="."):
    print "cmp_gene_variations_in_mean_and_std in"
    data_dict = {}
    control_stage_idx = get_stage_idx_from_stage_name(control_stage) + 1
    exp_stage_idx = get_stage_idx_from_stage_name(exp_stage) + 1
    for cancer_name in cancer_names:
        data_dict[cancer_name] = {}
        for stat_idx, stats_name in enumerate(stats_names):
            data_dict[cancer_name][stats_name] = {}

            temp_arr = [[],[],[]]

            input_file_path = input_dir + os.sep + cancer_name + "_" + stats_name + ".dat"
            input_file = open(input_file_path, "r")
            input_file.readline()
            line = input_file.readline()
            while line:
                line_contents = line.split(",")
                gene_name = line_contents[0]
                control_val = float(line_contents[control_stage_idx])
                exp_val = float(line_contents[exp_stage_idx])
                if math.fabs(exp_val - control_val) > stat_epsilons[stat_idx]:
                    if exp_val - control_val > 0:
                        temp_arr[0].append(gene_name)
                    else:
                        temp_arr[2].append(gene_name)
                else:
                    temp_arr[1].append(gene_name)
                line = input_file.readline()
            input_file.close()

            data_dict[cancer_name][stats_name]["gene_names"] = temp_arr
            data_dict[cancer_name][stats_name]["gene_counts"] = []
            gene_sum = 0
            for arr_item in temp_arr:
                gene_count = len(arr_item)
                gene_sum += gene_count
                data_dict[cancer_name][stats_name]["gene_counts"].append(gene_count)
            data_dict[cancer_name][stats_name]["gene_percentage"] = [float(item)/float(gene_sum)*100.0 for item in data_dict[cancer_name][stats_name]["gene_counts"]]
    result_names = ["gene_names", "gene_counts", "gene_percentage"]
    pre_name = control_stage + "_" + exp_stage + "_"
    for stats_name in stat_names:
        for result_name in result_names:
            outfile_path = out_dir + os.sep + pre_name + stats_name + "_" + result_name + ".dat"
            out_file = open(outfile_path, "w")
            header = ["cancer", "up", "stable", "down"]
            out_file.write(",".join(header) + "\n")
            for cancer_name in cancer_names:
                if result_name == "gene_names":
                    data_arr = [",".join(item) for item in data_dict[cancer_name][stats_name][result_name]]
                else:
                    data_arr = [str(round(item, 4)) for item in data_dict[cancer_name][stats_name][result_name]]
                out_file.write(cancer_name + "\t" + "\t".join(data_arr) + "\n")
            out_file.close()
    out_lei_classification_path = out_dir + os.sep + pre_name + "lei_classification_gene_names.dat"
    out_lei_classification = open(out_lei_classification_path, "w")
    stat_dict = {}
    up_down_names = ["up", "stable", "down"]
    out_lei_classification_count_path = out_dir + os.sep + pre_name + "lei_classification_count.dat"
    out_lei_count = open(out_lei_classification_count_path, "w")

    out_lei_classification_percentage_path = out_dir + os.sep + pre_name + "lei_classification_percentage.dat"
    out_lei_percentage = open(out_lei_classification_percentage_path, "w")

    header = ["cancer", "mean_up", "mean_down", "mean_s_std_up", "mean_s_std_s", "mean_s_std_down"]
    out_lei_classification.write(",".join(header) + "\n")
    out_lei_count.write(",".join(header) + "\n")
    out_lei_percentage.write(",".join(header) + "\n")

    for cancer_name in cancer_names:
        stat_dict[cancer_name] = {}
        for stats_name in stat_names:
            stat_dict[cancer_name][stats_name] = {}
            for up_idx, up_down_name in enumerate(up_down_names):
                stat_dict[cancer_name][stats_name][up_down_name] = data_dict[cancer_name][stats_name]["gene_names"][up_idx] #array
        stat_dict[cancer_name]["mean_up"] = stat_dict[cancer_name]["mean"]["up"]
        stat_dict[cancer_name]["mean_down"] = stat_dict[cancer_name]["mean"]["down"]
        mean_s_std_up = []
        mean_s_std_s = []
        mean_s_std_down = []
        for gene_name in stat_dict[cancer_name]["mean"]["stable"]:
            if gene_name in stat_dict[cancer_name]["std"]["up"]:
                mean_s_std_up.append(gene_name)
            elif gene_name in stat_dict[cancer_name]["std"]["stable"]:
                mean_s_std_s.append(gene_name)
            else:
                mean_s_std_down.append(gene_name)
        stat_dict[cancer_name]["mean_s_std_up"] = mean_s_std_up
        stat_dict[cancer_name]["mean_s_std_s"] = mean_s_std_s
        stat_dict[cancer_name]["mean_s_std_down"] = mean_s_std_down

        class_names = ["mean_up", "mean_down", "mean_s_std_up", "mean_s_std_s", "mean_s_std_down"]
        tot_gene_count = 0
        for class_name in class_names:
            tot_gene_count += len(stat_dict[cancer_name][class_name])
        data_arr = [",".join(stat_dict[cancer_name][class_name]) for class_name in class_names]
        out_lei_classification.write(cancer_name + "\t" + "\t".join(data_arr) + "\n")
        count_arr = [str(len(stat_dict[cancer_name][class_name])) for class_name in class_names]
        percent_arr = [str(round(float(len(stat_dict[cancer_name][class_name]))*100.0/float(tot_gene_count), 4)) for class_name in class_names]
        out_lei_count.write(cancer_name + "\t" + "\t".join(count_arr) + "\n")
        out_lei_percentage.write(cancer_name + "\t" + "\t".join(percent_arr) + "\n")
    out_lei_classification.close()
    out_lei_count.close()
    out_lei_percentage.close()
    print "cmp_gene_variations_in_mean_and_std success"
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
if __name__ == '__main__':
    all_cancer_profiles = []
    stat_names = ["mean","std"]
    stat_epsilons = [0.05, 0.05]
    out_stage_list = ["normal","i"]

    for cancer_name in cancer_names:
        print "now start %s" % cancer_name
        data_path = data_dir + os.sep+ cancer_name + os.sep
        pickle_filepath = pickle_dir + os.sep + cancer_name + ".pkl"

        #local scripts
        # filenames = os.listdir(data_path)
        # uuids = get_exist_uuids_from_filenames(filenames)
        # temp_profile_list = gene_and_cancer_stage_profile_of_dna_methy(cancer_name,data_path, pickle_filepath,uuids, load=False, whole_genes= True)

        #server script
        temp_profile_list = gene_and_cancer_stage_profile_of_dna_methy(cancer_name,data_path, pickle_filepath,uuid_dict[cancer_name], load=True, whole_genes= True)
        # save_cancer_std_and_mean_of_all_genes(cancer_name, temp_profile_list, [normal_keyword, "i"],out_dir="mean_std_data")
    # cmp_gene_variations_in_mean_and_std(cancer_names, stat_names, stat_epsilons, normal_keyword, "i", input_dir ="mean_std_data", out_dir="stat")
        # merged_stage_scatter_and_box_plot(cancer_name, temp_profile_list, overwritten=False)
        new_profile_list = convert_origin_profile_into_merged_profile(temp_profile_list)
        save_gene_methy_data(cancer_name, new_profile_list, out_stage_list, out_xy=False, out_all_stage=False)
        out_dir = tsv_dir + os.sep + cancer_name
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        dump_data_into_tsv_according_to_cancer_type_and_stage(cancer_name, uuid_dict[cancer_name], out_dir, new_profile_list)


        #all_cancer_profiles.append(new_profile_list[0])
    #plot_mean_and_var(all_cancer_profiles)