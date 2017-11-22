

cancer_markers = {"BRCA":'rs', "COAD":'gp', "LIHC":'bo',"LUAD":'kx',"LUSC":'c*'}

#checked



#checked

#checked



#checked


#checked

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





if __name__ == '__main__':
    dna_dir = data_dir + os.sep + "GRCh38"
    file_pre = "Homo_sapiens.GRCh38.dna.chromosome."
    chr_list = range(1, 23)
    sequence_rtn = get_all_dna_sequences(dna_dir, file_pre, chr_list)

    all_cancer_profiles = []
    stat_names = ["mean","std"]
    stat_epsilons = [0.05, 0.05]
    out_stage_list = ["normal","i"]

    for cancer_name in cancer_names:
        if cancer_name in ["LGG","OV","GBM","LAML", "PRAD","UCEC","SARC", "UVM","CESC", "DLBC"]:
            continue
        print "now start %s" % cancer_name
        data_path = data_dir + os.sep+ cancer_name + os.sep
        pickle_filepath = pickle_dir + os.sep + cancer_name + ".pkl"
        #local scripts
        # filenames = os.listdir(data_path)
        # uuids = get_exist_uuids_from_filenames(filenames)


        #server script
        temp_profile_list = gene_and_cancer_stage_profile_of_dna_methy(cancer_name,data_path, pickle_filepath,uuid_dict[cancer_name], load=True, whole_genes= True)
        # save_cancer_std_and_mean_of_all_genes(cancer_name, temp_profile_list, ["normal", "i"],out_dir="mean_std_data")
    # cmp_gene_variations_in_mean_and_std(cancer_names, stat_names, stat_epsilons, "normal", "i", input_dir ="mean_std_data", out_dir="stat")




        #all_cancer_profiles.append(new_profile_list[0])
    #plot_mean_and_var(all_cancer_profiles)