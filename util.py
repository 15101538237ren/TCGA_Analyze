# -*- coding:utf-8 -*-
import os, math, re, json, time, pickle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

manifest_path = "data/manifest.txt"
tumor_suppressed_gene_file = "data/TSG.txt"
figure_path = "figure"
json_file_path = "data/metadata.json"
normal_keyword = "normal"
tumor_stages = ["normal","i","ia","ib","ii","iia","iib","iic","iii","iiia","iiib","iiic","iv","iva","x","not reported"]
#tumor_stages = ["normal","i","ia","iia","iib","iiia","iiib","iiic","x","not reported"]
tumor_stages_xaxis = {}

for idx, item in enumerate(tumor_stages):
    tumor_stages_xaxis[item] = idx + 1

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
    return genes

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

TSG = read_genenames(tumor_suppressed_gene_file)

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

                    if len(obj["cases"]):
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
def plot_for_each_gene(cancer_name, gene_name,x, y, box_data, y_means, c, xrange, xticks):
    plt.clf()
    plt.cla()
    fig, ax = plt.subplots()
    plt.xticks(xrange, xticks)
    ax.scatter(x, y, color=c, s=2.0)
    ax.boxplot(box_data,sym='')
    for xidx, y_mean in enumerate(y_means):
        if y_mean:
            ax.text(xidx + 0.8, y_mean, str(round(y_mean,2)),color='red',fontsize=12)
    ax.set_xlim([0, len(tumor_stages) - 0.5])
    ax.set_ylim([0, 1.0])
    ax.set_title(gene_name + " methylation for different cancer stage")
    figure_dir = figure_path+ os.sep + cancer_name
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)
    plt.savefig(figure_dir +os.sep + gene_name.lower() + '.png')

#checked
def gene_and_cancer_stage_profile_of_dna_methy(cancer_name, data_path, pickle_filepath,uuids, load=False):
    if not load:
        profile = {}
        x_profile = {}
        y_profile = {}
        for gene in TSG:
            profile[gene] = []
            x_profile[gene] = []
            y_profile[gene] = []
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
            while line:
                line_contents = line.split("\t")
                try:
                    gene_symbols = line_contents[5].split(";")
                    positions_to_tss = line_contents[8].split(";")
                    beta_val = -1.0 if line_contents[1] == "NA" else float(line_contents[1])
                    for idx, gene_symbol in enumerate(gene_symbols):
                        if (gene_symbol in TSG) and (-1500 <= int(positions_to_tss[idx]) <= 1000) and beta_val > 0.0:
                            profile[gene_symbol][tumor_stages_xaxis[uuid_to_stage[uuid]] - 1].append(beta_val)
                            x_profile[gene_symbol].append(tumor_stages_xaxis[uuid_to_stage[uuid]])
                            y_profile[gene_symbol].append(beta_val)
                            #one gene only add once for each cpg
                            break
                except IndexError, e:
                    print "line_contents :",
                    print line_contents

                line=now_file.readline()
            now_file.close()
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
        print "start load pickle file %s" % (pickle_filepath)
        profile = pickle.load(pickle_file)
        x_profile = pickle.load(pickle_file)
        y_profile = pickle.load(pickle_file)
        pickle_file.close()
        print "load pickle file %s finished" % (pickle_filepath)
    # for gene in TSG:
    #     xs = range(1,len(tumor_stages)+1)
    #     y_means = [np.array(arr).mean() for arr in profile[gene]]
    #     y_means = y_means[0:-1]
    #     plot_for_each_gene(cancer_name, gene, x_profile[gene], y_profile[gene], profile[gene], y_means,"blue",xs,tumor_stages)
if __name__ == '__main__':
    base_dir = "data"
    cancer_names = ["BRCA","COAD","LIHC","LUAD","LUSC"]
    for cancer_name in cancer_names:
        print "now analysing %s" % cancer_name
        data_path = base_dir + os.sep+ cancer_name + os.sep
        pickle_filepath = base_dir + os.sep + cancer_name + ".pkl"
        filenames = os.listdir(data_path)
        uuids = get_exist_uuids_from_filenames(filenames)
        gene_and_cancer_stage_profile_of_dna_methy(cancer_name,data_path, pickle_filepath,uuids, load=False)
