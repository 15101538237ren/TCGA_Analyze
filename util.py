# -*- coding:utf-8 -*-
import os, math, re, json
import numpy as np

manifest_path = "data/manifest.txt"
tumor_suppressed_gene_file = "data/TSG.txt"
data_path = "data/BRCA/"
json_file_path = "data/metadata.json"
normal_keyword = "normal"
tumor_stages = ["i","ia","ib","ii","iia","iib","iic","iii","iiia","iiib","iiic","iv","x","not reported"]
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

def get_exist_uuids_from_filenames(filenames):
    uuids = []

    for filename in filenames:
        uuids.append(filename_to_uuid[filename])

    print "get_exist_uuids_from_filenames called"
    return uuids
def connect_uuid_to_cancer_stage(uuids, json_file_path):
    stage_to_uuids = {}
    uuid_to_stage = {}

    json_obj = json.load(open(json_file_path,'r'))

    for uuid in uuids:
        filename = uuid_to_filename[uuid]
        sample_type = filename.split(".")[5].split("-")[3]
        tumor_type = int(sample_type[0 : -1])
        normal = 1 if tumor_type > 9 else 0

        if not normal:
            # if cancer, detail stage classification
            pass
        else:
            uuid_to_stage[uuid] = normal_keyword
            if normal_keyword not in stage_to_uuids.keys():
                stage_to_uuids[normal_keyword] = [uuid]
            else:
                stage_to_uuids[normal_keyword].append(uuid)
        print normal
    return [uuid_to_stage, stage_to_uuids]
def gene_and_cancer_stage_profile_of_dna_methy(uuids):
    profile = {}
    for gene in TSG:
        profile[gene] = {}
        for tumor_stage in tumor_stages:
            profile[gene][tumor_stage] = []
    [uuid_to_stage, stage_to_uuids] = connect_uuid_to_cancer_stage(uuids, json_file_path)
    for uuid in uuids:
        file_path = uuid_to_filename[uuid]
        now_file = open(data_path + file_path,'r')
        now_file.readline()
        line = now_file.readline()
        while line:
            line_contents = line.split("\t")
            gene_symbols = line_contents[5].split(";")
            positions_to_tss = line_contents[8].split(";")
            beta_val = -1.0 if line_contents[1] == "NA" else float(line_contents[1])
            for idx, gene_symbol in enumerate(gene_symbols):
                if (gene_symbol in TSG) and (-1500 <= int(positions_to_tss[idx]) <= 1000) and beta_val > 0.0:
                    profile[gene_symbol][uuid_to_stage[uuid]].append(beta_val)
                    #one gene only add once for each cpg
                    break
            line=now_file.readline()
        now_file.close()
if __name__ == '__main__':
    filenames = os.listdir(data_path)
    uuids = get_exist_uuids_from_filenames(filenames)
    #connect_uuid_to_cancer_stage(uuids, json_file_path)
    #gene_and_cancer_stage_profile_of_dna_methy(uuids)
