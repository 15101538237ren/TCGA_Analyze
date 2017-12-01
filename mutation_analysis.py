# -*- coding:utf-8 -*-
import os, json, requests
from methylation_analysis import read_whole_genenames
from expression_analysis import read_tab_seperated_file_and_get_target_column,write_tab_seperated_file_for_a_list
cancer_names = ["BRCA", "COAD", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "THCA"]
tumor_stage_convert = {"normal":"normal","i":"i","ia":"i","ib":"i","ii":"ii","iia":"ii","iib":"ii","iic":"ii","iii":"iii","iiia":"iii","iiib":"iii","iiic":"iii","iv":"iv","iva":"iv","ivb":"iv","ivc":"iv","x":"x","not reported":"not reported"}
merged_stage = ["normal","i","ii","iii","iv","x","not reported"]
merged_stage_n = ["i","ii","iii","iv","x"]

snv_base_dir = "snv"
run_needed = "run_needed"
tumor_suppressed_gene_file = run_needed + os.sep + "gene_with_protein_product.tsv"
[TSG, alias_dict] = read_whole_genenames(tumor_suppressed_gene_file)


#输入TCGA的submitter_id列表, 和查询大小(最大300), 返回一个[dict], dict的key:value分别为submitter_id : stage_name
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

#从突变文件(.maf)获取所有癌症各自的submitter_id列表,并将列表写入到output_cancer_dir的cancer_name_submitter_ids.txt文件中, 文件含有两列,第一列为index, 第二列为submitter_id
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

            write_tab_seperated_file_for_a_list(outfile_path, submitter_dict[cancer_name].keys(), index_included= True)

#1. 读取get_maf_submitter_ids函数产生的*_submitter_ids.txt文件, 从而获得该癌症的所有submitter_id列表;
#2. 通过query_stage_of_an_submitter_id函数, 查询query_size大小的submitter_id列表对应的癌症阶段列表(stages)
#3. *_submitter_ids.txt文件中submitter_id列表对应的癌症阶段写到cancer_name_stages.txt中,第一列为idx,对应submitter_ids.txt的idx, 第二列为癌症阶段名称(stage_name)
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
            submitter_ids = read_tab_seperated_file_and_get_target_column(1, input_path)
            print cancer_name
            len_submitters = len(submitter_ids)
            print len_submitters
            for i in range(0, len_submitters, query_size):
                sub_submitter_ids = submitter_ids[i:i+query_size]
                stages_of_subids_dict = query_stage_of_an_submitter_id(sub_submitter_ids, query_size)
                for k,v in stages_of_subids_dict[0].items():
                    submitter_id_stage_dict[k] = tumor_stage_convert[v]
            print submitter_id_stage_dict

            with open(output_path,"w") as outfile:
                for gidx, submitter_id in enumerate(submitter_ids):
                    if submitter_id not in submitter_id_stage_dict.keys():
                        submitter_id_stage_dict[submitter_id] = "not reported"
                outfile.write("\n".join([str(gidx+1)  + "\t" + submitter_id_stage_dict[submitter_id] for gidx, submitter_id in enumerate(submitter_ids)])) #+ "\t" + submitter_id
                print "write %s successful" % input_path

#1. 从input_cancer_dir中, 输入cancer_name_submitter_ids.txt和cancer_name_stages.txt文件,获取该癌症所有的submitter_id列表和癌症阶段列表
#2. 建立字典对应关系,方便用submitter_id查对应癌症阶段, 即代码中的submitter_id_to_stage_dict
#3. 返回submitter_id_to_stage_dict, 和所有submitter_id列表(按从cancer_name_submitter_ids.txt中读到的顺序)
def input_submitter_id_and_its_stages(input_cancer_dir, cancer_name):
    submitter_id_to_stage_dict = {}
    submitter_id_input_filepath = os.path.join(input_cancer_dir, cancer_name + "_submitter_ids.txt")
    stages_input_filepath = os.path.join(input_cancer_dir, cancer_name + "_stages.txt")
    submitter_ids = []
    stages =[]
    with open(submitter_id_input_filepath,"r") as submitter_id_input_file:
        line = submitter_id_input_file.readline()
        while line:
            submitter_id = line.split("\t")[1]
            if submitter_id.endswith("\n"):
                submitter_id = submitter_id[0:-1]
            submitter_ids.append(submitter_id)
            line = submitter_id_input_file.readline()
    with open(stages_input_filepath,"r") as stages_input_file:
        line = stages_input_file.readline()
        while line:
            stage = line.split("\t")[1]
            if stage.endswith("\n"):
                stage = stage[0:-1]
            stages.append(stage)
            line = stages_input_file.readline()
    for sidx in range(len(submitter_ids)):
        submitter_id_to_stage_dict[submitter_ids[sidx]] = stages[sidx]
    return [submitter_id_to_stage_dict, submitter_ids]

# 运行生成突变数据的流水线: 首先生成gene_id文件,方便对数据文件的gene_id索引
# 从.maf突变数据中,统计会产生编码区mRNA突变或者翻译后氨基酸改变的基因突变数量
# 输出各个癌症阶段对应的基因突变数据, cancer_name_stage_mutation_count.dat, 第一列为gene_index, 第一行为submitter_id的index(对应见cancer_name_submitter_ids.txt文件), 内容为某基因在某样本中影响翻译的基因突变数量
def dna_mutation_data_transform_pipline():
    outdir ="snv_dat"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    gene_idx_path = os.path.join(outdir, "gene_idx.txt")
    if not os.path.exists(gene_idx_path):
        write_tab_seperated_file_for_a_list(gene_idx_path, TSG, index_included=True)

    colum_idxs = [0, 8, 15] #要提取maf文件的列编号
    mutation_expections = ["Missense_Mutation", "Translation_Start_Site", "Splice_Region", "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Frame_Shift_Del", "Frame_Shift_Ins"]
    for cancer_name in cancer_names:
        cancer_dir = os.path.join(snv_base_dir, cancer_name)
        #只对有突变数据的癌症做分析
        if os.path.exists(cancer_dir):
            data_dict = {stage:{} for stage in merged_stage}
            output_cancer_dir = os.path.join(outdir, cancer_name)
            if os.path.exists(output_cancer_dir):
                [submitter_id_to_stage_dict, submitter_ids] = input_submitter_id_and_its_stages(output_cancer_dir, cancer_name)

                for cancer_stage in merged_stage:
                    for submitter_id in submitter_ids:
                        data_dict[cancer_stage][submitter_id] = {gene: 0 for gene in TSG}

                file_names = os.listdir(cancer_dir)
                for file_name in file_names:
                    if not file_name.startswith("."):
                        file_path = os.path.join(cancer_dir,file_name)
                        snv_file = open(file_path, "r")
                        print "start %s" % file_path
                        # pass head 6 lines
                        for i in range(6):
                            snv_file.readline()
                        line = snv_file.readline()
                        while line:
                            line_contents = line.split("\t")
                            [gene_name, variation_classification, bar_code] = [line_contents[i] for i in colum_idxs]
                            submitter_id = "-".join(bar_code.split("-")[0:3])
                            stage = submitter_id_to_stage_dict[submitter_id]

                            if gene_name in TSG and (variation_classification in mutation_expections):
                                data_dict[stage][submitter_id][gene_name] += 1
                            line = snv_file.readline()
                        print "end %s" % file_path

                for cancer_stage in merged_stage_n:
                    with open(os.path.join(output_cancer_dir, cancer_name + "_" + cancer_stage+"_mutation_count.dat"), "w") as data_file:
                        header = "\t".join([str(item) for item in range(len(submitter_ids) + 1)]) + "\n"
                        data_file.write(header)
                        data_str = []
                        for gidx, gene in enumerate(TSG):
                            arr = [gidx + 1]
                            arr.extend([data_dict[cancer_stage][submitter_id][gene] for submitter_id in submitter_ids])
                            data_str.append("\t".join([str(item) for item in arr]))
                        data_file.write("\n".join(data_str))

if __name__ == '__main__':
    # get_maf_submitter_ids()
    # get_submitter_id_stages()
    # dna_mutation_data_transform_pipline()
    pass