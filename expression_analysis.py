# -*- coding:utf-8 -*-
import os, json

from util import read_whole_genenames
cancer_names = ["BRCA", "COAD", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "THCA"]

tumor_stage_convert = {"normal":"normal","i":"i","ia":"i","ib":"i","ii":"ii","iia":"ii","iib":"ii","iic":"ii","iii":"iii","iiia":"iii","iiib":"iii","iiic":"iii","iv":"iv","iva":"iv","ivb":"iv","ivc":"iv","x":"x","not reported":"not reported"}
merged_stage = ["normal","i","ii","iii","iv","x","not reported"]
merged_stage_n = ["i","ii","iii","iv","x"]

rna_base_dir = "rna"
rna_out_dir = "rna_dat"
run_needed = "run_needed"
metadata_dir = os.path.join(rna_base_dir, "metadata")
rna_data_dir =  os.path.join(rna_base_dir, "rna_expression")
htsep_file_end = ".htseq.counts"

tumor_suppressed_gene_file = os.path.join(run_needed, "gene_with_protein_product.tsv")
[TSG, alias_dict] = read_whole_genenames(tumor_suppressed_gene_file)

gene_idx_path = os.path.join(rna_out_dir, "gene_idx.txt")
if not os.path.exists(gene_idx_path):
    with open(gene_idx_path,"w") as gene_idx_file:
        gene_idx_file.write("\n".join([str(gidx+1) + "\t" + gene for gidx, gene in enumerate(TSG)]))

#获得每种癌症的htseq-count数据的文件名列表
def obtain_htseq_count_filelist():
    filename_dict = {}
    for cancer_name in cancer_names:
        cancer_data_dir = os.path.join(rna_data_dir, cancer_name)

        output_cancer_dir = os.path.join(rna_out_dir, cancer_name)
        if not os.path.exists(output_cancer_dir):
            os.makedirs(output_cancer_dir)

        htseq_file_names = []
        filenames = os.listdir(cancer_data_dir)
        for filename in filenames:
            if filename.endswith(htsep_file_end):
                htseq_file_names.append(filename)
        filename_dict[cancer_name] = htseq_file_names

        outfile_path = os.path.join(output_cancer_dir, cancer_name + "_htseq_filelist.txt")
        write_tab_seperated_file_for_a_list(outfile_path, htseq_file_names, index_included=True)

    print "obtain htseq_count filelist successful!"
    return [filename_dict]

#向target_file_path中写target_list的值, 如果index_included=True,则第一列为自增索引, 第二列为value
def write_tab_seperated_file_for_a_list(target_file_path, target_list, index_included=True):
    with open(target_file_path,"w") as outfile:
        if index_included:
            outfile.write("\n".join([str(gidx+1) + "\t" + str(value) for gidx, value in enumerate(target_list)]))
        print "write %s successful" % target_file_path

#读取tab分隔的文件(input_file_path) 的第target_col_index, 返回该列的所有值到一个list
def read_tab_seperated_file_and_get_target_column(target_col_index, input_file_path):
    ret_value_list = []
    with open(input_file_path, "r") as input_file:
        line = input_file.readline()
        while line:
            value = line.split("\t")[target_col_index][0:-1]
            ret_value_list.append(value)
            line = input_file.readline()
    return ret_value_list

#
def obtain_stage_info_from_metadata():
    #make output_dir
    for cancer_name in cancer_names:
        output_cancer_dir = os.path.join(rna_out_dir, cancer_name)
        if not os.path.exists(output_cancer_dir):
            os.makedirs(output_cancer_dir)

        input_path = os.path.join(output_cancer_dir, cancer_name + "_htseq_filelist.txt")
        output_path = os.path.join(output_cancer_dir, cancer_name + "_stages.txt")
        metadata_file_path = os.path.join(metadata_dir, "rna_" + cancer_name + "_metadata.json")
        htseq_file_names = read_tab_seperated_file_and_get_target_column(1, input_path)
        file_name_to_stage = {}
        json_obj = json.load(open(metadata_file_path,'r'))
        stage_list = []
        for htseq_file_name in htseq_file_names:
            htseq_file_name_extend = htseq_file_name + ".gz"
            stage_save = "not reported"
            for obj in json_obj:
                if obj["file_name"] != htseq_file_name_extend:
                    continue
                else:
                    entity_submitter_id = obj["associated_entities"][0]["entity_submitter_id"]
                    sample_type = entity_submitter_id.split("-")[3]
                    tumor_type = int(sample_type[0 : -1])
                    normal = 1 if tumor_type > 9 else 0

                    if normal:
                        stage_save = "normal"
                    else:
                        if "cases" in obj.keys():
                            if len(obj["cases"]):
                                if "diagnoses" in obj["cases"][0].keys():
                                    if len(obj["cases"][0]["diagnoses"]):
                                        stage = obj["cases"][0]["diagnoses"][0]["tumor_stage"]
                                        if stage != "not reported":
                                            stage_save = tumor_stage_convert[stage.split(" ")[1]]
                    break
            file_name_to_stage[htseq_file_name] = stage_save
            stage_list.append(stage_save)
        write_tab_seperated_file_for_a_list(output_path, stage_list, index_included=True)
if __name__ == '__main__':
    #obtain_htseq_count_filelist()
    obtain_stage_info_from_metadata()