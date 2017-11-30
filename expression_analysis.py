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

# 获得每种癌症的htseq-count数据的文件名列表
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

# 向target_file_path中写target_list的值, 如果index_included=True,则第一列为自增索引, 第二列为value
def write_tab_seperated_file_for_a_list(target_file_path, target_list, index_included=True):
    with open(target_file_path,"w") as outfile:
        if index_included:
            outfile.write("\n".join([str(gidx+1) + "\t" + str(value) for gidx, value in enumerate(target_list)]))
        else:
            outfile.write("\n".join([str(value) for value in target_list]))
        print "write %s successful" % target_file_path

# 读取tab分隔的文件(input_file_path) 的第target_col_index, 返回该列的所有值到一个list
def read_tab_seperated_file_and_get_target_column(target_col_index, input_file_path):
    ret_value_list = []
    with open(input_file_path, "r") as input_file:
        line = input_file.readline()
        while line:
            line_contents = line.split("\t")
            if target_col_index== len(line_contents)-1:
                value = line_contents[target_col_index][0:-1]
            else:
                value = line_contents[target_col_index]
            ret_value_list.append(value)
            line = input_file.readline()
    return ret_value_list

# 从meta_data文件中获取所有htseq_count文件对应的癌症主分期列表,并将其输出到文件cancer_name_stages.txt中
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

# 所有htseq-count文件的第一列的Ensembl基因列表顺序和数量均相同,要注意最后5行为统计信息,不是基因名
def generate_ensembl_gene_ids_in_a_htseq_file():
    first_htseq_file_path = ""
    for cancer_name in cancer_names:
        cancer_data_dir = os.path.join(rna_data_dir, cancer_name)
        filenames = os.listdir(cancer_data_dir)
        flag = False
        for filename in filenames:
            if filename.endswith(htsep_file_end):
                flag = True
                first_htseq_file_path = os.path.join(cancer_data_dir, filename)
                break
        if flag:
            break

    ensembl_gene_ids = read_tab_seperated_file_and_get_target_column(0, first_htseq_file_path)
    outfile_path = os.path.join(rna_out_dir, "ensembl_gene_ids.txt")
    write_tab_seperated_file_for_a_list(outfile_path, ensembl_gene_ids[0:-5], index_included=False)

# 根据基因列表文件和gtf注释文件,将gene_id与gene_name进行match,之后输出到gene_symbols.txt中,该文件的索引与gene_id的索引一一对应
def obtain_gene_symbols_and_emsembl_ids_from_gtf(gene_ids_filepath, gtf_filepath, outfile_path):
    ensembl_gene_ids = read_tab_seperated_file_and_get_target_column(0, gene_ids_filepath)
    ensembl_gene_id_to_gene_symbol_dict = {}

    gtf_file = open(gtf_filepath , "r")
    for i in range(5):
        gtf_file.readline()
    line = gtf_file.readline()

    while line:
        line_contents = line.split("\t")
        feature = line_contents[2]
        if feature == "gene":
            groups = line_contents[-1].split(";")[0:-1]
            group_dict = {}
            for item in groups:
                [k , v] = item.strip().split(" ")
                group_dict[k] = v.replace('\"',"")
            if "gene_type" in group_dict.keys() and group_dict["gene_type"] == "protein_coding":
                if ("gene_id" in group_dict.keys()) and ("gene_name" in group_dict.keys()):
                    ensembl_gene_id_to_gene_symbol_dict[group_dict["gene_id"]] = group_dict["gene_name"]
        line = gtf_file.readline()
    gene_symbols = []
    for ensembl_gene_id in ensembl_gene_ids:
        if ensembl_gene_id in ensembl_gene_id_to_gene_symbol_dict.keys():
            gene_symbol = ensembl_gene_id_to_gene_symbol_dict[ensembl_gene_id]
        else:
            gene_symbol = "-"
        gene_symbols.append(gene_symbol)
    print "finished handling gtf file, now writing gene_symbols to %s" % outfile_path
    write_tab_seperated_file_for_a_list(outfile_path, gene_symbols, index_included=True)

def connect_whole_genome_to_ensembl_gene_symbol_index(whole_gene_idx_filepath , ensembl_gene_symbol_index_path, out_correspondent_index_path):
    whole_genome = read_tab_seperated_file_and_get_target_column(1, whole_gene_idx_filepath)
    ensembl_gene_symbols = read_tab_seperated_file_and_get_target_column(1, ensembl_gene_symbol_index_path)

    #利用基因名查在全基因组中的index
    whole_gene_index_dict = {gene_name: (gindex + 1) for gindex, gene_name in enumerate(whole_genome)}

    #给出全基因组中所有基因对应ensembl_gene_symbols文件中该基因的行号, -1代表没有该基因对应的gene_symbol
    ret_index_dict = {(gindex + 1): -1 for gindex, gene_name in enumerate(whole_genome)}

    for gindex2, gene_symbol in enumerate(ensembl_gene_symbols):
        if (gene_symbol != "-") and (gene_symbol in alias_dict.keys()):
            gene_name = alias_dict[gene_symbol]
            if gene_name in whole_gene_index_dict.keys():
                origin_gene_idx_in_whole_genome_file = whole_gene_index_dict[gene_name]
                ret_index_dict[origin_gene_idx_in_whole_genome_file] = gindex2
    ordered_dict= sorted(ret_index_dict.items(), key=lambda d:d[0])

    with open(out_correspondent_index_path, "w") as out_correspondent_index_file:
        ltw = []
        for (origin_index, gene_symbol_index) in ordered_dict:
            ltw.append("\t".join([str(origin_index), str(gene_symbol_index)]))
        out_correspondent_index_file.write("\n".join(ltw))
    print "connect_whole_genome_to_ensembl_gene_symbol_index successful"
#some global variables
gene_idx_path = os.path.join(rna_out_dir, "gene_idx.txt")
if not os.path.exists(gene_idx_path):
    with open(gene_idx_path,"w") as gene_idx_file:
        gene_idx_file.write("\n".join([str(gidx+1) + "\t" + gene for gidx, gene in enumerate(TSG)]))

gene_symbols_filepath = os.path.join(rna_out_dir, "gene_symbols.txt")
emsembl_ids_filepath = os.path.join(rna_out_dir, "ensembl_gene_ids.txt")
gtf_filepath = os.path.join("data","GRCh38","gencode.v22.annotation.gtf")
if not os.path.exists(gene_symbols_filepath):
    obtain_gene_symbols_and_emsembl_ids_from_gtf(emsembl_ids_filepath,gtf_filepath, gene_symbols_filepath)

correspondent_index_path = os.path.join(rna_out_dir, "correspondent_index.txt")
if not os.path.exists(correspondent_index_path):
    connect_whole_genome_to_ensembl_gene_symbol_index(gene_idx_path, gene_symbols_filepath,correspondent_index_path)
if __name__ == '__main__':
    pass

