library(hash)
library(fitdistrplus)
library(dplyr)

gene_names = read.table("gene_with_protein_product_new.tsv", header=FALSE, stringsAsFactors = FALSE) #GeneNames类型是list,使用mode函数查看
#print(head(gene_names))
#headGe = head(gene_names)
#print(mode(headGe[[1]]))
#print(mode(gene_names))
#print(mode(gene_names[[1]][1]))
#print(head(gene_names))
#print(length(gene_names))

sig_level = 0.01

onco_names = read.table("Oncogene.txt", header=FALSE, stringsAsFactors = FALSE)
tsg_names = read.table("TSG.txt", header=FALSE, stringsAsFactors = FALSE)

onco_hash = hash(keys=onco_names[[1]], values=rep(1, length(onco_names[[1]])))
tsg_hash = hash(keys=tsg_names[[1]], values=rep(2, length(tsg_names[[1]])))
#print(head(onco_names[[1]]))
gene_category = rep(3, length(gene_names[[1]]))

for(i in 1:length(gene_names[[1]]))
{
    gene_name = gene_names[[1]][i]
    #print(gene_name)
    if(has.key(gene_name, onco_hash))
    {
        #print(onco_hash[[gene_name]]) 
        gene_category[i] = onco_hash[[gene_name]]
    }
    else if(has.key(gene_name, tsg_hash))
    {
        #print(tsg_hash[[gene_name]])
        gene_category[i] = tsg_hash[[gene_name]]
    }

}
#print(gene_category)
#print(length(gene_category))  len = 19114

methy_data_dir = "./methy_data"

gene_num = length(gene_names[[1]])
cancer_category_num = 5
file_name_info_list = c("gene", "stage", "cancer")
stage_name_list = c("normal", "i")
cancer_name_list = c("BRCA", "COAD", "LIHC", "LUAD", "LUSC")
gene_category_list = c("Onco", "TSG", "Other")

get_info_idx = function(info)
{
    if(info == "gene")
        return (1)
    else if(info == "stage")
	    return (2)
    else if(info == "cancer")
        return (3)
}

get_category_idx = function(category)
{
    if(category == "Onco")
        return (1)
    else if(category == "TSG")
        return (2)
    else if(category == "Other")
        return (3)
}

get_cancer_idx = function(cancer_name)
{
    if(cancer_name == "BRCA")
        return (1)
    else if(cancer_name == "COAD")
        return (2)
    else if(cancer_name == "LIHC")
        return (3)
    else if(cancer_name == "LUAD")
        return (4)
    else if(cancer_name == "LUSC")
        return (5)
}

get_info_list = function(file_name)
{
    str_list = strsplit(file_name, "[.]")
    #print(str_list[[1]])
    file_name_pre = str_list[[1]][1]
    #print(file_name_pre)
    info_list = strsplit(file_name_pre, "_")
    #print(info_list)
    return (info_list[[1]])
}

get_gene_name = function(file_name)
{
    info_list = get_info_list(file_name)
    gene_idx = get_info_idx("gene")
    gene_name = info_list[gene_idx]
    return (gene_name)
}

get_stage_name = function(file_name)
{
    info_list = get_info_list(file_name)
    stage_idx = get_info_idx("stage")
    stage_name = info_list[stage_idx]
    return (stage_name)
}

get_cancer_name = function(file_name)
{
    info_list = get_info_list(file_name)
    cancer_idx = get_info_idx("cancer")
    cancer_name = info_list[cancer_idx]
    return (cancer_name)
}

beta_hypothesis_test = function(x, shape1, shape2)
{
    beta_test_res = ks.test(x, "pbeta", shape1, shape2)
    #print(beta_test_res)
    return (beta_test_res$p.value)
}

transfer_cancer_idx_to_name = function(series)
{
    cancer_idx = as.numeric(series["CancerNames"])
    #print(cancer_idx)
    return (cancer_name_list[cancer_idx])
}

transfer_category_idx_to_name = function(series)
{
    category_idx = as.numeric(series["GeneCategory"])
    return (gene_category_list[category_idx])
}


gene_names_in_frame = rep(gene_names[[1]], times=cancer_category_num)
cancer_names_in_frame = sort(rep(seq(1,length(cancer_name_list)), times=gene_num))
gene_category_in_frame = rep(gene_category, times=cancer_category_num)
#print(gene_names_in_frame)
#print(mode(cancer_names_in_frame))
#print(cancer_names_in_frame)
#print(mode(gene_category_in_frame))
data_frame_len = length(gene_names_in_frame)

#print(mode(gene_names_in_frame))

cancer_dataframe = data.frame(GeneNames=gene_names_in_frame, CancerNames=cancer_names_in_frame,SubNums=rep(0,data_frame_len), 
    TotalNums=rep(0,data_frame_len), Ratio=rep(0,data_frame_len), GeneCategory=gene_category_in_frame)

cancer_dataframe$GeneNames = as.character(cancer_dataframe$GeneNames)
#print(head(cancer_dataframe))
#print(cancer_dataframe["GeneNames"])
#   print(str(cancer_dataframe))


for(i in 1:nrow(cancer_dataframe))
{
    gene_name = cancer_dataframe[i, "GeneNames"]
    cancer_name_idx = cancer_dataframe[i, "CancerNames"]
    cancer_name = cancer_name_list[cancer_name_idx]
    noraml_file_name = sprintf("%s_%s_%s.dat", gene_name, stage_name_list[1], cancer_name) #normal idx = 1
    i_th_file_name = sprintf("%s_%s_%s.dat", gene_name, stage_name_list[2], cancer_name) #i_th idx = 2
    #print(noraml_file_name)
    #print(i_th_file_name)
    if(file.exists(methy_data_dir))
    {
        normal_file_path = sprintf("%s/%s", methy_data_dir, noraml_file_name)
        i_th_file_path = sprintf("%s/%s", methy_data_dir, i_th_file_name)
        if(file.exists(normal_file_path) && file.exists(i_th_file_path))
        {
            normal_methy = read.table(normal_file_path)[[1]]
            i_th_methy = read.table(i_th_file_path)[[1]]
            
            if(sum(is.na(normal_methy)) > 0 || sum(is.na(i_th_methy)) > 0)
            {
                next
            }
            #print(head(i_th_methy))
            #print(length(i_th_methy))
            fit_res = fitdist(normal_methy, "beta")
            fit_estimate_res = fit_res$estimate
            fit_shape1 = fit_estimate_res["shape1"] #p
            fit_shape2 = fit_estimate_res["shape2"] #q
            #plot(fit_res)
            #print(fit_res)

            i_th_array = array(i_th_methy, c(length(i_th_methy), 1))
            #print(i_th_array)
            p_val_list = apply(i_th_array, 1, beta_hypothesis_test, shape1=fit_shape1, shape2=fit_shape2)
            #print(length(p_val_list))
            sub_num = sum(p_val_list >= sig_level)
            #print(sub_num)
            cancer_dataframe[i, "SubNums"] = sub_num
            cancer_dataframe[i, "TotalNums"] = length(p_val_list)
            cancer_dataframe[i, "Ratio"] = sub_num / length(p_val_list)
            if(i %% 100 == 0)
            {
                cat("i = ", i, "is finished\n")
            }
        }
    }
    
}

print("YES")

cancer_dataframe = arrange(cancer_dataframe, GeneCategory, desc(Ratio))

for(i in 1:length(cancer_name_list))
{
    sub_dataframe = cancer_dataframe[which(cancer_dataframe$CancerNames==i), ]
    sub_dataframe$CancerNames = apply(sub_dataframe, 1, transfer_cancer_idx_to_name)
    sub_dataframe$GeneCategory = apply(sub_dataframe, 1, transfer_category_idx_to_name)
    output_file_name = sprintf("%s.csv",cancer_name_list[i])
    #output_dataframe = sub_dataframe[c("GeneNames", "CancerNames", "SubNums", "TotalNums", "Ratio", "GeneCategory")]
    write.csv(output_dataframe, file=output_file_name, row.names=F, quote=F)
}
