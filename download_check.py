import os, re, shutil

#check the stop line from the download dir and manifest file
def check_download_interrupt(manifest, download_dir):
    now_file = open(manifest,'r')
    now_file.readline()
    line = now_file.readline()
    lineno = 1
    while line:
        lineno += 1
        line_content = line.split("\t")
        file_path = line_content[1]
        if not os.path.exists(download_dir + os.sep + file_path):
            print "line %d, file %s" % (lineno, file_path)
            break
        else:
            line = now_file.readline()
    now_file.close()
def exclude_data_list_from_path(manifest_path, base_dir, outfile_path):
    now_file = open(manifest_path,'r')
    now_file.readline()
    outfile = open(outfile_path,"w")

    str_pattern = r'([^\t]+)\t([^\t]+)'
    cancer_pattern = r'jhu-usc.edu_([^\.]+)*'

    line = now_file.readline()
    while line:
        match_p = re.search(str_pattern, line)
        if match_p:
            uuid = match_p.group(1)
            file_name = match_p.group(2)
            cancer_name = re.search(cancer_pattern, file_name).group(1)
            if not os.path.exists(base_dir + os.sep + cancer_name + os.sep + file_name):
                outfile.write(line)
        line = now_file.readline()
    now_file.close()
    outfile.close()
    print "exclude_data_list_from_path success!"

def exclude_data_list_from_file(manifest_path, data_list_path, outfile_path):
    d_file = open(data_list_path, "r")
    files_list = d_file.readlines()
    existed_file_names = [file_name[0:-1] for file_name in files_list]
    outfile = open(outfile_path,"w")

    m_file = open(manifest_path,'r')
    line = m_file.readline()
    outfile.write(line)

    line = m_file.readline()
    while line:
        line_content = line.split("\t")
        file_path = line_content[1]
        if file_path not in existed_file_names:
            outfile.write(line)
        line = m_file.readline()
    m_file.close()
    outfile.close()
    print "exclude_data_list_from_file successful"
def move_files_according_to_cancer_type(src_dir, dst_dir):
    filenames = os.listdir(src_dir)
    cancer_pattern = r'jhu-usc.edu_([^\.]+)*'
    for filename in filenames:
        file_path = src_dir + os.sep + filename
        cancer_name = re.search(cancer_pattern, filename).group(1)
        dst_subdir = dst_dir + os.sep + cancer_name
        if not os.path.exists(dst_subdir):
            os.makedirs(dst_subdir)
        if not os.path.exists(dst_subdir + os.sep + filename):
            shutil.move(file_path, dst_subdir)
    print "moved successfully!"
if __name__ == '__main__':
    #check_download_interrupt("10_caner_manifests.tsv", "downloads")
    # exclude_data_list_from_file("9_type_cancer_1116.tsv","data_list.txt","9_type_cancer_1116_new.tsv")
    # exclude_data_list_from_path("run_needed/10_cancer_manifest.tsv", "data", "10_cancer_manifest_new.tsv")
    move_files_according_to_cancer_type("remind_files","data")