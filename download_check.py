import os

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
if __name__ == '__main__':
    #check_download_interrupt("10_caner_manifests.tsv", "downloads")
    exclude_data_list_from_file("10_caner_manifests.tsv","data_list.txt","10_caner_manifests_new.tsv")