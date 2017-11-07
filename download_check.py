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