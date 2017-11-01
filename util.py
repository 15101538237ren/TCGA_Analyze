# -*- coding:utf-8 -*-
import os, math, re, json
import numpy as np

manifest_path = "data/manifest.txt"
data_path = "data/BRCA"
json_file_path = "data/metadata.json"

normal_keyword = "normal"
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

def get_exist_uuids_from_filenames(filenames):
    [_, filename_to_uuid] = connect_filename_to_uuid()

    uuids = []

    for filename in filenames:
        uuids.append(filename_to_uuid[filename])

    print "get_exist_uuids_from_filenames called"
    return uuids
def connect_uuid_to_cancer_stage(uuids, json_file_path):
    [uuid_to_filename, _] = connect_filename_to_uuid()
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
if __name__ == '__main__':
    filenames = os.listdir(data_path)
    uuids = get_exist_uuids_from_filenames(filenames)
    connect_uuid_to_cancer_stage(uuids, json_file_path)