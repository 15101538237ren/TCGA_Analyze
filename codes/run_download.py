#!/usr/bin/python
#encoding:utf-8
#this script download data from GDC Portal by using UUID and rename them.
#This script was inspired by [python download file - OPEN](http://www.open-open.com/lib/view/open1430982247726.html)

import urllib
import os

cancer_names = ["COAD", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "THCA"]

for cancer_name in cancer_names:
	infile = "../manifests/rna_" + cancer_name + "_manifest.txt"
	outdir = "../rna/" +cancer_name + "/"
	print "now %s" % cancer_name
	#outdir exists or not
	if(not os.path.exists(outdir)):
		print outdir,"not exists, create"
		os.makedirs(outdir)

	tsv = open(infile)
	tsv.readline()


	for line in tsv:
		arr = line.split('\t')
		UUID, filename = arr[0], arr[1]
		print filename, "is downloading..."
		url = os.path.join('https://gdc-api.nci.nih.gov/data/', UUID)
		local = os.path.join(outdir, filename)
		urllib.urlretrieve(url,local)
	tsv.close()
