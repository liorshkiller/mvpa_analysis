#!/usr/bin/python
import sys
import os

from OpenFMRIData import OpenFMRIData
from OpenFMRIAnalyzer import OpenFMRIAnalyzer

subnames    	= sys.argv[1:]
data_dir   	= os.environ.get('DATA_DIR') or '/home/user/data'
study_name 	= os.environ.get('STUDY_NAME') or 'LP'


op = OpenFMRIData(data_dir, study_name)
subcodes = []
for subname in subnames:
	subject_dir = op.subject_dir(subname=subname)
	subcodes.append(subject_dir.subcode())
analyzer = OpenFMRIAnalyzer(op,subcodes)
analyzer.analyze(mc_align_task=True, func_seg=True)

