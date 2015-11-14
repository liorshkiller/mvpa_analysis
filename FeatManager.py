#!/usr/bin/python

from mako.template import Template
import os
import tempfile
import subprocess
from glob import glob

class FeatManager(object):
	def __init__(self, op, fsf_file):
		self._op = op
		self._fsf_file = os.path.join(op.study_dir(), "models", fsf_file)
		self.__read_feat_template__(self._fsf_file)	
	
	def __read_feat_template__(self, fsf_file):
		with open(fsf_file, 'r') as fh:
			template_data = fh.read()

		self._fsf_template = Template(template_data)
	
	def __generate_feat_config__(self, run_name, subject_code, config_file):
		base_dir = self._op.study_dir()
		sub_name = "sub{:0>3d}".format(int(subject_code))
		run_name = run_name

		fsf_config = self._fsf_template.render(base_dir=base_dir, sub_name=sub_name, run_name=run_name)

		with open(config_file, 'w') as config_file_handler:
				config_file_handler.write(fsf_config)
	
	def __execute_feat__(self, subject_code, run_name):
		if os.path.exists(os.path.join(self._op.study_dir(), "sub{:0>3d}".format(int(subject_code)), 'model', 'model001', '{}.feat'.format(run_name))):
			return

		config_file = os.path.join(self._op.study_dir(), "sub{:0>3d}".format(int(subject_code)), 'model', 'model001', '{}.fsf'.format(run_name)) 
		self.__generate_feat_config__(run_name, subject_code, config_file)

		subprocess.call('feat {}'.format(config_file), shell=True)

	def __get_task_runs__(self, subject_code, task_number):
		sub_name  = "sub{:0>3d}".format(int(subject_code))
		task_name = "task{:0>3d}".format(int(task_number))

		task_runs = glob(os.path.join(self._op.study_dir(), sub_name, 'BOLD', '{}_*'.format(task_name)))
		run_names = [task_run.split('/')[-1] for task_run in task_runs]

		return run_names

	def execute(self, subject, task):
		for run_name in self.__get_task_runs__(subject, task):
			self.__execute_feat__(subject, run_name)


if __name__ == "__main__":
	import OpenFMRIData
	from threading import Thread
	op = OpenFMRIData.OpenFMRIData('/home/user/data', 'LP')

	FM1= FeatManager(op, 'model001/task001.fsf')
	FM2= FeatManager(op, 'model001/task002.fsf')
	FM_gfeat= FeatManager(op, 'model001/task001_gfeat.fsf')


	#x = FM.__execute_feat__('task001_run002', 2)
	#x = FM.__execute_feat__('task001_run001', 2)
	#FM.__get_task_runs__(1, 1)
	all_subject_dirs = op.all_subjects_dirs_with_raw()
	for subject in all_subject_dirs:
		#Thread(target=FM2.execute, args=(subject.subcode(),2)).start()
		#Thread(target=FM1.execute, args=(subject.subcode(),1)).start()
		FM_gfeat.__execute_feat__(subject.subcode(),'task001')
