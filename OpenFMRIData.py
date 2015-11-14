#!/usr/bin/python

import os
import json

from glob import glob
from SubjectDir import SubjectDir

class OpenFMRIData(object):
    def __init__(self, data_dir, study_name):
        self._data_dir                = data_dir
        self._raw_data_dir            = os.path.join(data_dir, 'raw', study_name)
        self._behavioural_dir	      = os.path.join(data_dir, 'behavioural', study_name)

        self._study_name              = study_name
        self._study_dir               = os.path.join(self._data_dir, self._study_name)

        self._subject_mapping_file    = os.path.join(self._study_dir, 'mapping_subject.json')

        self._task_order              = None
        self._subject_mapping         = None
        self._task_mapping            = None

        self.load_subject_mapping()
        self.load_task_order()
        self.load_behavioural_mapping()

    def load_behavioural_mapping(self):
        behaviuoral_mapping_file = os.path.join(self._behavioural_dir, 'task_mapping.txt')
        if os.path.exists(behaviuoral_mapping_file):
            with open(behaviuoral_mapping_file, 'r') as fh:
                task_data = fh.read().splitlines()

            self._task_mapping = dict()

            task_pairs = [(task.split('\t')[0], task.split('\t')[1]) for task in task_data]
            for task_name, task_sequence in task_pairs:
                self._task_mapping[task_name] = task_sequence

    def load_task_order(self):
        task_order_file = os.path.join(self._study_dir, 'task_order.txt')
        if os.path.exists(task_order_file):
            with open(task_order_file, 'r') as fh:
                self._task_order = fh.read().splitlines()

    def load_subject_mapping(self):
        if os.path.isfile(self._subject_mapping_file):
            mapping_file = open(self._subject_mapping_file, 'r')

            self._subject_mapping = json.load(mapping_file)
            mapping_file.close()
        else:
            self._subject_mapping = dict()
            self.__write_subject_mapping__()

    def all_subjects_list(self):
        if len(self.mapping_json().values()) > 0:
            return self.mapping_json().values()
        else:
            return map(lambda x: int(x.split('sub')[1]),glob(os.path.join(self._study_dir,"sub*")))

    def all_subjects_dirs_with_raw(self):
        dirs = []
        for subname in self.mapping_json().keys():
            dirs.append(self.subject_dir(subname=subname))
        return dirs


    def mapping_json(self):
        return self._subject_mapping

    def data_dir(self):
        return self._data_dir

    def study_dir(self):
        return self._study_dir

    def behavioural_dir(self):
        return self._behavioural_dir

    def raw_study_dir(self):
        return self._raw_data_dir

    def __get_latest_subject_directory__(self):
        directory_prefix = "sub"

        study_dirlist = glob("{}/{}*".format(self._study_dir, directory_prefix))

        latest_dir_sequence = (max([int(x.replace("{}/{}".format(self._study_dir, directory_prefix), "")) for x in study_dirlist] + [0])+1)

        return latest_dir_sequence

    def __write_subject_mapping__(self):
        mapping_file = open(self._subject_mapping_file, 'w')
        json.dump(self._subject_mapping, mapping_file)
        mapping_file.close()

    def __subcode_to_dir_format__(self, code):
        return os.path.join(self.study_dir(), "sub{:0>3d}".format(int(code)))

    def subject_dir(self, **kwargs):
        """Expecting either subcode (i.e. 1), or subname (i.e. AzOr)"""

        subject_code 	= None
        raw_dir      	= None
        behavioural_dir = None

        if 'subcode' in kwargs:
            subject_code = kwargs['subcode']
        elif 'subname' in kwargs:
            # load mapping from pickle, if subject name does not exist, create and populate pickle
            subject_name = kwargs['subname']


                # Check whether subject raw directory exists before adding the mapping
	    raw_dirs = glob('{}/*{}*'.format(self.raw_study_dir(), subject_name))

            if len(raw_dirs) == 0:
                raise BaseException("No subject by the name of {}".format(subject_name))

            raw_dir = raw_dirs[0]

            behavioural_dirs = glob('{}/*{}*'.format(self.behavioural_dir(), subject_name))

            if len(behavioural_dirs) == 0:
                raise BaseException("No behavioural data for subject {}".format(subject_name))
            behavioural_dir = behavioural_dirs[0]

            if subject_name in self._subject_mapping:
                subject_code = self._subject_mapping[subject_name]
            else:
                subject_code = self.__get_latest_subject_directory__()

                self._subject_mapping[subject_name] = subject_code
	        self.__write_subject_mapping__()

        return SubjectDir(subject_code, self.__subcode_to_dir_format__(subject_code), raw_dir, behavioural_dir, self._task_order, self._task_mapping)

    def handle_bias_field_map(self):
	for name, code in self._subject_mapping:
	    self.subject_dir(subname=name)
		
def test():
    o = OpenFMRIData('/home/user/data', 'LP')

    print "{:<40s}{}".format("data_dir", o.data_dir())
    print "{:<40s}{}".format("study_dir", o.study_dir())
    print "{:<40s}{}".format("raw_study_dir", o.raw_study_dir())
    print "{:<40s}{}".format("behavioural_dir", o.behavioural_dir())

    subject_dir_by_name = o.subject_dir(subname='CoMi')
    print "{:<40s}{}".format("subject(subname='CoMi').path()", subject_dir_by_name.path())
    print "{:<40s}{}".format("subject(subname='CoMi').raw_path()", subject_dir_by_name.raw_path())
    print "{:<40s}{}".format("subject(subname='CoMi').behavioural_path()", subject_dir_by_name.behavioural_path())
    print "{:<40s}{}".format("subject(subname='CoMi').create_bias_fieldmap()", subject_dir_by_name.create_bias_fieldmap())

    subject_dir = o.subject_dir(subcode=1)
    print "{:<40s}{}".format("subject(subcode=1).path()", subject_dir.path())
    print "{:<40s}{}".format("subject(subcode=1).anatomical_dir()", subject_dir.anatomical_dir())

    print "{:<40s}{}".format("mapping_json", o.mapping_json())
    print "{:<40s}{}".format("__get_latest_subject_directory__", o.__get_latest_subject_directory__())

    print "{:<40s}{}".format("task_order", o._task_order)

if __name__ == "__main__":
    test()
