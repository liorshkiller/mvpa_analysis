#!/usr/bin/python

import os
from glob2 import glob
import shutil
import pandas as pd


class SubjectDir(object):
    def __init__(self, subject_code, path, raw_path=None, behavioural_path=None, task_order=None, task_mapping=None):
        self._path = path
        self._raw_path = raw_path
        self._behavioural_path = behavioural_path
        self._task_order = task_order
        self._task_mapping = task_mapping
        self._subject_code = subject_code

        self._subdirs = {'functional': 'BOLD',
                         'anatomical': 'anatomy',
                         'model': 'model',
                         'masks': 'masks',
                         'fieldmap': 'fieldmap'}

        if (not os.path.isdir(self._path)):
            if raw_path is None:
                raise Exception("Cannot create new subject directory from subcode")
            else:
                self.__create_subject_dirtree__()

        self.__load_dirtree__()

    def __load_dirtree__(self):
        self._dir_tree = {'functional': dict(),
                          'model': dict()}

        task_pairs = [(directory.split('/')[-1].split('_')[0], directory) for directory in
                      glob(os.path.join(self.functional_dir(), 'task[0-9][0-9][0-9]_run[0-9][0-9][0-9]'))]

        for task, run_dir in task_pairs:
            if task in self._dir_tree['functional']:
                self._dir_tree['functional'][task].append(run_dir)
            else:
                self._dir_tree['functional'][task] = [run_dir]

    def dir_tree(self, sub_dir=None):
        if sub_dir:
            return self._dir_tree[sub_dir]
        else:
            return self._dir_tree

    def __repr__(self):
        return "<Subject Dir: sub{:0>3d}>".format(self._subject_code)

	def remove_volumes_from_bold(self, model_num, bold_suffix, task_to_trim, num_of_volumes):
	    from nipype.interfaces.nipy.preprocess import Trim
        for run_dir in self.dir_tree('functional')[task_to_trim]:
            file = glob(run_dir + "/bold{}.nii.gz".format(bold_suffix))[0]
            out_file = file.replace('.nii.gz', '_{}trim.nii.gz'.format(num_of_volumes))
            if os.path.isfile(out_file):
                print "removing volumes already ran!"
                return
            trim = Trim()
            trim.inputs.in_file = file
            trim.inputs.begin_index = num_of_volumes  # remove X first volumes
            trim.inputs.out_file = out_file
            res = trim.run()

    def remove_volumes_from_model(self, model_num, bold_suffix, task_to_trim, num_of_volumes):
        for onset_file in glob(os.path.join(self.model_dir(model_num), 'onsets', '{}*'.format(task_to_trim), '*.txt')):
            data = pd.read_csv(onset_file, '\t', names=['onset', 'duration', 'w'])
            data['onset'] = data['onset'] - (num_of_volumes) * 2  # TODO: fix! assuming TR of 2
            data.to_csv(onset_file, sep='\t', index=False, header=False)

    def __create_subject_dirtree__(self):
        num_models = 3 + 1  # we need the +1 because we don't start from index 0

        onset_dirs = ["model{:0>3d}/onsets".format(directory) for directory in range(1, num_models)]

        dir_tree = {
            self._path: self._subdirs.values(),
            self.functional_dir(): self._task_order,
            self.model_dir(): [os.path.join(onset_dir, task) for onset_dir in onset_dirs for task in self._task_order],
            self.masks_dir(): ['anatomy'] + self._task_order
        }

        for path, directories in dir_tree.iteritems():
            for directory in directories:
                if not os.path.isdir(os.path.join(path, directory)):
                    os.makedirs(os.path.join(path, directory))

        self.__dcm_to_nii__()
        self.__create_evs__()

    def __create_evs__(self, trim_volumes = 0, TR = 2, dummy=False):
        for file_regex, task_sequence in self._task_mapping.iteritems():
            run_files = sorted(glob(os.path.join(self.behavioural_path(), '{}*.csv'.format(file_regex)))
                               )
            ### From here on out this function is specific to my experiment. parsing logs from psychopy
            print "{}->{}".format(file_regex, run_files)
            for run_number, run_file in enumerate(run_files):
                df = pd.read_csv(run_file)
                conds = {}

                if 'keyPressed' in df.columns:
                    conds['catch'] = df[df.keyPressed.notnull()]
                    df = df[-(df.stim1.str.contains('atch') | df.stim2.str.contains('atch') | df.stim3.str.contains(
                        'atch') | df.stim4.str.contains('atch')) & df.keyPressed.isnull()]

                df['cond'] = df.stim1.str.split('\\').str[1]
                conds.update(dict((cond, df[df.cond.eq(cond)]) for cond in df.cond.unique()))
                for condition_number, key in enumerate(sorted(conds.keys())):
                    value = conds[key]
                    col_idx = df.columns.get_loc("start stim")
                    value = value.ix[:, col_idx:col_idx + 1]
                    value['start stim'] = value['start stim'].astype(int)
                    value['start stim']= value['start stim'] - (trim_volumes) * TR
                    if 'MVPA' in run_file:
                        value['length'] = 4
                    # TODO: Replace this shit!
                    else:
                        value['length'] = 12
                    # TODO: Replace this shit!

                    value['value'] = 1
                    conds[key] = value

                    # run_name = "{}_run{:0>3d}".format(task_sequence, run_number+1)
                    condition_name = "cond{:0>3d}.txt".format(condition_number + 1)
                    print ">>> condition mapping: {}->{}".format(key, condition_number)
                    # condition mapping for one beta per condition
                    value.to_csv(os.path.join(self.model_dir(1), 'onsets', task_sequence, condition_name), sep='\t',
                                 index=False, header=False)
                    # TODO: add condition mapping for one beta per trial

    def __dcm_convert__(self, source_directory, target_directory, target_filename, rename_prefix, erase=False):
        cmd = "dcm2nii -o {} {} > /dev/null".format(target_directory, source_directory)
        os.system(cmd)

        nii_file = glob("{}/{}*".format(target_directory, rename_prefix))[0]
        os.rename(nii_file, os.path.join(target_directory, '{}.nii.gz'.format(target_filename)))

        if erase:
            for file_name in glob("{}/*".format(target_directory)):
                if target_filename not in file_name:
                    os.remove(file_name)

    def create_bias_fieldmap(self):
        if os.path.exists(os.path.join(self.fieldmap_dir(),'phase.nii.gz')):
            return False
        fieldmap_dirs = sorted(glob("{}/*gre_field_mapping*".format(self._raw_path)))
        magnitude_dir = fieldmap_dirs[0]
        phase_dir = fieldmap_dirs[1]
        mag1_dir = os.path.join(magnitude_dir, 'mag1')
        mag2_dir = os.path.join(magnitude_dir, 'mag2')
        if(not os.path.exists(mag1_dir)):
            second_fieldmap_files = sorted(glob("{}/*_1*.dcm".format(magnitude_dir)))
            all_fieldmap_files = glob("{}/*.dcm".format(magnitude_dir))
            first_fieldmap_files = sorted(list(set(all_fieldmap_files) - set(second_fieldmap_files)))

            os.mkdir(mag1_dir)
            os.mkdir(mag2_dir)
            for file in first_fieldmap_files:
                shutil.move(file,mag1_dir)
            for file in second_fieldmap_files:
                shutil.move(file,mag2_dir)
        if not os.path.exists(self.fieldmap_dir()):
            os.mkdir(self.fieldmap_dir())
        self.__dcm_convert__(mag1_dir, self.fieldmap_dir(), 'mag1', '20')
        self.__dcm_convert__(mag2_dir, self.fieldmap_dir(), 'mag2', '20')
        self.__dcm_convert__(phase_dir, self.fieldmap_dir(), 'phase', '20')
        return True


    def __dcm_to_nii__(self, dummy=False):
        raw_anatomical = glob("{}/*MPRAGE*".format(self._raw_path))[0]
        raw_functional_dirs = sorted(glob("{}/*ep2*".format(self._raw_path)))
        raw_functional_dirs = filter(lambda x: "ignore" not in x, raw_functional_dirs)
        print "Converting DCM to NII"

        if dummy:
            print "self.__dcm_convert__(raw_anatomical, self.anatomical_dir(), 'highres001', 'co', True)"
        else:
            self.__dcm_convert__(raw_anatomical, self.anatomical_dir(), 'highres001', 'co', True)

        if dummy:
            print "creating functional niis"
        else:
            for raw_functional_directory, run_name in zip(raw_functional_dirs, self._task_order):
                self.__dcm_convert__(raw_functional_directory,
                                     os.path.join(self.functional_dir(), run_name),
                                     'bold',
                                     '',
                                     False)

    def path(self):
        return self._path

    def raw_path(self):
        return self._raw_path

    def behavioural_path(self):
        return self._behavioural_path

    def masks_dir(self):
        return os.path.join(self._path, self._subdirs['masks'])

    def anatomical_dir(self):
        return os.path.join(self._path, self._subdirs['anatomical'])

    def fieldmap_dir(self):
        return os.path.join(self._path, self._subdirs['fieldmap'])

    def anatomical_brain_nii(self):
        return os.path.join(self.anatomical_dir(), 'highres001_brain.nii.gz')

    def anatomical_nii(self):
        return os.path.join(self.anatomical_dir(), 'highres001.nii.gz')

    def model_dir(self, modelnum=None):
        if modelnum:
            return os.path.join(self._path, self._subdirs['model'], "model{:0>3d}".format(modelnum))
        return os.path.join(self._path, self._subdirs['model'])

    def functional_dir(self):
        return os.path.join(self._path, self._subdirs['functional'])

    def subcode(self):
        return self._subject_code


def test():
    pass


if __name__ == "__main__":
    test()
