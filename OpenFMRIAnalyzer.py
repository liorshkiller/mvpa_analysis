#!/usr/bin/python

import os
import nibabel
import shutil
import subprocess
import signal
import pre_proc_fromexample as pp
from OpenFMRIData import OpenFMRIData

from nipype.interfaces import fsl


class OpenFMRIAnalyzer(object):
    def __init__(self, fmri_data, conf, subjects=[]):
        self._fmri_data = fmri_data
        self.subjects_list = []

	self.conf = conf

        self.__load_subjects__(subjects)

    def __load_subjects__(self, subjects):
        if len(subjects) == 0:
            subjects = self._fmri_data.all_subjects_list()
        for subject in subjects:
            self.subjects_list.append(self._fmri_data.subject_dir(subcode=subject))

    def generate_functional_gm_masks(self, subject):
        print ">>> Creating functional gray matter masks"
        mask_name = 'grey.nii.gz'
        gm_mask = os.path.join(subject.masks_dir(), 'anatomy', mask_name)

        for task, directories in subject.dir_tree('functional').iteritems():
            if task in self.conf.mvpa_tasks:
                output_dir = os.path.join(subject.masks_dir(), task)
                output_file = os.path.join(output_dir, mask_name)
                if os.path.exists(output_file):
                    continue
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                task_dir = os.path.join(subject.functional_dir(), task)
                reg_dir = os.path.join(task_dir, "reg")
                gm2func_mask = fsl.preprocess.ApplyXfm()
                gm2func_mask.inputs.in_matrix_file = os.path.join(reg_dir, 'highres2example_func.mat')
                gm2func_mask.inputs.reference = os.path.join(task_dir, 'refvol.nii.gz')
                gm2func_mask.inputs.in_file = gm_mask
                gm2func_mask.inputs.out_file = output_file
                gm2func_mask.run()
                #gen_mask = fsl.utils.ImageMaths(in_file=output_file ,op_string = '-thr {} -bin'.format(0.5), out_file=output_file)
                #print gen_mask.cmdline
                #gen_mask.run()

    def __register_functional_file__(self, subject, func_file, brain_image, directory):
        reg_dir = os.path.join(directory, 'reg')
        if os.path.isfile(os.path.join(reg_dir, 'example_func2standard_warp.nii.gz')):
            print "registration for {} already performed".format(directory)
            return
        print "working on {}".format(directory)
        if not os.path.isdir(reg_dir):
            os.mkdir(reg_dir)
        else:
            shutil.rmtree(reg_dir)
        log_file = os.path.join(directory, 'log_reg')
        cmd = 'mainfeatreg -F 6.00 -d {} -l {} -i {} -h {} -w BBR -x 90 > /dev/null'.format(directory, log_file,
                                                                                            func_file, brain_image)
        #		cmd = 'mainfeatreg -F 6.00 -d {} -l {} -i {} -h {} -w 6 -x 90  > /dev/null'.format(directory,log_file, func_file, brain_image)
        subprocess.call(cmd, shell=True)
        anat_reg_dir = os.path.join(subject.anatomical_dir(), 'reg')
        highres2mni_mat = os.path.join(anat_reg_dir, 'highres2standard.mat')
        highres2standard_warp = os.path.join(anat_reg_dir, 'highres2standard_warp.nii.gz')
        example_func2highres_mat = os.path.join(reg_dir, 'example_func2highres.mat')
        example_func2standard_warp = os.path.join(reg_dir, 'example_func2standard_warp.nii.gz')

        standard_image = fsl.Info.standard_image('MNI152_T1_2mm_brain.nii.gz')
        convert_warp = fsl.utils.ConvertWarp(reference=standard_image,
                                             premat=example_func2highres_mat,
                                             warp1=highres2standard_warp,
                                             out_file=example_func2standard_warp)
        convert_warp.run()
        apply_warp = fsl.preprocess.ApplyWarp(ref_file=standard_image,
                                              in_file=func_file,
                                              field_file=example_func2standard_warp,
                                              out_file=os.path.join(reg_dir, 'example_func2standard.nii.gz'))
        apply_warp.run()

    def functional_registration(self, subject):
        print ">>> Functional Registration"
        brain_image = subject.anatomical_brain_nii()
        
        for task, directories in subject.dir_tree('functional').iteritems():
            if task in self.conf.mvpa_tasks:
                task_dir = os.path.join(subject.functional_dir(), task)
                ref_vol = os.path.join(task_dir, "refvol.nii.gz")
                input_file = ref_vol.replace(".nii.gz","_restore.nii.gz")
                fast = fsl.FAST(in_files=ref_vol,
                            out_basename=os.path.join(task_dir, 'refvol'),
                            bias_lowpass=10,
                            output_biascorrected=True,
                            output_biasfield=True,
                            img_type=2,
                            bias_iters=10,
                            no_pve=True,
                            iters_afterbias=1)
                try:
                    if not os.path.exists(input_file):
                        fast = fast.run()
                except:
                    pass

                self.__register_functional_file__(subject, input_file, brain_image, task_dir)

    def anatomical_registration(self, subject):
        print ">>> Anatomical registration"
        brain_image = subject.anatomical_brain_nii()
        reg_dir = os.path.join(subject.anatomical_dir(), 'reg')
        out_file = os.path.join(reg_dir, 'highres2standard.nii.gz')
        out_mat_file = os.path.join(reg_dir, 'highres2standard.mat')
        standard_image = fsl.Info.standard_image('MNI152_T1_2mm_brain.nii.gz')

        if not os.path.isfile(out_mat_file):
            print ">>> FLIRT"
            if not os.path.exists(reg_dir):
                os.mkdir(reg_dir)
            flirt = fsl.FLIRT(in_file	= brain_image,
                              reference	 = standard_image,
                              out_file 	 = os.path.join(reg_dir ,"flirt.nii.gz")  ,  # out_file,
                              out_matrix_file= out_mat_file,
                              cost		 = 'corratio',
                              dof		 = 12,
                              searchr_x	 = [-90, 90],
                              searchr_y	 = [-90, 90],
                              searchr_z	 = [-90, 90],
                              interp	 ='trilinear')
            flirt.run()

        anatomical_head = os.path.join(subject.anatomical_dir() ,'highres001.nii.gz')
        output_fielf_coeff = os.path.join(reg_dir, 'highres2standard_warp.nii.gz')
        output_jacobian = os.path.join(reg_dir, 'highres2highres_jac')
        standard_head = fsl.Info.standard_image('MNI152_T1_2mm.nii.gz')
        standard_mask = fsl.Info.standard_image('MNI152_T1_2mm_brain_mask_dil.nii.gz')
        if not os.path.isfile(output_fielf_coeff):
            print ">>> FNIRT"
            fnirt = fsl.FNIRT(warped_file	 = out_file,
                              # in_fwhm        = [10,4,2,2],

                              in_file	 = anatomical_head,
                              affine_file	 = out_mat_file,
                              fieldcoeff_file= output_fielf_coeff,
                              jacobian_file	 = output_jacobian,
                              config_file	 = 'T1_2_MNI152_2mm',
                              ref_file	 = standard_head,
                              refmask_file	 = standard_mask
                              )
            fnirt.run()
            cmd = 'fslview {} {} -t 0.5 '.format(standard_image ,out_file)
            pro = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   shell=True, preexec_fn=os.setsid)
        inv_mat_file = os.path.join(reg_dir, 'standard2highres.mat')
        if not os.path.exists(inv_mat_file):
            invt = fsl.ConvertXFM()
            invt.inputs.in_file = out_mat_file
            invt.inputs.invert_xfm = True
            invt.inputs.out_file = inv_mat_file
            invt.run()

    def fieldmap_preperation(self,subject):
        mag_file = os.path.join(subject.fieldmap_dir(), "mag2.nii.gz")
        mag_brain = mag_file.replace("nii.gz","_brain.nii.gz")
        bet = fsl.BET(in_file=mag_file,
                      out_file=mag_brain,
                      robust=True,
                      frac=0.5)
        bet.run()
        prepare = fsl.PrepareFieldmap()
        prepare.inputs.in_phase = os.path.join(subject.fieldmap_dir(),"phase.nii")
        prepare.inputs.in_magnitude = mag_brain
        prepare.inputs.output_type = "NIFTI_GZ"
        prepare.inputs.out_fieldmap = os.path.join(subject.fieldmap_dir(),"fieldmap.nii")
        prepare.cmdline
        res = prepare.run()


    def __gm_mask_from_functional__(self, input_file, mask_dir):
        gm_mask_name = os.path.join(mask_dir, 'grey.nii.gz')
        if os.path.isfile(gm_mask_name):
            return
        out_basename = os.path.join(mask_dir, 'seg')
        input_file_brain = input_file.replace('.nii.gz' ,'_brain.nii.gz')
        bet = fsl.BET(in_file=input_file,
                      out_file=input_file_brain,
                      mask=True,
                      robust=True,
                      frac=0.3)
        bet.run()
        fast = fsl.FAST(in_files		= input_file_brain,
                        out_basename		= out_basename,
                        img_type		= 2,
                        number_classes		= 3,
                        hyper			= 0.1,
                        output_biascorrected	= True,
                        output_biasfield	= True,
                        bias_iters		= 5,
                        iters_afterbias		= 2,
                        segments		= True)
        try:
            result = fast.run()
            gm_pve_file = result.outputs.partial_volume_files[0]
        except:
            gm_pve_file = '{}_pve_0.nii.gz'.format(out_basename)
        try:
            os.rename(gm_pve_file ,gm_mask_name)
            cmd = 'fslview {} {} -l Red '.format(input_file_brain,gm_mask_name)
            pro = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   shell=True, preexec_fn=os.setsid)
        except:
            pass

    def standard_mask2native_space(self, mask_file, task_dir, ref_vol,mask_dir):
        out_mask_file = os.path.join(mask_dir, os.path.basename(mask_file))
        print out_mask_file
        if os.path.exists(out_mask_file):
            return
        warp_file = os.path.join(task_dir,'reg', 'standard2example_func_warp.nii.gz')
        invwarp = fsl.utils.InvWarp(reference=ref_vol,
                                    warp=os.path.join(task_dir,'reg', 'example_func2standard_warp.nii.gz'),
                                    inverse_warp=warp_file)
        invwarp.run()
        apply_warp = fsl.preprocess.ApplyWarp(ref_file=ref_vol,
                                              in_file=mask_file,
                                              field_file=warp_file,
                                              out_file=out_mask_file)
        apply_warp.run()

    def warp_standard_mask(self, subject):
        print ">>> Functional Segmentation"
        standard_mask = self.conf.standard_mask

        for task, directories in subject.dir_tree('functional').iteritems():
            task_dir = os.path.join(subject.functional_dir(), task)
            ref_vol  = os.path.join(task_dir, "refvol.nii.gz")
            mask_dir = os.path.join(subject.masks_dir(), task)
            if task in self.conf.mvpa_tasks:
                task_dir = os.path.join(subject.functional_dir(), task)
                ref_vol  = os.path.join(task_dir, "refvol.nii.gz")
                mask_dir = os.path.join(subject.masks_dir(), task)
                if not os.path.isdir(mask_dir):
                    os.mkdir(mask_dir)
                self.standard_mask2native_space(standard_mask, task_dir, ref_vol, mask_dir)

    def functional_segmentation(self, subject):
        print ">>> Functional Segmentation"

        for task, directories in subject.dir_tree('functional').iteritems():
            if task in self.conf.mvpa_tasks:
                task_dir = os.path.join(subject.functional_dir(), task)
                ref_vol  = os.path.join(task_dir, "refvol_restore.nii.gz")
                mask_dir = os.path.join(subject.masks_dir(), task)
                if not os.path.isdir(mask_dir):
                    os.mkdir(mask_dir)
                self.__gm_mask_from_functional__(ref_vol, mask_dir)

    def segmentation(self, subject):
        print ">>> Segmentation"
        mask_dir = os.path.join(subject.masks_dir() ,'anatomy')
        gm_mask_name = os.path.join(mask_dir ,'grey.nii.gz')
        if os.path.exists(gm_mask_name):
            return
        if not os.path.exists(mask_dir):
            os.makedirs(mask_dir)
        brain_image = os.path.join(subject.anatomical_dir() ,"highres001_brain.nii.gz")
        fast = fsl.FAST(in_files=brain_image, out_basename=os.path.join(subject.masks_dir() ,'anatomy' ,'seg'), img_type=1, number_classes=3, hyper=0.4 ,segments=True)

        try:
            result = fast.run()
            gm_pve_file = result.outputs.partial_volume_files[1]
        except:
            gm_pve_file = os.path.join(subject.masks_dir() ,'anatomy' ,'seg_pve_1.nii.gz')
            gm_seg_file = os.path.join(subject.masks_dir() ,'anatomy' ,'seg_seg_1.nii.gz')

        # TODO: bug in fast result output parsing!!!
        if False:
            cmd = 'fslview {} {} -l Red -t 0.1 -b 0,0.1'.format(brain_image ,gm_pve_file)
            pro = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   shell=True, preexec_fn=os.setsid)
            thr = float(raw_input("Insert a GM thershold for the mask: default is 0\n")) or 0.0
            os.killpg(pro.pid ,signal.SIGTERM)
            gen_mask = fsl.utils.ImageMaths(in_file=gm_pve_file ,op_string = '-thr {} -bin'.format(thr), out_file=gm_mask_name)
            gen_mask.run()
        else:
            os.rename(gm_seg_file ,gm_mask_name)


    def estimate_bias_field(self, subject):
        print ">>> Bias field estimation"

        anat_filename = os.path.join(subject.anatomical_dir(), 'highres001.nii.gz')
        restore_file = os.path.join(subject.anatomical_dir(), 'highres001_restore.nii.gz')

        if os.path.isfile(restore_file):
            return anat_filename

        try:
            fast = fsl.FAST(in_files=anat_filename,
                            out_basename=os.path.join(subject.anatomical_dir(), 'highres001'),
                            bias_lowpass=10,
                            output_biascorrected=True,
                            output_biasfield=True,
                            img_type=1,
                            bias_iters=5,
                            no_pve=True,
                            iters_afterbias=1)

            fast = fast.run()

            # TODO need to remove temporary files created by this process

            os.rename(anat_filename ,anat_filename.replace('.nii.gz' ,'_pre_restore.nii.gz'))
            shutil.copy(fast.outputs.image_restored ,anat_filename)
            return anat_filename
        except:
            os.rename(anat_filename ,anat_filename.replace('.nii.gz' ,'_pre_restore.nii.gz'))
            shutil.copy(restore_file, anat_filename)
            return anat_filename

    def extract_brain(self, subject):
        # Check whether brain has already been extracted
        brain_image = os.path.join(subject.anatomical_dir(), 'highres001_brain.nii.gz')

        if os.path.isfile(brain_image):
            return

        # Estimate bias field - removed for now
        # input_image = self.estimate_bias_field(subject)
        input_image = os.path.join(subject.anatomical_dir(), 'highres001.nii.gz')

        print ">>> Brain Extraction"

        f = 0.5
        g = -0.1

        bet = fsl.BET(in_file=input_image,
                      out_file=brain_image,
                      mask=True,
                      robust=True,
                      frac=f,
                      vertical_gradient=g)

        result = bet.run()

        is_ok = 'n'
        while 'n' in is_ok:
            cmd = 'fslview {} {} -l Green'.format(input_image, result.outputs.out_file)
            pro = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, preexec_fn=os.setsid)

            is_ok = raw_input("Is this ok? [y]/n\n") or 'y'
            if 'n' in is_ok:
                bet.inputs.frac = float \
                    (raw_input("Set fraction: default is previous ({})\n".format(bet.inputs.frac)) or bet.inputs.frac)
                bet.inputs.vertical_gradient = float(raw_input("Set gradient: default is previous ({})\n".format
                                                               (bet.inputs.vertical_gradient)) or bet.inputs.vertical_gradient)
                result = bet.run()
            os.killpg(pro.pid ,signal.SIGTERM)
        os.rename(os.path.join(subject.anatomical_dir() ,'highres001_brain_mask.nii.gz')
                  ,os.path.join(subject.masks_dir() ,'anatomy' ,'brain.nii.gz'))


    def __motion_correct_file__(self, input_file, output_file ,subject ,directory ,refvol=None, use_example_pp=True):
        # Check whether motion correction has already been completed
        if os.path.isfile(output_file):
            return

        print "{}".format(input_file)

        if use_example_pp:
            pp.preproc.inputs.inputspec.func = input_file
            pp.preproc.inputs.inputspec.struct = os.path.join(subject.anatomical_dir() ,'highres001.nii.gz')
            pp.preproc.base_dir = directory
            if refvol:
                pp.preproc.get_node('realign').inputs.ref_file = refvol
            pp.preproc.run()
            # TODO: copy motion correction photos as well
            shutil.copy \
                (os.path.join(directory ,'preproc' ,'maskfunc2' ,'mapflow' ,'_maskfunc20' ,'bold_dtype_mcf_mask.nii.gz')
                 ,output_file)
            cmd = "eog {}".format \
                (os.path.join(directory ,'preproc' ,'realign' ,'mapflow' ,'_realign0' ,'bold_dtype_mcf.nii.gz_rot.png'))
            subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
            cmd = "eog {}".format(os.path.join(directory ,'preproc' ,'realign' ,'mapflow' ,'_realign0'
                                               ,'bold_dtype_mcf.nii.gz_trans.png'))
            subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, preexec_fn=os.setsid)

        else:
            mcflt = fsl.MCFLIRT(in_file=input_file, out_file=output_file, save_plots=True)
            result = mcflt.run()

            pmp = fsl.PlotMotionParams(in_file = result.outputs.par_file ,in_source='fsl')

            pmp.inputs.plot_type = 'rotations'
            pmp.run()
            pmp.inputs.plot_type = 'translations'
            pmp.run()

    def motion_correction(self, subject):
        print ">>> Motion correction"
        
        for task, directories in subject.dir_tree('functional').iteritems():
            if task in self.conf.mvpa_tasks:
                bold_files = [os.path.join(directory, 'bold.nii.gz') for directory in directories]
                bold_files = sorted(bold_files)
                mcf_files = [os.path.join(directory, 'bold_mcf.nii.gz') for directory in directories]
                if all(map(lambda x: os.path.isfile(x), mcf_files)):
                    print ">>> Motion Correction has already been performed"
                    continue

                mid_run_num = len(bold_files ) /2
                task_dir = os.path.join(subject.functional_dir() ,task)
                if not os.path.isdir(task_dir):
                    os.mkdir(task_dir)
                refvol = os.path.join(task_dir ,'refvol.nii.gz')
                extract_first = fsl.ExtractROI(in_file	= bold_files[mid_run_num],
                                               roi_file	= refvol,
                                               t_min	= 0,
                                               t_size	= 1)
                result = extract_first.run()
                for directory in sorted(directories):
                    input_file = os.path.join(directory, 'bold.nii.gz')
                    output_file = input_file.replace('.nii.gz', '_mcf.nii.gz')
                    self.__motion_correct_file__(input_file, output_file,
                                                 subject   , directory  , refvol)

