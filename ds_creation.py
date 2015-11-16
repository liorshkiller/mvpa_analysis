
from mvpa2.suite import *
from os.path import join as _opj

def create_betas_per_trial_with_pymvpa(study_path, subj, conf, mask_name, flavor,TR):
	dhandle = OpenFMRIDataset(study_path)
	model = 1
	task = 1
	# Do this for other tasks as well. not only the first
	mask_fname = _opj(study_path,"sub{:0>3d}".format(subj),"masks",conf.mvpa_tasks[0],"{}.nii.gz".format(mask_name))
	print mask_fname
	run_datasets = []
	for run_id in dhandle.get_task_bold_run_ids(task)[subj]:
		if type(run_id) == str:
			continue
		all_events = dhandle.get_bold_run_model(model, subj, run_id)
		run_events = []
		i=0
		for event in all_events:
			if event['task']==task:
				event['condition'] = '{}-{}'.format(event['condition'],i)
				run_events.append(event)
				i+=1

		# load BOLD data for this run (with masking); add 0-based chunk ID
		run_ds = dhandle.get_bold_run_dataset(subj, task, run_id,
											  flavor = flavor,
											  chunks=run_id -1,
											  mask=mask_fname)
		# convert event info into a sample attribute and assign as 'targets'
		run_ds.sa.time_coords = run_ds.sa.time_indices*TR
		print run_id

		run_ds.sa['targets'] = events2sample_attr(
					run_events, run_ds.sa.time_coords, noinfolabel='rest')
		# additional time series preprocessing can go here
		poly_detrend(run_ds, polyord=1, chunks_attr='chunks')
		zscore(run_ds, chunks_attr='chunks', param_est=('targets', ['rest']), dtype='float32')
		glm_dataset = fit_event_hrf_model(run_ds, run_events,time_attr='time_coords',condition_attr='condition')
		glm_dataset.sa['targets'] = [x[:x.find('-')] for x in glm_dataset.sa.condition]
		glm_dataset.sa.condition = glm_dataset.sa['targets']
		glm_dataset.sa['chunks'] = [run_id -1] * len(glm_dataset.samples)
		run_datasets.append(glm_dataset)
	return vstack(run_datasets, 0)

def detrend(ds):
	poly_detrend(ds, polyord=1, chunks_attr='chunks')
	zscore(ds, chunks_attr='chunks', dtype='float32')
	return ds

def create_betas_per_run_with_pymvpa(study_path, subj, conf, mask_name, flavor):
	of = OpenFMRIDataset(study_path)
	mask_fname = _opj(study_path,"sub{:0>3d}".format(subj),"masks",conf.mvpa_tasks[0],"{}.nii.gz".format(mask_name))
	ds = of.get_model_bold_dataset(
    model_id=1, subj_id=subj,
    flavor=flavor,
    mask= mask_fname,
    # preproc_img=smooth,
    preproc_ds=detrend,
    modelfx=fit_event_hrf_model,
    time_attr='time_coords',
    condition_attr='condition')
	return ds