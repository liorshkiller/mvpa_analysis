
from os.path import join as _opj
from mvpa2.mappers.detrend import poly_detrend
from mvpa2.mappers.zscore import zscore
from mvpa2.datasets.sources.openfmri import OpenFMRIDataset
from mvpa2.datasets.sources.openfmri import _model2id, _cond2id, _taskrun

from mvpa2.datasets.eventrelated import events2sample_attr
from mvpa2.datasets.eventrelated import fit_event_hrf_model
from mvpa2.base.dataset import vstack
import numpy as np
from collections import Counter

def get_bold_run_model(dhandle, model, subj, run):
    """Return the stimulation design for a particular subject/task/run.

    Parameters
    ----------
    model : int
      Model identifier.
    subj : int
      Subject identifier.
    run : int
      Run ID.

    Returns
    -------
    list
      One item per event in the run. All items are dictionaries with the
      following keys: 'condition', 'onset', 'duration', 'intensity',
      'run', 'task', 'trial_idx', 'ctrial_idx', where the first is a
      literal label, the last four are integer IDs, and the rest are
      typically floating point values. 'onset_idx' is the index of the
      event specification sorted by time across the entire run (typically
      corresponding to a trial index), 'conset_idx' is analog but contains
      the onset index per condition, i.e. the nth trial of the respective
      condition in a run.
    """

    conditions = dhandle.get_model_conditions(model)
    events = []
    ev_fields = ('onset', 'duration', 'intensity','id')

    # get onset info for specific subject/task/run combo
    for cond in conditions:
        task_id = cond['task']
        def _load_hlpr(fname):
            return np.recfromtxt(fname, names=ev_fields)


        evdata = np.atleast_1d(dhandle._load_subj_data(
            subj,
            ['model', _model2id(model), 'onsets',
             _taskrun(1, run), '%s.txt' % _cond2id(cond['id'])],
            _load_hlpr)
                )
        for i, ev in enumerate(evdata):
            evdict = dict(zip(ev_fields,
                              [ev[field] for field in ev_fields]))
            evdict['task'] = task_id
            evdict['condition'] = cond['name']
            evdict['run'] = run
            evdict['conset_idx'] = i
            events.append(evdict)
    events = sorted(events, key=lambda x: x['onset'])
    for i, ev in enumerate(events):
        ev['onset_idx'] = i
    return events

def create_betas_per_trial_with_pymvpa(study_path, subj, conf, mask_name, flavor, TR):
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
											  chunks = run_id -1,
											  mask = mask_fname)
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

def create_betas_per_trial_with_pymvpa_roni(study_path, subj, conf, mask_name, flavor, TR):
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

		#all_events = dhandle.get_bold_run_model(model, subj, run_id)
		all_events = get_bold_run_model(dhandle,2, subj, run_id)
		run_events = []
		i=0
		for event in all_events:
			if event['task']==task:
				event['condition'] = '{}-{}'.format(event['condition'],event['id'])
				run_events.append(event)
				i+=1

		# load BOLD data for this run (with masking); add 0-based chunk ID
		run_ds = dhandle.get_bold_run_dataset(subj, task, run_id,
											  flavor = flavor,
											  chunks = run_id -1,
											  mask = mask_fname)
		# convert event info into a sample attribute and assign as 'targets'
		run_ds.sa.time_coords = run_ds.sa.time_indices*TR
		run_ds.sa['targets'] = events2sample_attr(
					run_events, run_ds.sa.time_coords, noinfolabel='rest')
		# additional time series preprocessing can go here
		poly_detrend(run_ds, polyord=1, chunks_attr='chunks')
		zscore(run_ds, chunks_attr='chunks', param_est=('targets', ['rest']), dtype='float32')
		glm_dataset = fit_event_hrf_model(run_ds, run_events,time_attr='time_coords',condition_attr='condition')
		glm_dataset.sa['targets'] = [x[:x.find('-')] for x in glm_dataset.sa.condition]
		glm_dataset.sa['id'] = [x[x.find('-')+1:] for x in glm_dataset.sa.condition]
		glm_dataset.sa.condition = glm_dataset.sa['targets']
		glm_dataset.sa['chunks'] = [run_id -1] * len(glm_dataset.samples)

		# If a trial was dropped (the subject pressed on a button) than the counter trial from the
		# other condition should also be dropped
		for pair in conf.conditions_to_compare:
			cond_bool = np.array([c in pair for c in glm_dataset.sa['condition']])
			sub_dataset = glm_dataset[cond_bool]
			c = Counter(sub_dataset.sa.id)
			for value in  c:
				if c[value] < 2:
					id_bool = np.array([value in cond_id for cond_id in glm_dataset.sa['id']])
					glm_dataset = glm_dataset[np.bitwise_not(np.logical_and(id_bool,cond_bool))]

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