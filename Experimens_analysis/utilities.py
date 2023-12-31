# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:58:34 2023

@author: Roni
"""




def load_experiments(collection,exp_date,exp_name,pert_time,loadir):
    """
    

    Parameters
    ----------
    collection (ExperimentCollection object): object of all experiments
    exp_date (str): date of the experiment
    exp_name (str): name of the experiment
    pert_time (list of ints): pertubation start and end
    loadir (str): loading directory.

    Returns
    -------
    None.

    """
    # Add experiments to the collection
    [collection.add_experiment(loadir, exp_date, exp_name, pert_time) for exp_date,exp_name,pert_time in zip(exp_date,exp_name,pert_time)]
    
    for experiment_id, experiment in collection.iterate_experiments():
        experiment = collection.get_experiment(experiment_id)
        experiment.load_data()

    

def clean_data(collection,input_del_short_mov ,manual_del ,input_filter_manouvering):
    """
    Clean data for each experiment in the collection

    Parameters
    ----------
    collection (ExperimentCollection object): object of all experiments
    input_del_short_mov (list of ints 1 X 2): input_del_short_mov[0] = length of movie
                                              input_del_short_mov[1] = min(time) (None - do nothing)
    manual_del (dict (str:list of str)): date of experiment: [name of moves]
    input_filter_manouvering (list 1 X 4) : input_filter_manouvering[0] (str) = property of manouver
                                            input_filter_manouvering[1] (list of int) = time vector to sample
                                            input_filter_manouvering[2] (float) = value of manouver to filter by 
                                            input_filter_manouvering[3] (bool) = plot filtered movies

    Returns
    -------
    None.

    """
 
    for experiment_id, experiment in collection.iterate_experiments():
        experiment = collection.get_experiment(experiment_id)
        # Delete movies based on maneuvering, movie length, and manually
        experiment.clean_data(input_del_short_mov ,manual_del ,input_filter_manouvering)
        
        
def manipulate_data(collection,features_to_keep,features_to_keep_mean):
    """
    Manipulate data for each experiment:
        -project location, velocity and acceleration
        - calculate the mean features (by stroke)
        - keep only 'features_to_keep' and 'features_to_keep_mean'

    Parameters
    ----------
    collection (ExperimentCollection object): object of all experiments
    features_to_keep (list of str): properties of features to keep.
    features_to_keep_mean (list of str): properties of features to keep for mean data.

    Returns
    -------
    None.

    """
    # Manipulate data for each experiment in the collection
    for experiment_id, experiment in collection.iterate_experiments():
        experiment = collection.get_experiment(experiment_id)
        experiment.property_projection()
        experiment.meanfeature()
        experiment.data = experiment.data[features_to_keep]
        experiment.mean_data = experiment.mean_data[features_to_keep_mean]
        
        
def calculate_interest_points(collection):
    """
    Calculate interest points for each experiment in the collection
        - calculates the idx of min\max points
        - calculate the difference between value in t = 0 and othe point

    Parameters
    ----------
    collection (ExperimentCollection object): object of all experiments

    Returns
    -------
    None.

    """
    for experiment_id, experiment in collection.iterate_experiments():
        experiment = collection.get_experiment(experiment_id)
        experiment.calculate_interest_idx(time_min=50, time_max=150)
        experiment.calculate_delta('projected_x_dot', 'minimal_velocity', 'delta_v')
        experiment.calculate_delta('pitch smth', 'maximal_pitch', 'delta_pitch')
 
        