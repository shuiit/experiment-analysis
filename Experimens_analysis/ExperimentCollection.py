from Experiment import Experiment

# -*- coding: utf-8 -*-
"""
Created on Sun Jun  4 14:05:35 2023

@author: Roni
"""

class ExperimentCollection():
    """
    A class for loading multiple experiments.
    """
    
    def __init__(self):
        """
        initilize ExperimentCollection class

        Returns
        -------
        None.

        """
        self.experiments = {}

    def add_experiment(self,loadir,exp_date,exp_name,pert_time = [None, None]):
        """
        add an experiment to the collection of experiments

        Parameters
        ----------
        loadir (str): loading diractory
        exp_date (str): date of experiment
        exp_name (str): name of experiment
        pert_time (list of ints, optional): The default is [None, None].

        Returns
        -------
        None.

        """
        experiment = Experiment(loadir,exp_date,exp_name,pert_time)   # load experimens
        self.experiments[exp_name] = experiment
        
    def get_experiment(self, exp_name):
        """
        get experiment from the collection

        Parameters
        ----------
        exp_name (str): name of experiment

        Returns
        -------
        experiment from colecction (Experiment class)

        """
        return self.experiments.get(exp_name)
        
    def delete_experiment(self, exp_name):
        """
        Delete experiment from the collection

        Parameters
        ----------
        exp_name (str): name of the experiment

        Returns
        -------
        None.

        """
        if exp_name in self.experiments:
            del self.experiments[exp_name]
    
    def iterate_experiments(self):
        """
        Creates a generator to iterate on all experiments

        Yields
        ------
        exp_name (str): name of experiment
        experiment (Experiment object): object that has all the data of the experiment

        """
        for exp_name, experiment in self.experiments.items():
            yield exp_name, experiment

        