import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
import scipy.signal as signal
import itertools
import glob


# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 13:22:00 2023

@author: Roni
"""


class Experiment():
    """
    A class for loading and preprocessing data.
    """
    
    def __init__(self,load_path,exp_date,exp_name,pert_time = None):
        """
        Initializes the Experiment class.

        Parameters
        ----------
        load_path (str): Path to load the data from.
        exp_date (list of str, optional): experiment to load. If None, load all experiments.
        exp_name (str, optional): List of experiments to load. If None, load all experiments.
        pert_time (list, otional): time of pertubation ([t_0 t_end])
        mov_order (list, optional): list of movies names in the same order as in the data table
        

        Returns
        -------
        None.

        """

        self.exp_date = exp_date 
        self.exp_name = exp_name
        self.pert_time = pert_time
        self.mov_order = None
        self.mean_data = None
        
        self.load_path = load_path        
        self.data = None
        self.body_vectors = {}
        self.meanstroke = {}
        
        
        self.input_del_short_mov = None
        self.manual_del = None    
        self.input_filter_manouvering = None  
        self.interest_points = pd.DataFrame()
        
        
        


# =============================================================================
#   Load and clean
# =============================================================================
        
    def load_data(self):
        """
        Loads data from load_path into object
        
        This method loads data from files located in the 'load_path/data_analysis' directory. It iterates through the files in the directory and performs the following actions:: 
        1. If 'exp_date' is specified, only files whose experiment names are in 'exp_date' will be processed.
        2. If a file contains 'body_vec' in its name, its contents are loaded into the 'body_vectors' attribute.
        2. If a file contains 'DF_all' in its name, its contents (angles,[X Y Z]...) its contents are loaded into the 'data' attribute. 

        Returns
        -------
        None.
        
        Note:
        - this function modifies the data and body_vectors in place

        """
       
        name_of_path = self.load_path + 'data_analysis'
        file_pattern = f"{name_of_path}/*{self.exp_date}*.pkl" 
        matching_files = glob.glob(file_pattern)

        for nm in matching_files:
            with open(nm, 'rb') as file:
                if 'body_vec' in nm:
                    self.body_vectors = (pickle.load(file)) 
                if 'DF_all' in nm:
                    self.data = pickle.load(file)


    def clean_data(self,input_del_short_mov = None,manual_del = None,input_filter_manouvering = None):
        """
        This method calls movie filtering functions
        - del_short_mov: deletes movies shorter than threashold
        - del_mov_manually: deletes movies according to user input
        - input_filter_manouvering: deletes movies by manouvers threshold

        Parameters
        ----------
        
        input_del_short_mov (list, optional): Input parameters for deleting short movies. Defaults to [None, None].
            - input_del_short_mov[0] (int): Minimal length of a movie to keep. If None, keep all movies.
            - input_del_short_mov[1] (int): Minimum time of a movie to keep.


        manual_del (dict, optional): Movies to manually delete. (date of experiment: name of movies in a list example: {'2022_02_03':['mov6','mov7','mov65'],'2022_03_10':['mov37']})

        
        input_project_on_zero_bodaxs (list of ints, optional): Input parameters for properties projection on XY plane. Defaults to [211,3,400].
            input_project_on_zero_bodaxs[0] : Window for savgol_filter to plot projected velocities. 
            input_project_on_zero_bodaxs[1] : polynom degree for savgol_filter to plot projected velocities. 
            input_project_on_zero_bodaxs[2] : Frame of referance axis.
            input_project_on_zero_bodaxs[3] : plot flag - plot the filtered movies
        
        

        Returns
        -------
        None.

        """
       
        self.input_del_short_mov = input_del_short_mov
        
        self.manual_del = manual_del
        
        self.input_filter_manouvering = input_filter_manouvering
        
        
        
        if self.input_del_short_mov[0] is not None:
            self.del_short_mov()
        if self.manual_del is not None:
            self.del_mov_manually()
        if self.input_filter_manouvering['manouver_prop'] is not None:
            self.filter_manouvering()
        self.mov_table_order()
        
    def mov_table_order(self):
        
        """
        Get a vecotor of movies according to the order in the data table

        Returns
        -------
        None.
        
        Note:
        - This function modifies 'mov_order' attribute in-place.

        """
        self.mov_order,mov_idx = np.unique(self.data['mov'], return_index=True)
        self.mov_order = self.mov_order[np.argsort(mov_idx)]
    
    def del_short_mov(self):
        """
        Deletes rows (movies) from the data based on specified conditions.

        This method deletes rows from the 'data' attribute of the current object based on the following conditions:
        1. If 'min_time_del' is not None, rows with 'time [ms]' values greater than 'min_time_del' are removed.
        2. Rows with 'mov' values that appear less than 'short_mov_len' times are deleted.

        
        Returns
        -------
        None.
        
        Note:
        - This function modifies the 'data' attribute in-place.

        """
        min_time_del = self.input_del_short_mov[1]
        short_mov_len = self.input_del_short_mov[0]
        
            
        # delete movies that have minimum time grater than min_time_del
        # Delete movies that have a minimum time greater than min_time_del
        if min_time_del is not None:
            self.data = self.data.groupby('mov').filter(lambda x: x['time [ms]'].min() <= min_time_del) 
    
        # Filter out rows with 'mov' values that appear less than short_mov_len times
        self.data = self.data.groupby('mov').filter(lambda x: len(x) >= short_mov_len).reset_index(drop=True)
        self.data.drop('index', axis=1, errors='ignore', inplace=True)
        
 
 
    def del_mov_manually(self):
        
        """
        Manually delete rows from data based on user input 'manual_del'

        Returns
        -------
        None.
        
        Note:
        - This function modifies the 'data' attribute in-place.

        """
        if self.exp_date in self.manual_del:
            manual_del_values = self.manual_del[self.exp_date] 
            self.data = self.data.query('mov not in @manual_del_values')
        else:
            pass
        

        
    def find_location(self,group):
        """
        Find rows based on a group of data speciefic conditions 
        
        Parameters
        ----------
        group : (pandas.DataFrame): A DataFrame containing the data group.

        Returns
        -------
        bool
            True if the absolute sum of the specified property values in the group,
            averaged over time, exceeds the threshold. False otherwise.
            
        Notes:
        ------
        - This method calculates the sum of the specified property values in the
          group for each time step and checks if the average exceeds the threshold.
        - The property and threshold values are set as attributes in the class instance.


        """
        manouver_prop = self.input_filter_manouvering['manouver_prop']
        time_of_maneuvers = self.input_filter_manouvering['time_of_maneuvers']
        prop_threshold = self.input_filter_manouvering['prop_threshold']
        
        return abs(sum(group.query('`time [ms]` == @time_of_maneuvers')[manouver_prop]))/len(time_of_maneuvers) > prop_threshold
        
        
       
       
    
    def filter_manouvering(self,pltflg = False,col_plot = 'pitch smth'):
        """
        Filters rows based on specific conditions
        
        This method calls 'find_location' for each movie in 'data' and keeps rows that are assigned with False value 

        Parameters
        ----------
        pltflg (bool, optional): For debugging, True - plot the deleted rows. The default is False.
        col_plot (str, optional): column to plot. The default is 'pitch smth'.

        Returns
        -------
        None.
        
        Notes:
        ------
        - This function modifies the 'data' attribute in-place.
        

        """           
        rows_to_filter = self.data.groupby('mov').apply(self.find_location)
        
        
        # plot the deleted movies - use for debugging
        if self.input_filter_manouvering['plot_delted']: 
            plt.figure()
            for mvnm in rows_to_filter[rows_to_filter].index: 
                mv = self.data.groupby('mov').get_group(mvnm)
                plt.plot(mv['time [ms]'],mv[col_plot])
                    
        self.data = self.data.query('mov not in @rows_to_filter.keys()[@rows_to_filter == True]').reset_index(drop=True)

# =============================================================================
#   Stroke avereged
# =============================================================================
    def max_min_wing_angles(self,mean_data,group_mov,wing):
        """
        Calculates the max and min values of the wing angke per stroke
       
        Parameters
        ----------
        mean_data (pandas data frame): The data avereged by strokes.
        group_mov (pandas group): The data grouped by movies.
        wing (str): 'lw'/'rw' - left or right wing to update.
       
        Returns
        -------
        mean_data pandas data frame
            the data frame avereged by strokes with new columns for max/min wing angles.
       
        """
        mean_data[wing]['max_phi_' + wing] = group_mov.max().reset_index()['maxmin_val1' + wing]
        mean_data[wing]['min_phi_' + wing] = group_mov.min().reset_index()['maxmin_val1' + wing]
        mean_data[wing]['max_theta_' + wing] = group_mov.max().dropna().reset_index()['theta_' + wing + ' smth']
        mean_data[wing]['min_theta_' + wing] = group_mov.min().dropna().reset_index()['theta_' + wing + ' smth']
        mean_data[wing]['max_psi_' + wing] = group_mov.max().dropna().reset_index()['psi_' + wing + ' smth']
        mean_data[wing]['min_psi_' + wing] = group_mov.min().dropna().reset_index()['psi_' + wing + ' smth']
        return mean_data[wing].dropna()
        
     
     
    def meanfeature(self):
        """
        A method that calculates the mean stroke value of the data
       
        Returns
        -------
        None.
        
        Notes:
        ------
            - This function modifies the 'data' attribute in-place.
       
        """
        wings = ['lw', 'rw']
        mean_data = {}
        # loop on wings, calculate the mean values for each wing
        for wing in ['lw', 'rw']:  
            dropwng = 'lw' if wing == 'rw' else 'rw'
            data_wng = self.data.filter(regex=f'^(?!.*{dropwng}).*$', axis=1) 
            
            if wing == 'lw':
                lw_to_keep = list(filter(lambda string: 'lw' in string ,data_wng.keys())) + ['mov']
                data_wng = data_wng[lw_to_keep]

            grouped_by_mov = data_wng.groupby(['mov','Full_strk' + wing])
            mean_data[wing] = grouped_by_mov.mean().reset_index()              
            mean_data[wing] = self.max_min_wing_angles(mean_data,grouped_by_mov,wing)

        self.mean_data = pd.merge(mean_data['rw'],mean_data['lw'], left_on=['Full_strkrw','mov'], right_on=['Full_strklw','mov'], how='inner')


                   

        

         
    
# =============================================================================
#   Manipulate
# =============================================================================
    def property_projection(self):
        """
        input_project_on_zero_bodaxs (list of ints, optional): Input parameters for properties projection on XY plane. Defaults to [211,3,400].
            input_project_on_zero_bodaxs[0] : Window for savgol_filter to plot projected velocities. 
            input_project_on_zero_bodaxs[1] : polynom degree for savgol_filter to plot projected velocities. 
            input_project_on_zero_bodaxs[2] : Frame of referance axis.
            

        Returns
        -------
        None.

        """
        """
        call data manipulation functions. 
        - project_on_zero_bodaxs: project location, velocity and acceleration to zero bodi axis
        - move_prop_to_zero: substract zero value from property
        

        Returns
        -------
        None.

        """
        self.project_on_zero_bodaxs()
        self.move_prop_to_zero('yaw') # move yaw[0] to zero
        
        
    def project_on_zero_bodaxs(self,input_project_on_zero_bodaxs = [211,3,400]):
        """
        Project the data onto the XY projected zero-body-axis reference.
    
        This method calculates and projects the zero-body-axis and projects the data on it for each movement and axis. The projected data includes position, velocity, and acceleration in both X and Y axes.
    
        Returns
        -------
        None.
        
        Notes:
        ------
            - This function modifies the 'data' attribute in-place.
        """
        
        savgol_win = input_project_on_zero_bodaxs[0]
        savgol_poly = input_project_on_zero_bodaxs[1]
        reframe = input_project_on_zero_bodaxs[2]
        
        xy_on_xy = {'projected_x':[],'projected_x_dot':[],'projected_x_dot_dot':[],'projected_y':[],'projected_y_dot':[],'projected_y_dot_dot':[]}
        dt = (self.data['time [ms]'].iloc[1] - self.data['time [ms]'].iloc[0])/1000
        
        grouped_mov = self.data.groupby('mov')

        for mov_name,movie_group in grouped_mov:      
            # check if zero time exist, use it to calculate the referance axes
            reframe = movie_group.query('`time [ms]` == 0').index if movie_group.query('`time [ms]` == 0').index.size > 0 else reframe
           
            for axis_name in ['x','y']:
                body_axis = self.body_vectors[mov_name][axis_name.upper()][reframe - movie_group.index[0]] # x axis in t = 0/c
                returned_xy_on_xy = self.vector_on_xy(['X smth','Y smth','Z smth'],body_axis,movie_group)
                xy_on_xy['projected_' + axis_name].append(returned_xy_on_xy)
                xy_on_xy['projected_' + axis_name + '_dot'].append(signal.savgol_filter(returned_xy_on_xy/dt, savgol_win, savgol_poly,deriv = 1))
                xy_on_xy['projected_' + axis_name + '_dot_dot'].append(signal.savgol_filter(returned_xy_on_xy/dt**2, savgol_win, savgol_poly,deriv = 2))
  
    
        for column_name, column_data in xy_on_xy.items():
            self.data[column_name] = pd.DataFrame(np.hstack(column_data))
                    
    

    def vector_on_xy(self,props_to_project,new_axis,mov):
        """
        This method calculates the projection of an axis on XY plane and then calculates the projection of 'props_to_project' on this axis

        Parameters
        ----------
        props_to_project (list): A list of the properties to project [px,py,pz].
        new_axis (np.array (1 X 3)): A new axis to calculate the properties on.
        mov  (str): name of movie.

        Returns
        -------
        prop_on_xy : np array
            The property projected on XY plane on the new axis.

        """
        x_project_on_z = new_axis*[0,0,1] # dot product between X body and Z lab multiplied by Zlab
        x_axis_on_xy = (new_axis - x_project_on_z)/LA.norm(new_axis-x_project_on_z ,axis = 1)[np.newaxis].T # project the new axis to XY plane
        prop_on_xy = np.sum(np.multiply(np.tile(x_axis_on_xy,(len(mov),1)),mov[props_to_project].to_numpy()),axis = 1) # dot product of the props and the new projected body vector 
        return prop_on_xy
      
    
    def move_prop_to_zero(self,prop):
        """
        This function substracts the value of property in t = 0 for each movie

        Returns
        -------
        None.
        
        Notes:
        ------
            - This function modifies the 'data' attribute in-place.

        """
        self.data[prop + '_zer'] = self.data.groupby(['mov'])[prop].apply(lambda x: x.sub(x.iloc[0])).reset_index()['yaw']
        
   
    
    def find_zero_cross(self,group,prop,time_min = 50,time_max = 150):
        """
        This method finds the indices where a property change sign. 

        Parameters
        ----------
        group (pandas group): a group from a grouped pandas array (grouped by 'mov')
        prop (str): property to calculate the zero crossing
        time_min (int, optional) : minimal time to search for crossing. The default is 50.
        time_max (int, optional) : maximal time to search for crossing. The default is 150.

        Returns
        -------
        idx_zero_cross (pandas series) : the indices of the zero crossing for each movie

        """
        asign = np.sign(group[prop])
        signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
        idx = np.where(signchange)[0] if len(signchange) > 0 else None     
        zero_cross = group.iloc[idx]
        zero_cross = zero_cross[(zero_cross['time [ms]'] < time_max ) & (zero_cross['time [ms]'] > time_min )]
        idx_zero_cross = zero_cross.index[0] if len(zero_cross.index) > 0 else None
        
        return idx_zero_cross

    def calculate_delta(self,prop,prop_idx,delta_name):
        """
        This method calculates the change in a property (min/max(prop) - prop(t=0))


        Parameters
        ----------
        prop (str) : property to subtract.
        prop_idx (str): property of index to subtract (minimal velocity/ maximal pitch...)
        delta_name (str): name of new column

        Returns
        -------
        None.
        
        Notes:
        ------
            - This function modifies the 'interest_points' attribute in-place.

        """
        
        
        self.interest_points[delta_name] = self.data[prop].iloc[self.interest_points['time_0']].values - self.data[prop].iloc[self.interest_points[prop_idx]].values



    def calculate_interest_idx(self,**kwargs):
        """
        This method runs multiple methods to calculate interest points:
            - 'find_zero_cross' (Change in sign is defined between 'time_max' and 'time_min')
            - 'idxmin' to get the indices of sign change and minimum value of prop. 
            - finds the indices of t = 0
           
       

        Parameters
        ----------
        time_min (float, optional): minimum time to search. The default is 50. (passed to find_zero_cross)
        time_max (float, optional): maximal time to search. The default is 150.

        Returns
        -------
        None.
        
        Notes:
        ------
            - This function modifies the 'idx_zero_cross' attribute in-place.

        """
        
        self.interest_points['minimal_velocity'] = self.data.groupby('mov')['projected_x_dot'].idxmin() # find min to use when there are no sign changes. No time constrain
        self.interest_points['maximal_velocity'] = self.data.groupby('mov')['projected_x_dot'].idxmax() # find min to use when there are no sign changes. No time constrain
        self.interest_points['minimal_acc'] = self.data.groupby('mov')['projected_x_dot_dot'].idxmin() # find min to use when there are no sign changes. No time constrain
        self.interest_points['maximal_acc'] = self.data.groupby('mov')['projected_x_dot_dot'].idxmax() # find min to use when there are no sign changes. No time constrain
     
        self.interest_points['maximal_pitch'] = self.data.groupby('mov')['pitch smth'].idxmax() # find min to use when there are no sign changes. No time constrain
        self.interest_points['zero_velocity']  = self.data.groupby('mov').apply(self.find_zero_cross,'projected_x_dot',**kwargs)  
        self.interest_points['time_0'] = self.data.query('`time [ms]` == 0').index
        
  
 
   
 
    
    @staticmethod
    def get_mov(data,mov_name):
        return data[data['mov'].isin(mov_name)]
         
    
