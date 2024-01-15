import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px
import matplotlib.cm
from numpy import linalg 

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from plotly.subplots import make_subplots
import h5py
from scipy.signal import argrelextrema, savgol_filter,find_peaks
from scipy.spatial.transform import Rotation as R



pio.renderers.default='browser'

class ManipulatedMovie():
    def __init__(self,experiment,mov_name,pertubation = False):    
        self.mov = {}
        self.data = {dict_name:np.array(experiment[mov_name][dict_name]) for dict_name in experiment[mov_name].keys()}
        self.header = {dict_name:self.get_header(experiment[mov_name][dict_name]) for dict_name in experiment[mov_name].keys()}
        self.name = mov_name
        self.ref_frame =experiment[mov_name].attrs['ref_frame']
        self.pertubation = pertubation
    

    def get_header(self,dataset):
        return { header: idx for idx,header in enumerate(dataset.attrs['header'])}
    
    
    def add_to_header(self, string_to_add,dict_name):
        [self.header[dict_name].update({name:len(self.header[dict_name])}) for name in string_to_add]


    def min_max_point(self,prop,wing_body = 'body'):
        idx_time = self.header[wing_body]['time']
        idx_prop = self.header[wing_body][prop]

        time_min_v = self.data[wing_body][np.argmin(self.data[wing_body][:,idx_prop]),idx_time]
        time_max_pitch = self.data[wing_body][np.argmax(self.data[wing_body][:,idx_prop]),idx_time]
        return np.vstack((time_min_v,time_max_pitch))
    
    def zero_velocity(self,prop,wing_body = 'body'):
        
        idx_time = self.header[wing_body]['time']
        idx_prop = self.header[wing_body][prop]

        zero_v = np.where(np.diff(np.sign(self.data[wing_body][:,idx_prop]))<0)[0]
        if len(zero_v) > 0:
            return self.data[wing_body][zero_v,idx_time]
        
    
    def mean_props(self,prop1,prop2,wing_body,header_name):
        mean_prop = (self.get_prop(prop1,wing_body)  + self.get_prop(prop2,wing_body) )/2
        self.data[wing_body] = np.hstack((self.data[wing_body], mean_prop[np.newaxis,:].T))
        self.add_to_header([header_name],wing_body)

    def calculation_for_3d_traj(self, color_prop = 'pitch'):
        data = {}
        plot_cofnig = {'fly_samples':150,'traj_samples':20,'size_x':1,'size_y':1/3,'delta_y_on_x':3/4}

        vectors = {prop_name.split('_')[0] : self.data['vectors'][:,self.header['vectors'][prop_name]:self.header['vectors'][prop_name] + 3] for prop_name in ['X_x_body','Y_x_body','Z_x_body']}
        cm_idx = self.header['body']['CM_real_x_body']
        pitch_idx = self.header['body'][color_prop]

        data['time'] = self.get_prop('time','body')

        data['cm'] = self.data['body'][:,cm_idx:cm_idx + 3]*1000
        data[color_prop] = self.data['body'][:,pitch_idx]
        body_x_vector = vectors['X']*plot_cofnig['size_x'] + data['cm'] # define the size of Xbody vecotr
        body_y_vector = vectors['Y'][::plot_cofnig['fly_samples'],:]*plot_cofnig['size_y']
        delta_y_on_x = (vectors['X'][::plot_cofnig['fly_samples'],:])*plot_cofnig['delta_y_on_x'] + data['cm'][::plot_cofnig['fly_samples'],:] # define the size of Ybody vecotr
        idx_end_pertubation = np.where((data['time'] <(self.pertubation + 1)) & (data['time'] >(self.pertubation - 1)) )[0][0] if self.pertubation != False else False
        
        
        data['start_pert_endpert'] = [0,self.ref_frame,idx_end_pertubation] if self.pertubation != False else [0,self.ref_frame]
        data['x_vector'] = self.disconnect_line_add_none(data['cm'][::plot_cofnig['fly_samples'],:],body_x_vector[::plot_cofnig['fly_samples'],:])
        data['y_vector'] = self.disconnect_line_add_none(-body_y_vector + delta_y_on_x,body_y_vector+delta_y_on_x)
        return data,plot_cofnig


        

    def get_prop(self,prop,wing_body):
        return self.data[wing_body][:,self.header[wing_body][prop]]
    

    
    @staticmethod
    def disconnect_line_add_none(array1,array2):
        """
        Combines two 2D arrays such that the rows alternate between the two arrays,
        and inserts `None` (represented as `np.nan`) every two rows.

        Args:
            array1 (numpy.ndarray): The first input array with shape (n, 3).
            array2 (numpy.ndarray): The second input array with shape (n, 3).

        Returns:
            numpy.ndarray: The combined array with alternating rows from `array1` and `array2`,
            and `None` (np.nan) inserted every two rows.

        """
        combined_array = np.empty((2 * array1.shape[0], array1.shape[1]))
        combined_array[::2] =array1[:array1.shape[0]]
        combined_array[1::2] =  array2[:array1.shape[0]]
        combined_array = np.insert(combined_array,range(2, combined_array.shape[0], 2),np.nan,axis = 0)
        return combined_array

            

        
    
