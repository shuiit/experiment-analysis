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
from Plotters import Plotters



pio.renderers.default='browser'

class Movie():
    def __init__(self,experiment,mov_name,pertubation = False):    
        self.mov = {}
        self.data = {dict_name:np.array(experiment[mov_name][dict_name]) for dict_name in experiment[mov_name].keys()}
        self.header = {dict_name:self.get_header(experiment[mov_name][dict_name]) for dict_name in experiment[mov_name].keys()}
        self.name = mov_name
        self.ref_frame =experiment[mov_name].attrs['ref_frame']
        self.pertubation = pertubation
        self.pertubation_name = experiment.attrs['dark_pert']
        
    

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
    
    def pqr_pqr_dot(self,angles_data):
        # calculate the body angular acceleratrion and velocity: pqr and pqr_dot

        # calculate the body angular acceleratrion and velocity: pqr and pqr_dot
        pitch=angles_data[:,0];yaw=angles_data[:,1];roll = angles_data[:,2]
        pitch_dot=angles_data[:,3];yaw_dot = angles_data[:,4];roll_dot = angles_data[:,5]
        pitch_dot_dot=angles_data[:,6];yaw_dot_dot = angles_data[:,5];roll_dot_dot = angles_data[:,6]
        
        p = roll_dot-yaw_dot*np.sin(np.radians(pitch))
        q = pitch_dot*np.cos(np.radians(roll))+yaw_dot*np.sin(np.radians(roll))*np.cos(np.radians(pitch));
        r = -pitch_dot*np.sin(np.radians(roll))+yaw_dot*np.cos(np.radians(roll))*np.cos(np.radians(pitch));
        
        
        p_dot = roll_dot_dot-yaw_dot_dot*np.sin(np.radians(pitch))-yaw_dot*pitch_dot*np.cos(np.radians(pitch));
        q_dot = (pitch_dot_dot*np.cos(np.radians(roll))-pitch_dot*np.sin(np.radians(roll))*roll_dot
                        +yaw_dot_dot*np.sin(np.radians(roll)*np.cos(pitch)+
                        +yaw_dot*np.cos(np.radians(roll))*np.cos(np.radians(pitch))*roll_dot
                        -yaw_dot*np.sin(np.radians(roll))*np.sin(np.radians(pitch))*pitch_dot))
                    
                    
        r_dot = (-pitch_dot_dot*np.sin(np.radians(roll))-pitch_dot*np.cos(np.radians(roll))*roll_dot
                        +yaw_dot_dot*np.cos(np.radians(roll))*np.cos(np.radians(pitch))
                        -yaw_dot*np.sin(np.radians(roll))*np.cos(np.radians(pitch))*roll_dot
                        -yaw_dot*np.cos(np.radians(roll))*np.sin(np.radians(pitch))*pitch_dot)
        return np.vstack((p,q,r,p_dot,q_dot,r_dot)).T

    def calculate_pqr_update_data_header(self):
        angles_data = self.data['body'][:,self.header['body']['pitch_body']:self.header['body']['roll_body_dot_dot']]
        pqr_pqr_dot_data = self.pqr_pqr_dot(angles_data)

        self.data['body'] = np.hstack((self.data['body'], pqr_pqr_dot_data))
        pqr_header = [pqr + deriv for deriv in ['','_dot','_dot_dot'] for pqr in  ['p','q','r']]
        self.add_to_header(pqr_header,'body')


    def plot_3d_traj_movie(self,color_prop):
        data,plot_cofnig = self.calculation_for_3d_traj(color_prop = color_prop)
        ploter = Plotters(data,False,self.pertubation)
        return ploter.plot_3d_traj(data,plot_cofnig,self.name,self.pertubation_name,color_prop = color_prop )
        

    def get_prop(self,prop,wing_body):
        return self.data[wing_body][:,self.header[wing_body][prop]]
    
    def get_idx_of_time(self,t):
        time = self.get_prop('time','body')
        return np.where(time == t)[0]

    

    def mean_prop_stroke(self,prop_name,wing_body_prop,wing_body_mean_save,phi_idx_to_mean = 'phi_rw_min_idx'):
        data = self.get_prop(prop_name,wing_body_prop)
        stroke_idx = self.get_prop(phi_idx_to_mean,'wing')
        mean_strok_idx = self.header[wing_body_mean_save]['mean_idx']
        min_idx = np.unique(stroke_idx[stroke_idx != None]).astype(int)[1:-1]
        val = np.intersect1d((min_idx[:-1]+min_idx[1:])/2, self.data[wing_body_mean_save][:,mean_strok_idx], assume_unique=False, return_indices=True)
        data[data == None] = np.nan
        mean_stroke = [np.nanmean(data[idx0:idx1],axis = 0) for idx0,idx1 in zip(min_idx[:-1],min_idx[1:])]
        self.data[f'{wing_body_mean_save}'] = np.hstack((self.data[f'{wing_body_mean_save}'],np.array(mean_stroke)[val[1]][np.newaxis,:].T))
        self.add_to_header([f'{prop_name}'],f'{wing_body_mean_save}')

    def plot_prop(self,prop,wing_body,color,name,fig,prop_x = 'time',t0 = False,t1 = False,**kwargs):

        t0_idx =  self.get_idx_of_time(t0) if (t0 != False & len(self.get_idx_of_time(t0))>0) else [0]
        t1_idx =  self.get_idx_of_time(t1) if (t1 != False & len(self.get_idx_of_time(t1))>0) else [int(-1)]

        data_y = self.get_prop(prop,wing_body)
        data_x = self.get_prop(prop_x,wing_body)
        ploter = Plotters(data_y,False,self.pertubation)
        return ploter.plot_prop_movie_v2(data_x[t0_idx[0]:t1_idx[0]],data_y[t0_idx[0]:t1_idx[0]],color,name,fig = fig,**kwargs)

    
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

            

        
    
