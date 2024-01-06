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

class Movie():
    def __init__(self,experiment,mov_name):    
        self.mov = {}
        self.data = {dict_name.split('_')[0]:np.array(experiment[mov_name][dict_name]) for dict_name in ['wing_angles','body_angles','raw_vectors']}
        self.header = {dict_name.split('_')[0]:self.get_header(experiment[mov_name][dict_name]) for dict_name in ['wing_angles','body_angles','raw_vectors']}

        self.data['vectors'] = self.data.pop('raw')
        self.header['vectors'] = self.header.pop('raw')

        self.name = mov_name
       
        self.dt = np.diff(self.data['body'][:,self.header['body']['time']])[0]/1000
        self.body_savgol_win = 211
        self.body_savgol_poly = 4
        self.ref_frame =  np.where(self.data['body'][:,self.header['body']['time']] == 0)[0][0] if len( np.where(self.data['body'][:,self.header['body']['time']] == 0)[0]) > 0 else 0
    
    def get_header(self,dataset):
        return { header: idx for idx,header in enumerate(dataset.attrs['header'])}
    
    
    def add_to_header(self, string_to_add,dict_name):
        [self.header[dict_name].update({name:len(self.header[dict_name])}) for name in string_to_add]

    def get_strokes(self):

        max_stroke_idx = np.array([self.wing_stroke(self.data['wing'][:,self.header['wing'][wing_name]]) for wing_name in ['phi_rw','phi_lw']]).T
        min_stroke_idx = np.array([self.wing_stroke(-self.data['wing'][:,self.header['wing'][wing_name]]) for wing_name in ['phi_rw','phi_lw']]).T
        self.data['wing'] = np.hstack([self.data['wing'],max_stroke_idx,min_stroke_idx])
        self.data['body'] = np.hstack([self.data['body'],max_stroke_idx,min_stroke_idx])
        [self.add_to_header( ['phi_rw_max_idx','phi_lw_max_idx','phi_rw_min_idx','phi_lw_min_idx'],dict_name) for dict_name in ['wing','body']]

    def wing_stroke(self, data):
        peaks, _ = find_peaks(data, prominence=20)
        return self.create_stroke_column(peaks,data,peaks[0:-1])
    
    def property_projection(self,prop_to_project):
        frame = self.data['vectors'][self.ref_frame,:]
        x_idx = self.header['vectors']['X_x_body']
        ref_axes = np.array(frame[x_idx:x_idx + 3])
        return self.body_axes_in_t0(prop_to_project,ref_axes,self.ref_frame)
     

    def project_body_props_on_xy(self,prop_to_project,projection_name):
        x_idx = self.header['body'][prop_to_project[0]]
        vector_to_project = self.data['body'][:,x_idx:x_idx + 3]

        self.data['body'] = np.hstack([self.data['body'],self.property_projection(vector_to_project)])
        self.data['body'] = np.hstack([self.data['body'],savgol_filter(self.data['body'][:,-1]/self.dt, self.body_savgol_win, self.body_savgol_poly,deriv = 1)[np.newaxis,:].T])
        self.data['body'] = np.hstack([self.data['body'],savgol_filter(self.data['body'][:,-1]/self.dt**2, self.body_savgol_win, self.body_savgol_poly,deriv = 2)[np.newaxis,:].T])
        self.add_to_header( [f'{projection_name}',f'{projection_name}_dot',f'{projection_name}_dot_dot'],'body')

    def calculate_angles_frame_ref_axes(self):
        angles = np.array([self.data['body'][:,self.header['body'][angle_name]] for angle_name in ['yaw_body','pitch_body','roll_body']]).T
        rotation_matrices = [self.rotation_matrix(angle[0]*np.pi/180,-angle[1]*np.pi/180,angle[2]*np.pi/180) for angle in angles]
        ref_axes = rotation_matrices[self.ref_frame].T
        rotated_rotation_mat = [np.matmul(ref_axes,rotmat) for rotmat in rotation_matrices]
        angles = [self.angles_body(rotmat) for rotmat in list(rotated_rotation_mat)]
        self.data['body'] = np.hstack([self.data['body'],np.array(angles)])
        self.add_to_header( ['yaw_z_frame','pitch_y_frame','roll_x_frame'],'body')



    
    def body_axes_in_t0(self,prop_to_project,ref_axes,ref_frame_index):
        x_project_on_z = ref_axes * [0,0,1]
        x_axis_on_xy = (ref_axes - x_project_on_z)/linalg.norm(ref_axes-x_project_on_z )[np.newaxis].T # project the new axis to XY plane
        x_axis_on_xy = np.tile(x_axis_on_xy,(len(prop_to_project),1))
        prop_to_project = prop_to_project -  np.array(prop_to_project[ref_frame_index,:])
        return np.sum(x_axis_on_xy * prop_to_project,axis = 1)[np.newaxis,:].T # dot product of the props and the new projected body vector 
    
    @staticmethod
    def angles_body(dcm):
        """calculate the yaw, pitch and roll angles that correspond to each rotation matrix in BODY AXES
        ! if you define new rotation order/ axes, the calculation will be different

        Args:
            dcm (np.array): rotation matrix in camera axes

        Returns:
            yaw_z,pitch_y,roll_x (float): angles of rotation
        """

        yaw_z = np.arctan2(dcm[1,0],dcm[0,0])*180/np.pi
        pitch_y = np.arcsin(dcm[2,0])*180/np.pi
        roll_x = np.arctan2(dcm[2,1],dcm[2,2])*180/np.pi
        return yaw_z,pitch_y,roll_x
        
    @staticmethod
    
    def rotation_matrix(yaw,pitch,roll):
        roll_mat = np.array([[1,0,0],[0 ,np.cos(roll),-np.sin(roll)],[0, np.sin(roll), np.cos(roll)]])
        pitch_mat = np.array([[np.cos(pitch),0,np.sin(pitch)],[0, 1,0],[-np.sin(pitch), 0, np.cos(pitch)]])
        yaw_mat = np.array([[np.cos(yaw),-np.sin(yaw),0],[np.sin(yaw),np.cos(yaw),0],[0, 0, 1]])
        return yaw_mat @ pitch_mat @ roll_mat

    @staticmethod
    def create_stroke_column(peaks,data,repeat_value ):
        initial_part_stroke = [None] * peaks[0]
        ending_part_stroke = [None] * (len(data) - peaks[-1])
        full_strokes = np.repeat(repeat_value, np.diff(peaks))
        return np.concatenate([initial_part_stroke, full_strokes, ending_part_stroke])

        
        