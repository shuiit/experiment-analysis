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
        self.wings = pd.DataFrame(experiment[mov_name]['wing_angles'], columns = experiment[mov_name]['wing_angles'].attrs['header'])
        self.body = pd.DataFrame(experiment[mov_name]['body_angles'], columns = experiment[mov_name]['body_angles'].attrs['header'])
        self.vectors = pd.DataFrame(experiment[mov_name]['raw_vectors'], columns = experiment[mov_name]['raw_vectors'].attrs['header'])
        self.dt = np.diff(self.body['time'])[0]/1000
        self.body_savgol_win = 211
        self.body_savgol_poly = 4



    def get_strokes(self):
        max_stroke_idx = self.wings[['phi_rw','phi_lw']].apply(lambda x: self.wing_stroke(x.to_numpy())).add_suffix('_max_idx')
        min_stroke_idx = self.wings[['phi_rw','phi_lw']].apply(lambda x: self.wing_stroke(-x.to_numpy())).add_suffix('_min_idx')
        self.wings = pd.concat([self.wings,max_stroke_idx,min_stroke_idx],axis = 1)
        self.body = pd.concat([self.body,max_stroke_idx,min_stroke_idx],axis = 1)

    

    def wing_stroke(self, data):
        peaks, _ = find_peaks(data, prominence=20)
        return self.create_stroke_column(peaks,data,peaks[0:-1])
    
    def property_projection(self,prop_to_project):
        ref_frame = self.vectors.query('time == 0') if self.vectors.query('time == 0').index.size > 0 else self.vectors.iloc[0]
        ref_axes = np.array(ref_frame[['X_x_body','X_y_body','X_z_body']])
        return self.body_axes_in_t0(prop_to_project,ref_axes,ref_frame.index)
     

    def project_body_props_on_xy(self,prop_to_project,projection_name):
        self.body[f'{projection_name}'] = self.property_projection(self.body[prop_to_project])
        self.body[f'{projection_name}_dot'] =savgol_filter(self.body[projection_name]/self.dt, self.body_savgol_win, self.body_savgol_poly,deriv = 1)
        self.body[f'{projection_name}_dot_dot'] =savgol_filter(self.body[projection_name]/self.dt**2, self.body_savgol_win, self.body_savgol_poly,deriv = 2)


    def body_axes_in_t0(self,prop_to_project,ref_axes,ref_frame_index):
        x_project_on_z = ref_axes * [0,0,1]
        x_axis_on_xy = (ref_axes - x_project_on_z)/linalg.norm(ref_axes-x_project_on_z ,axis = 1)[np.newaxis].T # project the new axis to XY plane
        x_axis_on_xy = np.repeat(x_axis_on_xy,(len(prop_to_project)),axis = 0)
        prop_to_project = prop_to_project -  np.array(prop_to_project.iloc[ref_frame_index])
        return np.sum(x_axis_on_xy * prop_to_project,axis = 1) # dot product of the props and the new projected body vector 

    @staticmethod
    def create_stroke_column(peaks,data,repeat_value ):
        initial_part_stroke = [None] * peaks[0]
        ending_part_stroke = [None] * (len(data) - peaks[-1])
        full_strokes = np.repeat(repeat_value, np.diff(peaks))
        return np.concatenate([initial_part_stroke, full_strokes, ending_part_stroke])

        
        