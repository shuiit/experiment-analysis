import plotly.io as pio
from numpy import linalg 

import numpy as np
from plotly.subplots import make_subplots
import h5py
from scipy.signal import argrelextrema, savgol_filter,find_peaks
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt


pio.renderers.default='browser'

class ProcessH5Movie():
    def __init__(self,movie,mov_name,h5_file):    
        self.mov = {}
        self.data = {dict_name.split('_')[0]:np.array(movie[dict_name]) for dict_name in ['wing_angles','body_angles','vectors_raw']}
        self.header = {dict_name.split('_')[0]:self.get_header(movie[dict_name]) for dict_name in ['wing_angles','body_angles','vectors_raw']}
        h5_file.create_group(f'/{mov_name}/')

        self.name = mov_name
       
        self.dt = np.diff(self.data['body'][:,self.header['body']['time']])[0]/1000
        self.body_savgol_win = [211,511,1011]
        self.body_savgol_poly = [4,2,2]
        self.ref_frame =  np.where(self.data['body'][:,self.header['body']['time']] == 0)[0][0] if len( np.where(self.data['body'][:,self.header['body']['time']] == 0)[0]) > 0 else 0
        
        x_idx = self.header['vectors']['X_x_body']
        if np.isnan(self.data['vectors'][self.ref_frame,:] ).any() == True:
            vectors_exist = np.where(np.isnan(np.sum(self.data['vectors'][: self.ref_frame+50,x_idx:x_idx + 3],axis = 1)) == False)
            self.ref_frame = np.argmin(np.abs(vectors_exist-self.ref_frame))
            frame = self.data['vectors'][self.ref_frame,:] 


        h5_file.attrs['dt'] = self.dt
        h5_file[f'/{mov_name}/'].attrs['ref_frame'] = self.ref_frame

    def get_strokes(self):
        max_stroke_idx = np.hstack([self.wing_stroke(self.data['wing'][:,self.header['wing'][wing_name]]) for wing_name in ['phi_rw','phi_lw']])
        min_stroke_idx = np.hstack([self.wing_stroke(-self.data['wing'][:,self.header['wing'][wing_name]]) for wing_name in ['phi_rw','phi_lw']])
        self.data['wing'] = np.hstack([self.data['wing'],max_stroke_idx,min_stroke_idx])
        self.data['body'] = np.hstack([self.data['body'],max_stroke_idx,min_stroke_idx])
        [self.add_to_header( ['phi_rw_max_idx','phi_rw_max_val','phi_lw_max_idx','phi_lw_max_val',
                              'phi_rw_min_idx','phi_rw_min_val','phi_lw_min_idx','phi_lw_min_val'],dict_name) for dict_name in ['wing','body']]


    def project_body_props_on_xy(self,prop_to_project,projection_name):
        x_idx = self.header['body'][prop_to_project[0]]
        vector_to_project = self.data['body'][:,x_idx:x_idx + 3]
        projected_prop = self.property_projection(vector_to_project)
        savgol_derives,header = self.savgol_and_header(projected_prop.T,['','_dot','_dot_dot'],[projection_name])
        self.data['body'] = np.hstack([self.data['body'],savgol_derives])
        self.add_to_header(header,'body')


    def calculate_angles_frame_ref_axes(self):
        angles = np.array([self.data['body'][:,self.header['body'][angle_name]] for angle_name in ['yaw_body','pitch_body','roll_body']]).T
        rotation_matrices = [self.rotation_matrix(angle[0]*np.pi/180,-angle[1]*np.pi/180,angle[2]*np.pi/180) for angle in angles]
        ref_axes = rotation_matrices[self.ref_frame].T
        rotated_rotation_mat = [np.matmul(ref_axes,rotmat) for rotmat in rotation_matrices]
        angles = np.array([self.angles_body(rotmat) for rotmat in list(rotated_rotation_mat)])
        savgol_derives,header = self.savgol_and_header(angles.T,['','_dot','_dot_dot'],['yaw_z_frame','pitch_y_frame','roll_x_frame'])
        self.data['body'] = np.hstack([self.data['body'],savgol_derives])
        self.add_to_header(header,'body')


    def wing_stroke(self, data):
        peaks, _ = find_peaks(data, prominence=20)
        idx_column = self.create_stroke_column(peaks,data,peaks[0:-1])
        value_column = self.create_stroke_column(peaks,data,-data[peaks[0:-1]])
        return np.vstack((idx_column,value_column)).T
    
    def property_projection(self,prop_to_project):
        frame = self.data['vectors'][self.ref_frame,:] 
        x_idx = self.header['vectors']['X_x_body']
        ref_axes = np.array(frame[x_idx:x_idx + 3])
        return self.body_axes_in_t0(prop_to_project,ref_axes,self.ref_frame)
     
    def body_axes_in_t0(self,prop_to_project,ref_axes,ref_frame_index):
        x_project_on_z = ref_axes * [0,0,1]
        x_axis_on_xy = (ref_axes - x_project_on_z)/linalg.norm(ref_axes-x_project_on_z )[np.newaxis].T # project the new axis to XY plane
        x_axis_on_xy = np.tile(x_axis_on_xy,(len(prop_to_project),1))
        return np.sum(x_axis_on_xy * prop_to_project,axis = 1)[np.newaxis,:].T # dot product of the props and the new projected body vector
     
    def savgol_and_header(self,data,derives,prop_name):
        savgol_derives = np.hstack([savgol_filter(data/(self.dt**deriv), self.body_savgol_win[deriv], self.body_savgol_poly[deriv],deriv = deriv)[np.newaxis,:].T for deriv,name_deriv in enumerate(derives)])
        header = self.add_suffix_to_str(prop_name,derives)
        return np.squeeze(savgol_derives),header

    def calcate_mean_stroke(self,stroke_idx_name,name_to_mean = ['body','wing']):

        stroke_idx = self.data['body'][:,self.header['body'][stroke_idx_name]]
        min_idx = np.unique(stroke_idx[stroke_idx != None]).astype(int)[1:-1]
        {self.data.update({f'mean_{property_name}':self.mean_stroke(min_idx,self.data[property_name])}) for property_name in name_to_mean}
        [self.header.update({f'mean_{property_name}':self.header[property_name]}) for property_name in name_to_mean]
        [self.add_to_header(['mean_idx'],f'mean_{property_name}') for property_name in name_to_mean]


    def add_to_header(self, string_to_add,dict_name):
        [self.header[dict_name].update({name:len(self.header[dict_name])}) for name in string_to_add]



    @ staticmethod
    def mean_stroke(min_idx,data):
        data[data == None] = np.nan
        mean_stroke = [np.nanmean(data[idx0:idx1,:],axis = 0) for idx0,idx1 in zip(min_idx[:-1],min_idx[1:])]
        mean_idx = (min_idx[:-1]+min_idx[1:])/2
        return np.hstack((mean_stroke,np.array(mean_idx)[np.newaxis,:].T))



    @ staticmethod
    def get_header(dataset):
        return { header: idx for idx,header in enumerate(dataset.attrs['header'])}
    
 
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

    @staticmethod
    def add_suffix_to_str(string_in_table,new_suffix):
        header = []
        for suffix in new_suffix:
            header +=  [angle + suffix for angle in  string_in_table ]
        return header
                


        
        