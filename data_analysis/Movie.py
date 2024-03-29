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
import Plotters 



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
        self.dt = experiment.attrs['dt']
        
        self.smoothing_config = {'smooth_window_body_angles': 73*7,'smooth_window_body_cm': 73*7, 'smooth_poly_body':3, 'smooth_window_wing': 15}
        
    

    def get_header(self,dataset):
        return { header: idx for idx,header in enumerate(dataset.attrs['header'])}
    
    
    def add_to_header(self, string_to_add,dict_name):
        [self.header[dict_name].update({name:len(self.header[dict_name])}) for name in string_to_add]




    def smooth_and_derive(self,prop,derivs,wing_body,three_col = 3):
        idx = self.header[wing_body][prop]
        header = [list(self.header['body'].keys())[list(self.header['body'].values()).index(idx)] for idx in range(idx,idx + three_col)]
        data = self.get_prop(prop,wing_body,three_col=3)
        
        header = [f'{prop}{derive}_smth' for derive in derivs for prop in header]
        smoothed = np.vstack([savgol_filter(data.T/self.dt**deriv,self.smoothing_config['smooth_window_body_cm'],self.smoothing_config['smooth_poly_body'],deriv = deriv) for deriv,der in enumerate(derivs)]).T
        self.data[wing_body] = np.hstack((self.data[wing_body], smoothed))
        self.add_to_header(header,wing_body)

    def sub_two_props(self,prop1,prop2,wing_body,header):
        prop1 = self.get_prop(prop1,wing_body)
        prop2 = self.get_prop(prop2,wing_body)
        self.data[wing_body] = np.hstack((self.data[wing_body], prop1-prop2))
        self.add_to_header(header,wing_body)



    def sub_ref_frame(self,prop,wing_body):
        prop_to_sub = self.get_prop(prop,wing_body)
        sub_prop = prop_to_sub - prop_to_sub[self.ref_frame,:] 
        self.data[wing_body] = np.hstack((self.data[wing_body], sub_prop))
        self.add_to_header([f'{prop}_min_ref_frame'],wing_body)




    def v_size(self,prop_x,prop_y,wing_body):
        xv = self.get_prop(prop_x,wing_body)
        yv = self.get_prop(prop_y,wing_body)

        sub_prop = np.sqrt(xv**2 + yv**2)
        self.data[wing_body] = np.hstack((self.data[wing_body], sub_prop))
        self.add_to_header([f'amp_v'],wing_body)



    def t0_t1_idx(self,t0,t1):
        time = self.get_prop('time','body')[:,0]
        idx_t1 = np.where(t1 == time)[0][0] if t1 != -1 else -1
        idx_t0 = np.where(t0 == time)[0][0]
        
        return idx_t0,idx_t1
    
    def project_axes_xy(self,header,axis = 'X_x_body',wing_body = 'vectors',add_to_vectors= False):
        
        data = self.get_prop(axis,wing_body, three_col= 3)
        data_axis_on_xy = (data[:,0:2])/np.linalg.norm(data[:,0:2],axis = 1)[np.newaxis].T # project the new axis to XY plane
        self.data[wing_body] = np.hstack((self.data[wing_body], data_axis_on_xy))
        self.add_to_header(header,wing_body)
        
        if add_to_vectors != False:
            [self.from_wing_body_to_vectors(head,add_to_vectors[0],add_to_vectors[1]) for head in header ]



    def project_prop_all_axes(self,prop,wing_body = 'body',header_name = 'CM_dot',ax_to_proj = 'X_x_body',add_to_vectors= False,three_col = 3):
        data = self.get_prop(prop,wing_body, three_col= three_col) 
    
        projected_axes = self.get_prop(ax_to_proj,wing_body,three_col = three_col)


        projected = np.sum(projected_axes * data,axis = 1)[np.newaxis,:].T

        self.data[wing_body] = np.hstack((self.data[wing_body], projected))
        self.add_to_header([f'{header_name}_projected_all_axes'],wing_body)
        if add_to_vectors != False:
            self.from_wing_body_to_vectors(f'{header_name}_projected_all_axes',add_to_vectors[0],add_to_vectors[1])
    
        
    def project_prop(self,prop,wing_body = 'body',header_name = 'CM_dot',ax_to_proj = 'X_x_body',add_to_vectors= False,three_col = 3):
        data = self.get_prop(prop,wing_body, three_col= three_col) 
        projected_axes = self.get_prop(ax_to_proj,'vectors',three_col = three_col)[self.ref_frame,:]
        projected = np.sum(np.tile(projected_axes,(len(data),1)) * data,axis = 1)[np.newaxis,:].T

        self.data[wing_body] = np.hstack((self.data[wing_body], projected))
        self.add_to_header([f'{header_name}_projected'],wing_body)
        if add_to_vectors == True:
            self.from_wing_body_to_vectors(f'{header_name}_projected','body','vectors')
    
    def calc_angle_unwarp(self, prop,delta_ang = 180):
        
        sgn = np.sign(np.cross(prop[0],prop[1]))[:,-1]
        ang_mov  = sgn*np.arccos(np.sum(prop[0] * prop[1],axis = 1))*180/np.pi
        idx = np.where(np.isnan(ang_mov))[0]
        ang_mov[idx] = 0
        unwarped = np.unwrap(ang_mov + delta_ang,period = 360) - delta_ang 
        unwarped[idx] = np.nan
        return unwarped

    def delta_ang_all_time(self,prop1_name,prop2_name,header,three_col = 2,**kwargs):
        
        
        prop = [np.insert(self.get_prop(prop_name,'vectors',three_col=three_col),2,0,axis = 1) for prop_name in [prop1_name,prop2_name]]
        unwarped = self.calc_angle_unwarp(prop,**kwargs)

        self.data['vectors'] = np.vstack((self.data['vectors'].T, unwarped)).T
        self.add_to_header( [header],'vectors')
        self.from_wing_body_to_vectors(header,'vectors', 'body')


    def from_wing_body_to_vectors(self,prop,from_wbv, to_wbv):

        prop_to_add = self.get_prop(prop,from_wbv)
        frame_wing_body = self.get_prop('frames',from_wbv)
        frames_vector = self.get_prop('frames',to_wbv)
        rows_frames,frames_vector_ind, frame_body_ind = np.intersect1d(frames_vector,frame_wing_body, return_indices=True )
        nan_row = np.full([len(frames_vector)], np.nan)
        nan_row[frames_vector_ind] = prop_to_add[frame_body_ind,:].T
        self.data[to_wbv] = np.vstack((self.data[to_wbv].T,nan_row)).T
        self.add_to_header( [prop],to_wbv)
       
    def delta_ang_ref_frame(self,ref_frame_axis,header,axis = 'X_x_body_projected',three_col = 2,**kwargs):
        
        data_axis = self.get_prop(axis,'vectors', three_col= three_col)
        prop_norm = data_axis/np.linalg.norm(data_axis,axis = 1, keepdims=True)
        ref_axis = self.get_prop(ref_frame_axis,'vectors', three_col= three_col)[self.ref_frame,:]
        ref_axis = ref_axis/np.linalg.norm(ref_axis)
        prop = [np.insert(prop_norm,2,0,axis = 1),
                np.repeat([np.insert(ref_axis,2,0)],len(data_axis),axis = 0)]

        unwarped = self.calc_angle_unwarp(prop,**kwargs)
        self.data['vectors'] = np.vstack((self.data['vectors'].T, unwarped)).T
        self.add_to_header([header],'vectors')
        self.from_wing_body_to_vectors(header,'vectors','body')

    def get_min(self,prop,t1 = False,t0 = False):
        
        prop_to_min = self.get_prop(prop,'body')[:,0]
        idx_time = (0,-1 )if (t1 == False) | (t0 == False) else self.t0_t1_idx(t0,t1)
        min_idx = np.argmin(prop_to_min[idx_time[0]:idx_time[1]], prominence=0.1 )
        if len(max) > 0:
            time = self.get_prop('time','body')
            return time[min_idx + idx_time[0]]
        
    def add_time_m_t0(self, t0,header,wing_body = 'body'):
        time = self.get_prop('time','body')
        self.data[wing_body] = np.hstack((self.data[wing_body], time - t0))
        self.add_to_header([header],wing_body)
         
    def get_delta_prop(self,prop,wing_body, t_fin = 230, delta_frames = 300, time_prop = 'time'):
        
        prop = self.get_prop(prop,wing_body)
        time = self.get_prop(time_prop,wing_body)
        tend = np.where(time > t_fin)

        if (time > t_fin).any() == True:
            delta_v_ini = np.nanmean(prop[:delta_frames])
            delta_v_fin = np.nanmean(prop[tend[0][0]:tend[0][0] + delta_frames])
            return np.abs(delta_v_ini - delta_v_fin)

    def get_peaks_min_max(self,prop,case,t1 = False,t0 = False):

        idx_time = (0,-1 )if (t1 == False) | (t0 == False) else self.t0_t1_idx(t0,t1)
        acc = self.get_prop(prop,'body')[idx_time[0]:idx_time[1],0]
        v0 = acc[self.ref_frame]/2

        if case == 'peaks_max':
            idx = find_peaks(acc, prominence=0.1 )[0]
        if case == 'peaks_min':
            idx = find_peaks(-acc, prominence=0.05 )[0]
        if case == 'min':
            idx = [np.argmin(acc)]
        if case == 'max':
            idx = [np.argmax(acc)]
        if case == 'zero_v':
            idx =  np.where((np.diff(np.sign(acc))<0) | (np.diff(np.sign(-acc))<0))[0]

        if case == 'half_v':
            idx =  np.where((np.diff(np.sign(acc - v0))<0) | (np.diff(np.sign(-acc - v0))<0))[0]

        if len(idx) > 0:
            time = self.get_prop('time','body')[:,0]
            return self.data['body'][idx[0]+ idx_time[0],:]
        else:
            return self.data['body'][0,:]*np.nan

    def add_point_to_plot(self,interest_points,ydata,fig,color,wing_body = 'body',xdata = 'time',**kwargs):   
        idx_y = self.header[wing_body][ydata] 
        idx_x = self.header[wing_body][xdata] 
        Plotters.add_point_to_plot(interest_points[idx_x],interest_points[idx_y],self.name,color,fig,**kwargs) 


    def calc_force(self):
        
        ay = self.get_prop('CM_dot_dot_y_projected_all_axes','body')
        ax = self.get_prop('CM_dot_dot_x_projected_all_axes','body')
        az = self.get_prop('CM_real_z_body_dot_dot_smth','body')

        force = np.hstack((ax,ay,az))
        force_norm = (force.T/np.linalg.norm(force,axis = 1)).T
        force_size = np.sqrt(np.sum(force**2,axis = 1))

        self.data['body'] = np.vstack((self.data['body'].T, force_size)).T
        self.add_to_header(['force_size'],'body')
    
    def mean_props(self,prop1,prop2,wing_body,header_name):
        mean_prop = (self.get_prop(prop1,wing_body)  + self.get_prop(prop2,wing_body) )/2
        self.data[wing_body] = np.hstack((self.data[wing_body], mean_prop))
        self.add_to_header([header_name],wing_body)

    def calculation_for_3d_traj(self, color_prop = 'pitch',plot_cofnig = {'fly_samples':150,'traj_samples':20,'size_x':1,'size_y':1/3,'delta_y_on_x':3/4},forces = None):
        """Calulations for ploting a 3d trajectory

        Args:
            color_prop (str, optional):property to color the cm .
            plot_cofnig (dict, optional):fly_samples - delta sample of the fly axes .
                                        traj_samples - delta samle of cm 
                                        size_x - scale of body axis
                                        size_y - scale of y axis
                                        delta_y_on_x - location of y axis on x axis (make it look like a cross)

        Returns:
            data (dict): a dictionary with all relevant data
            plot_cofig (dict) : the configuration of the plot
        """
        data = {}

        vectors = {prop_name.split('_')[0] : self.get_prop(prop_name,'vectors',three_col=3) for prop_name in ['X_x_body','Y_x_body','Z_x_body']}
        data['cm'] = self.get_prop('CM_real_x_body','body',three_col=3)*1000
        data['time'] = self.get_prop('time','body')
        if isinstance(color_prop,str):
            data[color_prop] =  self.get_prop(color_prop,'body')[:,0]  
        else:
            data.update(color_prop)
        # define body x and y vecotrs
        body_x_vector = vectors['X']*plot_cofnig['size_x'] + data['cm']
        body_y_vector = vectors['Y'][::plot_cofnig['fly_samples'],:]*plot_cofnig['size_y']
        delta_y_on_x = (vectors['X'][::plot_cofnig['fly_samples'],:])*plot_cofnig['delta_y_on_x'] + data['cm'][::plot_cofnig['fly_samples'],:] # define the size of Ybody vecotr
        
        # disconnect the coordinates to get line for the vectors
        data['x_vector'] = self.disconnect_line_add_none(data['cm'][::plot_cofnig['fly_samples'],:],body_x_vector[::plot_cofnig['fly_samples'],:])
        data['y_vector'] = self.disconnect_line_add_none(-body_y_vector + delta_y_on_x,body_y_vector+delta_y_on_x)       
        if forces.all() != None:
            forces = forces*plot_cofnig['size_x'] + data['cm']
            data['forces'] = self.disconnect_line_add_none(data['cm'][::plot_cofnig['fly_samples'],:],forces[::plot_cofnig['fly_samples'],:])

        # index for pertubationmarker
        idx_end_pertubation = np.where((data['time'] <(self.pertubation + 1)) & (data['time'] >(self.pertubation - 1)) )[0][0] if self.pertubation != False else False
        data['start_pert_endpert'] = [0,self.ref_frame,idx_end_pertubation] if self.pertubation != False else [0,self.ref_frame]

        return data,plot_cofnig
    
    def pqr_pqr_dot(self,angles_data):
        # calculate the body angular acceleratrion and velocity: pqr and pqr_dot

        # calculate the body angular acceleratrion and velocity: pqr and pqr_dot
        pitch=-angles_data[:,0];yaw=angles_data[:,1];roll = angles_data[:,2]
        pitch_dot=-angles_data[:,3];yaw_dot = angles_data[:,4];roll_dot = angles_data[:,5]
        pitch_dot_dot=-angles_data[:,6];yaw_dot_dot = angles_data[:,5];roll_dot_dot = angles_data[:,6]
        
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
        pqr_header = [pqr + deriv for deriv in ['','_dot'] for pqr in  ['p','q','r']]
        self.add_to_header(pqr_header,'body')


    def plot_3d_traj_movie(self,color_prop,**kwargs):
        data,plot_cofnig = self.calculation_for_3d_traj(color_prop = color_prop,**kwargs)
        color_prop = list(color_prop.keys())[0] if isinstance(color_prop,dict) else color_prop
        return Plotters.plot_3d_traj(data,plot_cofnig,self.name,self.pertubation_name,color_prop = color_prop )
        

    def get_prop(self,prop,wing_body,three_col = 1):
        return self.data[wing_body][:,self.header[wing_body][prop]:self.header[wing_body][prop] + three_col].copy()
    
    def get_idx_of_time(self,t,wing_body = 'body'):
        time = self.get_prop('time',wing_body)
        return np.where(time >= t)[0] 

    
    def get_mean_prop_time(self,t0,t1,prop,wing_body,mean_delta,freq_ol = 280):
        
        freq = self.get_prop(prop,wing_body)
        idx_time = self.mean_time([t0,t1],wing_body, time = 'time')
        data_to_mean = freq[idx_time[0]:idx_time[1]]
        data_to_mean[data_to_mean>freq_ol] = np.nan
        return  np.nanmean(data_to_mean) - mean_delta,t0
    
    def mean_prop_time_vector(self,prop,delta_t,t_vec,sub_mean = True,**kwargs):
        
        mean_delta = self.get_mean_prop_time(0,0 + delta_t,prop,'mean_body',0,**kwargs) if sub_mean == True else [0] 
        return [self.get_mean_prop_time(t0,t0 + delta_t,prop,'mean_body',mean_delta[0],**kwargs) for t0 in t_vec[1:]]

    def mean_time(self,t0_vec,wing_body, time = 'time'):
        time = self.get_prop(time,wing_body)
        return [np.argmin(np.abs(time - t0)) for t0 in t0_vec]


    # def mean_prop_stroke(self,prop_name,wing_body_prop,wing_body_mean_save,phi_idx_to_mean = 'phi_rw_min_idx'):
    #     data = self.get_prop(prop_name,wing_body_prop)
    #     stroke_idx = self.get_prop(phi_idx_to_mean,'wing')
    #     mean_strok_idx = self.header[wing_body_mean_save]['mean_idx']
    #     min_idx = np.unique(stroke_idx[stroke_idx != None]).astype(int)[1:-1]
    #     val = np.intersect1d((min_idx[:-1]+min_idx[1:])/2, self.data[wing_body_mean_save][:,mean_strok_idx], assume_unique=False, return_indices=True)
    #     data[data == None] = np.nan
    #     mean_stroke = [np.nanmean(data[idx0:idx1],axis = 0) for idx0,idx1 in zip(min_idx[:-1],min_idx[1:])]
    #     self.data[f'{wing_body_mean_save}'] = np.hstack((self.data[f'{wing_body_mean_save}'],np.array(mean_stroke)[val[1]][np.newaxis,:].T))
    #     self.add_to_header([f'{prop_name}'],f'{wing_body_mean_save}')
 
    def plot_prop(self,prop,wing_body,color,name,fig,prop_x = 'time',t0 = False,t1 = False,**kwargs):

        t0_idx =  self.get_idx_of_time(t0,wing_body) if ((t0 != False) and (len(self.get_idx_of_time(t0,wing_body))>0)) else [0]
        t1_idx =  self.get_idx_of_time(t1,wing_body) if ((t1 != False) and (len(self.get_idx_of_time(t1,wing_body))>0)) else [int(-1)]

        data_y = self.get_prop(prop,wing_body)
        data_x = self.get_prop(prop_x,wing_body)
        
        return Plotters.plot_prop_movie(data_x[t0_idx[0]:t1_idx[0],0],data_y[t0_idx[0]:t1_idx[0],0],color,name,fig = fig,**kwargs)

    def calculate_freq(self, idx_prop, mean_wing_body):
        prop_idx = self.get_prop(idx_prop,mean_wing_body)
        self.data[mean_wing_body] = np.hstack((self.data[mean_wing_body],np.vstack([1/(np.diff(prop_idx.T)*1/16000).T,[np.nan]])))
        self.add_to_header([f'freq_{idx_prop}'],mean_wing_body)
    
    def calculate_phi_amp(self,wing,mean_wing_body):
        min_idx = self.get_prop(f'phi_{wing}_min_val',mean_wing_body)
        max_idx = self.get_prop(f'phi_{wing}_max_val',mean_wing_body)

        self.data[mean_wing_body] = np.hstack((self.data[mean_wing_body],max_idx - min_idx))
        self.add_to_header([f'amp_{wing}'],mean_wing_body)


    def mean_by_stroke(self,prop,mean_wing_body,wing_body):
        data = self.get_prop(prop,wing_body)
        min_idx = self.get_prop('phi_rw_min_idx',wing_body)
        mean_idx = self.get_prop('phi_rw_min_idx',mean_wing_body)
        self.data[mean_wing_body] = np.vstack((self.data[mean_wing_body].T,[np.nanmean(data[min_idx == val]) for val in mean_idx if np.isnan(val) == False])).T
        self.add_to_header([prop],mean_wing_body)


    
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

            

        
    
