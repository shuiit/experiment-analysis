import pandas as pd
import matplotlib.pyplot as plt
import h5py
import os
import sys
from scipy.signal import savgol_filter
import numpy as np

class PreProcess:

    def __init__(self,exp_dir,exp_name,pert_for_attr, smooth_window_body = 211,smooth_poly_body = 4,smooth_window_wing = 15 ,
                 smooth_poly_wing = 2,zscore_window = 200,zscore_threashold = 3,interp_method = 'cubic'):
        self.exp_dir = exp_dir
        self.exp_name = exp_name
        self.exp_path = f'{exp_dir}/{exp_name}'
        self.pert_for_attr = pert_for_attr
        self.exp_file  = h5py.File(f'{exp_dir}/{exp_name}.hdf5', "a")
        self.exp_file.attrs['dark_pert'] = self.pert_for_attr


        self.smooth_window_body = smooth_window_body
        self.smooth_poly_body = smooth_poly_body
        self.smooth_window_wing = smooth_window_wing
        self.smooth_poly_wing = smooth_poly_wing
        self.zscore_window =  zscore_window
        self.zscore_threashold = zscore_threashold
        self.interp_method = interp_method

    def load_csv(self,mov):
        mov_path = f'{self.exp_path}/{self.exp_name}_mov_{mov}'
        self.angles = pd.read_csv(f'{mov_path}_angles_cm.csv')
        self.fix_phi()
        self.vectors = pd.read_csv(f'{mov_path}_vectors.csv')
        self.ew_to_lab_rotmat = pd.read_csv(f'{mov_path}_ew_to_lab_rotmat.csv', header = None)
        self.exp_file.attrs['ew_to_lab_rotmat'] = self.ew_to_lab_rotmat

    def fix_phi(self):
        phi_col = self.angles.loc[:,self.angles.columns.str.contains('phi')].copy() 
        phi_col[phi_col> 180] = 360 - phi_col
        self.angles.loc[:,self.angles.columns.str.contains('phi')] = phi_col

    def threashold_zscore(self):
        zscore = self.angles.apply( lambda x: abs(self.zscore(x,self.zscore_window)))
        self.angles[zscore > self.zscore_threashold] = None
        self.angles.interpolate(method=self.interp_method)
        return self.angles.dropna().reset_index(drop=True)
    
    def filter_body_wing_and_derive(self):
        angles = self.threashold_zscore()
        self.body = self.filter_and_derive(angles.loc[:, self.angles.columns.str.contains('body')],angles['frames'],angles['time'],self.smooth_window_body,self.smooth_poly_body,derivetive_name = ['','_dot'])
        self.wing = self.filter_and_derive(angles.loc[:, ~self.angles.columns.str.contains('body')],angles['frames'],angles['time'],self.smooth_window_wing,self.smooth_poly_wing,derivetive_name = [''])
         



    def manual_clip_frames(self,mov,ax,ax_twin):
       
        ax.plot(self.wing['frames'],self.wing['phi_rw'],'-*')
        ax_twin.plot(self.body['frames'],self.body['pitch_body'],color = 'red')
        zero_frame = self.wing['frames'][self.wing['time']==0]

        if len(zero_frame)>0: ax_twin.axvline(x = zero_frame.iloc[0], color = 'red',linewidth = 5)


        ax_twin.autoscale()
        ax_twin.relim()
        ax_twin.set_title(f'mov{mov}')
        plt.pause(0.00001)
        x = plt.ginput(2)


        if len(x) > 0: self.save_mov_to_hdf(mov,x)


    def save_mov_to_hdf(self,mov,x):
        mov_group = self.exp_file.create_group(f'/mov{mov}/')
        frames_to_keep = (self.wing['frames']>=int(x[0][0])) & (self.wing['frames']<=int(x[1][0]))

        mov_group['wing_angles'] =  self.wing[frames_to_keep]
        mov_group['body_angles'] = self.body[frames_to_keep]
        mov_group['raw_vectors'] = self.vectors[self.vectors['frames'].isin(self.body['frames'])]
        mov_group['raw_angles'] = self.angles[self.vectors['frames'].isin(self.body['frames'])]



    @staticmethod
    def zscore(x, window):
        # calculate rolling zscore
        # inp: x = data
        #      window = size of window for rolling
        r = x.rolling(window=window)
        m = r.mean().shift(1)
        s = r.std(ddof=0).shift(1)
        z = (x - m) / s
        return z
    
    
    # def filter_and_derive(self,pandas_columns,frames,time,smooth_window,smooth_poly,deriv = 0):
    #    pandas_columns.apply( lambda x: savgol_filter(x,smooth_window,smooth_poly,deriv = deriv))
    #    return pandas_columns.assign(frames=frames, time=time)
    def cut_frames(self,dataframe_dict,frames_to_keep):
        return {name_key:dataframe_dict[name_key][frames_to_keep] for name_key in dataframe_dict.keys()}

    def filter_and_derive(self,pandas_columns,frames,time,smooth_window,smooth_poly,derivetive_name = ['','_dot','_dot_dot']):
       data_frame =  pd.concat([pandas_columns.apply( lambda x: savgol_filter(x,smooth_window,smooth_poly,deriv = deriv)).add_suffix(name) for deriv,name in enumerate(derivetive_name)])
       data_frame = data_frame.assign(frames=frames, time=time)
       return data_frame

