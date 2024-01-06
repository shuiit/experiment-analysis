import pandas as pd
import matplotlib.pyplot as plt
import h5py
import os
import sys
from scipy.signal import savgol_filter
import numpy as np

class PreProcess:

    def __init__(self,exp_dir,exp_name,pert_for_attr, smoothing_config,hdf5_file_name):
        self.exp_dir = exp_dir
        self.exp_name = exp_name
        self.exp_path = f'{exp_dir}/{exp_name}'
        self.pert_for_attr = pert_for_attr
        self.exp_file  = h5py.File(f'{exp_dir}/{hdf5_file_name}.hdf5', "a")
        self.exp_file.attrs['dark_pert'] = self.pert_for_attr
        self.smoothing_config = smoothing_config
        



    def load_csv(self,mov):
        """ Load csv files: _angles_cm.csv and _vectors.csv. save the easywand to lab rotation matrix
        to the experiment attribution. 

        Args:
            mov (int): numberof movie
        """
        mov_path = f'{self.exp_path}/{self.exp_name}_mov_{mov}'
        self.angles = pd.read_csv(f'{mov_path}_angles_cm.csv')
        self.vectors = pd.read_csv(f'{mov_path}_vectors.csv')
        self.fix_phi()

        self.ew_to_lab_rotmat = pd.read_csv(f'{mov_path}_ew_to_lab_rotmat.csv', header = None)
        self.exp_file.attrs['ew_to_lab_rotmat'] = self.ew_to_lab_rotmat
        self.dt = np.diff(self.angles['time'])[0]/1000
    def fix_phi(self):
        """ makes sure that the phi is not the 360 - phi degree. 
        """
        phi_col = self.angles.loc[:,self.angles.columns.str.contains('phi')].copy() 
        phi_col[phi_col > 180] = 360 - phi_col
        self.angles.loc[:,self.angles.columns.str.contains('phi')] = phi_col

    def threashold_zscore(self):
        """ Calculate the Z score of the data, threashold it and interpulate the angles. 

        Returns:
            dataFrame: the dataframe after interpulating, zscore threasholding, and droping Nan values
        """
        zscore = self.angles.apply( lambda x: abs(self.zscore(x,self.smoothing_config['zscore_window'])))
        self.angles[zscore > self.smoothing_config['zscore_threashold']] = None
        self.angles.interpolate(method=self.smoothing_config['interp_method'])
        return self.angles.dropna().reset_index(drop=True)
    
    def filter_body_wing_and_derive(self):
        """ Threashold using Z score and use savitzky golay filter to smooth and calculate derivatives
        """
        angles = self.threashold_zscore()
        body_columns = self.angles.columns.str.contains('body')
        cm_columns = self.angles.columns.str.contains('CM')
        body_angles = self.filter_and_derive(angles.loc[:, (body_columns) & (~cm_columns)],
                                           angles['frames'],angles['time'],'smooth_window_body',
                                           'smooth_poly_body',derivetive_name = ['','_dot','_dot_dot'])
        
        body_cm = self.filter_and_derive(angles.loc[:, (body_columns) & (cm_columns)],
                                           angles['frames'],angles['time'],'smooth_window_body',
                                           'smooth_poly_body',derivetive_name = [''],assign_time_frame = False)
        self.body = pd.concat([body_angles,body_cm],axis = 1)
        self.wing = self.filter_and_derive(angles.loc[:, ~self.angles.columns.str.contains('body')],
                                           angles['frames'],angles['time'],'smooth_window_wing',
                                           'smooth_poly_wing',derivetive_name = [''])
         

    def manual_clip_frames(self,mov,ax,ax_twin):
        """Manually mark the initial and ending frames of data to save

        Args:
            mov (int): number of movie
            ax (axes): axes object of the plot
            ax_twin (axes): axes object of the second plot (ploted on the same figure)
        """
       
        ax.plot(self.wing['frames'],self.wing['phi_rw'],'-*')
        ax_twin.plot(self.body['frames'],self.body['pitch_body'],color = 'red')
        zero_frame = self.wing['frames'][self.wing['time']==0]

        if len(zero_frame)>0: ax_twin.axvline(x = zero_frame.iloc[0], color = 'red',linewidth = 5)
        ax_twin.autoscale(),ax_twin.relim(),ax_twin.set_title(f'mov{mov}')
        plt.pause(0.00001)
        x = plt.ginput(2)
        if len(x) > 0: self.save_mov_to_hdf(mov,x)


    def save_mov_to_hdf(self,mov,x):
        """Save the movie to HDF5 file, save the data to a subgroup with the name of the movie,
        the dataset will get the name of datasets_name. save the heading of all datasets as an attributes

        Args:
            mov (int):  number of movie
            x (list of int): initial and ending frames
        """
        mov_group = self.exp_file.create_group(f'/mov{mov}/')
        frames_to_keep = (self.wing['frames']>=int(x[0][0])) & (self.wing['frames']<=int(x[1][0]))
        
        data_list = [self.wing[frames_to_keep],self.body[frames_to_keep],
                     self.vectors[self.vectors['frames'].isin(self.wing[frames_to_keep]['frames'])],
                     self.angles[self.angles['frames'].isin(self.wing[frames_to_keep]['frames'])]]
        datasets_name = ['wing_angles','body_angles','vectors_raw','angles_raw']
        [self.create_datasets(mov_group,sub_group_name,data) for sub_group_name,data in zip(datasets_name,data_list)]


    def filter_and_derive(self,pandas_columns,frames,time,smooth_window,smooth_poly,derivetive_name = ['','_dot','_dot_dot'],assign_time_frame = True):
        """ run savitky golay filter and save the dataframe with time and frame columns

        Args:
            pandas_columns (DataFrame): a dataframe to run the filter on (wing/body/angles...)
            frames (DataFrame): frames of the dataFrame 
            time (DataFrame): time of frames (ms)
            smooth_window (int): size of window to smooth (different if body (211) or wing (7))
            smooth_poly (int): degree of polinom (different if body (4) or wing (2))
            derivetive_name (list, optional): name of suffix (the length of the name list defines how many times to derive
             (example ['','_dot','_dot_dot'] will generate: smoothing the data, first and second derivative ). Defaults to ['','_dot','_dot_dot'].

        Returns:
            data_frame: smoothed dataframe with time and frame columns (not smoothed)
        """
        data_frame =  pd.concat([pandas_columns.apply( lambda x: savgol_filter(x/self.dt**deriv,self.smoothing_config[smooth_window],
                                                                              self.smoothing_config[smooth_poly],deriv = deriv)).add_suffix(name) for deriv,name in enumerate(derivetive_name)],axis = 1)
        if assign_time_frame == True: data_frame = data_frame.assign(frames=frames, time=time)
        return data_frame



    @staticmethod
    def zscore(x, window):
        """Calculate Z score (represents how much athe data is "standard deviation")

        Args:
            x (np.array/dataframe): the data to calculate z score for
            window (int): size of window to calculate the z score for

        Returns:
            z: z score
        """
        r = x.rolling(window=window)
        m = r.mean().shift(1)
        s = r.std(ddof=0).shift(1)
        z = (x - m) / s
        return z

    @staticmethod
    def create_datasets(group,sub_group_name,data):
        """create fields of datasets and save it in the HDF5 file, save the heading of the data as attribution

        Args:
            group (HDF5 object): a group of the hdf5 file
            sub_group_name (string): name of the group
            data (dataframe): dataframe that contains all data
        """
        group[sub_group_name] =  data
        group[sub_group_name].attrs['header'] = list(data.columns)


    