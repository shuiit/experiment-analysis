import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px
import matplotlib.cm
import numpy as np
import os
import plotly
import matplotlib.pyplot as plt
import pandas as pd
from plotly.subplots import make_subplots
import h5py
from Movie import Movie
import matplotlib.cm as colormap
from Plotters import Plotters

pio.renderers.default='browser'

class Experiment():
    def __init__(self,loadir,exp_name, movie_name_list = False, movie_length = 3000):    
        self.experiment = h5py.File(f'{loadir}/{exp_name}.hdf5', "r")
        self.experiment_name = exp_name
        self.pertubation_name = self.experiment_name.split('_')[-1]
        self.color_map = colormap.datad["tab10"]['listed']
        self.mov_names = list(self.experiment.keys()) if movie_name_list == False else movie_name_list
        self.pertubation = int(self.pertubation_name.split('ms')[0]) if self.pertubation_name.isnumeric() == True else False
        self.loadir =loadir
        time_idx = np.where(self.experiment[self.mov_names[0]]['body'].attrs['header'] == 'time')[0][0]
        self.exp_dict = {mov : Movie(self.experiment,mov,pertubation = self.pertubation) for mov in self.mov_names if self.del_initial_tim_and_length(self.experiment[mov],movie_length,time_idx) !=True}
        self.mov_names = list(self.exp_dict.keys())
        self.body_header = self.exp_dict[self.mov_names[0]].header['body']
        self.wing_header = self.exp_dict[self.mov_names[0]].header['wing']

        self.figures_path = f'{self.loadir}/figures/{self.experiment_name.split("manipulated_")[1]}'
        if not os.path.exists(self.figures_path): os.makedirs(self.figures_path)
        self.experiment.close()

    def pqr_movies(self, mov = False):
        mov_names = self.mov_names if mov == False else mov
        [self.get_mov(mov_name).calculate_pqr_update_data_header() for mov_name in mov_names]


    def plot_3d_traj_movies(self,color_prop,save_plot = False, mov = False):
        mov_names = self.mov_names if mov == False else mov
        for mov in  list(mov_names):
            fig = self.exp_dict[mov].plot_3d_traj_movie(color_prop)
            plotly.offline.plot(fig, filename=f'{self.figures_path}/traj_3d_{mov}.html',auto_open=False) if save_plot == True else fig.show()

    def mean_mean_props_movies(self,prop1,prop2,wing_body,header_name, mov = False):
        mov_names = self.mov_names if mov == False else mov
        [self.get_mov(mov_name).mean_props(prop1,prop2,wing_body,header_name) for mov_name in mov_names]

    
    def min_max_point_movies(self,prop,**kwargs):
        return np.hstack([self.exp_dict[mov_name].min_max_point(prop,**kwargs) for mov_name in self.mov_names])
    

    def zero_velocity_movies(self,prop):
        zero_v_list = [self.exp_dict[mov_name].zero_velocity(prop) for mov_name in self.mov_names]
        return [velocity for velocity in zero_v_list if np.min(velocity) != None]
    
    def add_mean_prop_movies(self,prop_name,wing_body_prop,wing_body_mean_save,mov = False,phi_idx_to_mean = 'phi_rw_min_idx'):
        mov_names = self.mov_names if mov == False else mov
        [self.get_mov(mov_name).mean_prop_stroke(prop_name,wing_body_prop,wing_body_mean_save,phi_idx_to_mean = phi_idx_to_mean) for mov_name in mov_names]

    def plot_prop_movies(self,prop,wing_body,fig,mov = False,case = 'plot_mov',prop_x = 'time',
                         add_horizontal_line = 0,color = False,legend = False,**kwargs):
        mov_names = self.mov_names if mov == False else mov
        legend = mov_names if legend == False else legend

        color = [self.color_map[idx%len(self.color_map)] for idx,mov_name in enumerate(mov_names)] if color == False else color 
        if 'plot_exp' == case: [self.get_mov(mov_name).plot_prop(prop,wing_body,color,legend,fig,showlegend = idx == 0,prop_x = prop_x,**kwargs) for idx,mov_name in enumerate(mov_names)]
        if 'plot_mov' == case: [self.get_mov(mov_name).plot_prop(prop,wing_body,color[idx],legend[idx],fig,prop_x = prop_x,**kwargs) for idx,mov_name in enumerate(mov_names)]

        if add_horizontal_line != None: fig.add_hline(y=add_horizontal_line, line_width=3, line_color="black")
        fig.add_vline(x=0, line_width=3, line_color="lime")
        if self.pertubation != False: fig.add_vline(x=self.pertubation, line_width=3, line_color="red")
        fig.update_layout( xaxis_title = prop_x, yaxis_title = prop)     

    def smooth_prop_movies(self,prop,derives,wing_body):
        [self.get_mov(mov_name).smooth_and_derive(prop,derives,wing_body) for  mov_name in self.mov_names]


    def project_prop_movies(self,prop_to_project,**kwargs):
        [self.get_mov(mov_name).project_prop(prop_to_project,'body',**kwargs) for  mov_name in self.mov_names]

    def get_mov(self,mov_name):
        return self.exp_dict[mov_name]
    
    def mean_time_series_prop(self,prop,wing_body,window = (73*7)//2,dt = 1/16000):
        
        t_list = np.arange(-40,400,dt*1000)
        cmlist = np.full((len(self.mov_names),len(t_list)),np.nan)
        for idx_mov,mov_name in enumerate(self.mov_names):
            mov = self.get_mov(mov_name)
            time = mov.get_prop('time',wing_body)[window:,0]

            cm_dot = mov.get_prop(prop,wing_body)[window:,0]
            # cm_dot = cm_dot - np.mean(cm_dot[window:mov.ref_frame + 1])
            time_inter,tlist_idx,time_idx = np.intersect1d(t_list,time,return_indices = True )
            cmlist[idx_mov,tlist_idx] = cm_dot[time_idx]
            mean_prop = np.nanmean(cmlist,axis = 0)
        return mean_prop,t_list
    

    
    @staticmethod
    def del_initial_tim_and_length(mov,movie_length,time_idx):
        if (mov['body'].shape[0] < movie_length) | (mov['body'][0,time_idx] > 0):
            return True

    
    
    

    
    