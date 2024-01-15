import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px
import matplotlib.cm
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from plotly.subplots import make_subplots
import h5py
from ManipulatedMovie import ManipulatedMovie
import matplotlib.cm as colormap
from Plotters import Plotters

pio.renderers.default='browser'

class Experiment():
    def __init__(self,loadir,exp_name, movie_name_list = False):    
        self.experiment = h5py.File(f'{loadir}/{exp_name}.hdf5', "r")
        self.experiment_name = exp_name
        self.pertubation_name = self.experiment_name.split('_')[-1]
        self.color_map = colormap.datad["tab10"]['listed']
        self.mov_names = list(self.experiment.keys()) if movie_name_list == False else movie_name_list
        self.pertubation = int(self.pertubation_name.split('ms')[0]) if self.pertubation_name.isalpha() == False else False

        self.exp_dict = {mov : ManipulatedMovie(self.experiment,mov,pertubation = self.pertubation) for mov in self.mov_names}
        self.body_header = self.exp_dict[self.mov_names[0]].header['body']
        self.wing_header = self.exp_dict[self.mov_names[0]].header['wing']


    def cases_plot_exp_mov(self,case,mov_name,x_name,y_name,fig,i,color,prop = 'body',showlegend = True,add_horizontal_line = True):
        if 'plot_exp' == case: color,name,mode = color,self.pertubation_name,'lines'
        if ('plot_exp' == case) | ('plot_mean' == case):showlegend = True if i == 2 else False
        if 'plot_mov' == case: color,name,mode = i%len(self.color_map),mov_name,'lines'
        if 'plot_mean' == case: color, name,mode = i%len(self.color_map),f'stroke_mean','markers'

        body_plot = Plotters(self.exp_dict[mov_name].data[prop],self.exp_dict[mov_name].header[prop])
        fig = body_plot.plot_prop_movie(x_name,y_name,self.color_map[color],name,fig = fig,showlegend = showlegend,mode = mode,add_horizontal_line = add_horizontal_line)


    def plot_movies(self,x_name,y_name, case = 'plot_mov',
                    fig = False, mov = False, color = 0,prop = 'body',add_horizontal_line = True):
        mov_names = self.mov_names if mov == False else mov
        fig = go.Figure() if fig == False else fig
        [self.cases_plot_exp_mov(case,mov_name,x_name,y_name,fig,i + 2,color,prop = prop,add_horizontal_line = add_horizontal_line) for i,mov_name in enumerate(list(mov_names))]        
        return fig
    
    
    
    def min_max_point(self,prop):
        return np.hstack([self.exp_dict[mov_name].min_max_point(prop) for mov_name in self.mov_names])
    

    def zero_velocity(self,prop):
        zero_v_list = [self.exp_dict[mov_name].zero_velocity(prop) for mov_name in self.mov_names]
        return [np.min(velocity) for velocity in zero_v_list if np.min(velocity) != None]
    

    def get_mov(self,mov_name):
        return self.exp_dict[mov_name]
    
    
    

    
    