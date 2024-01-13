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
        self.exp_dict = {mov : ManipulatedMovie(self.experiment,mov) for mov in self.mov_names}

    def cases_plot_exp_mov(self,case,mov_name,x_name,y_name,fig,i,color,prop = 'body',showlegend = True):
        if 'plot_exp' == case: color,name,mode = color,self.pertubation_name,'lines'
        if ('plot_exp' == case) | ('plot_mean' == case):showlegend = True if i == 2 else False
        if 'plot_mov' == case: color,name,mode = i%len(self.color_map),mov_name,'lines'
        if 'plot_mean' == case: color, name,mode = i%len(self.color_map),f'stroke_mean','markers'

        body_plot = Plotters(self.exp_dict[mov_name].data[prop],self.exp_dict[mov_name].header[prop])
        fig = body_plot.plot_prop_movie(x_name,y_name,self.color_map[color],name,fig = fig,showlegend = showlegend,mode = mode)


    def plot_movies(self,x_name,y_name, case = 'plot_mov', fig = False, mov = False, color = 0,prop = 'body'):
        mov_names = self.mov_names if mov == False else mov
        fig = go.Figure() if fig == False else fig
        [self.cases_plot_exp_mov(case,mov_name,x_name,y_name,fig,i + 2,color,prop = prop) for i,mov_name in enumerate(list(mov_names))]        
        return fig
    
