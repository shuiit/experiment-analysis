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

    def cases_plot_exp_mov(self,case,mov_name,x_name,y_name,fig,i,showlegend,color):
        if 'plot_exp' == case: 
            color,name = color,self.pertubation_name
            showlegend = True if i == 0 else False
        if 'plot_mov' == case: color,name = i%len(self.color_map),mov_name
        body_plot = Plotters(self.exp_dict[mov_name].data['body'],self.exp_dict[mov_name].header['body'])
        fig = body_plot.plot_prop_movie(x_name,y_name,self.color_map[color],name,fig = fig,showlegend = showlegend)


    def plot_movies(self,x_name,y_name, case = 'plot_mov', fig = False, mov = False,showlegend = True, color = 0):
        fig = go.Figure() if fig == False else fig
        [self.cases_plot_exp_mov(case,mov_name,x_name,y_name,fig,i,showlegend,color) for i,mov_name in enumerate(list(self.mov_names))]        
        return fig
    
