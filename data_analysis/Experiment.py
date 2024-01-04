import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px
import matplotlib.cm
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from plotly.subplots import make_subplots
import h5py

pio.renderers.default='browser'

class Experiment():
    def __init__(self,loadir,exp_name):    
        self.experiment = h5py.File(f'{loadir}/{exp_name}.hdf5', "r")
        self.experiment_name = exp_name

    # def get_data(self,mov_name,group_name,prop_xy):
    #     data = self.experiment[mov_name][group_name]
    #     propx_idx = [int(np.where(data.attrs['header'] == prop)[0]) for prop in prop_xy]
    #     return data[:,propx_idx[0]],data[:,propx_idx[1]]
    
    # def plot_prop_movie(self,fig,x,y,
    #                     color,name, size = 5,showlegend = True):
    #     traces = go.Scattergl(
    #         x=x,
    #         y=y,
    #         legendgroup = name,
    #         name = name,
    #         showlegend = showlegend,
    #         marker=dict(color = f'rgba{str(tuple(np.append(color,0.5)))}',size = size))
    #     fig.add_traces(traces)

    # def plot_mov(self,fig,mov_name,group_name,prop_xy,color,legend_name,**kwargs):
    #     x,y = self.get_data(mov_name,group_name,prop_xy)
    #     self.plot_prop_movie(fig,x,y,color,legend_name, **kwargs)

        