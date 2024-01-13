import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
import matplotlib.pyplot as plt



pio.renderers.default='browser'

# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 08:25:02 2023

@author: Roni
"""

class Plotters():
    def __init__(self,data,header):
        self.data = data
        self.header = header



    
    def plot_prop_movie(self,x_name,y_name,color,name,
                        marker_size = 7,fig = False,showlegend = True,line_width = 5,mode = 'lines'):
        fig = go.Figure() if fig == False else fig
        x_name_idx = self.header[x_name]
        y_name_idx = self.header[y_name]

        traces = go.Scattergl(
            x=self.data[:,x_name_idx],
            y=self.data[:,y_name_idx],
            legendgroup = name,
            name = name,
            showlegend = showlegend,
            mode=mode,
            line_width= line_width,
            marker=dict(color = f'rgba{str(tuple(np.append(color,0.5)))}',size = marker_size, symbol = 'circle'))
        fig.add_traces(traces)
        fig.update_layout( xaxis_title = x_name, yaxis_title = y_name)
        return fig

        

    