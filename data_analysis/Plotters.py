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
    def __init__(self,data,header,pertubation):
        self.data = data
        self.header = header
        self.pertubation = pertubation



    
    def plot_prop_movie(self,x_name,y_name,color,name,
                        marker_size = 7,fig = False,showlegend = True,
                        line_width = 5,mode = 'lines',add_horizontal_line = None):
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
        return fig
    

    def scatter_3d(self,fig,data,hover_data, name,**kwargs):

        fig.add_trace(go.Scatter3d(
        x=data[:,0],
        y=data[:,1],
        z=data[:,2],
        customdata = hover_data.astype(int),**kwargs,
        name = name,
        ))
        fig.update_traces( hovertemplate='<b>time</b>: %{customdata:,.f}<br>') 
        return fig
    
    def plot_3d_traj(self,data,plot_cofnig,mov_name,exp_name,color_prop):

        fig = go.Figure()

        scene=dict(camera=dict(eye=dict(x=1., y=1, z=1.25)), #the default values are 1.25, 1.25, 1.25
        xaxis=dict(nticks=10),
        yaxis=dict(nticks=10),
        zaxis=dict(nticks=10),
        aspectmode='data' #this string can be 'data', 'cube', 'auto', 'manual'
        )

        fig.update_scenes(scene)


        vector_kwargs = {'marker':dict(size=5,opacity=0.8,color = 'red'),'line':dict(color='black' ,width=8)}
        cm_kwargs = {'marker':dict(size=8,opacity=0.8,color =data[color_prop][::plot_cofnig['traj_samples']],colorbar=dict(thickness=10,title = color_prop),colorscale='jet',colorbar_x = -0.1),'mode' : 'markers'}
        pertubations = {'marker':dict(size=13,opacity=0,color = ['pink','green','red']),'mode' : 'markers'}

        start_pert_endpert = data['start_pert_endpert'] 

        fig = self.scatter_3d(fig,data['cm'][::plot_cofnig['traj_samples'],:],
                        data['time'][::plot_cofnig['traj_samples']],name = 'CM',**cm_kwargs)

        fig = self.scatter_3d(fig,data['x_vector'],data['time'][::plot_cofnig['fly_samples']],name = 'X vector',**vector_kwargs)
        fig = self.scatter_3d(fig,data['y_vector'],data['time'][::plot_cofnig['fly_samples']],name = 'Y vector',**vector_kwargs)
        fig = self.scatter_3d(fig,data['cm'][start_pert_endpert,:],data['time'][start_pert_endpert],name =  'pertubations',**pertubations)
        fig.update_layout( title = f'{mov_name} {exp_name}',coloraxis_colorbar_y=0.5,coloraxis_colorscale='jet' )
        return fig
   









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

            

        