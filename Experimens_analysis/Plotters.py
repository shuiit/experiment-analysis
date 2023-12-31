import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px
import matplotlib.cm
import numpy as np
import matplotlib.pyplot as plt
from Experiment import Experiment
import pandas as pd
from plotly.subplots import make_subplots


pio.renderers.default='browser'

# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 08:25:02 2023

@author: Roni
"""

class Plotter():
    def __init__(self,experiment):
        
        self.data = experiment.data
        self.body_vectors = experiment.body_vectors
        self.mean_data = experiment.mean_data
        self.exp_name = experiment.exp_name
        self.pert_time = experiment.pert_time
        self.interest_points = experiment.interest_points
        self.mov_order = experiment.mov_order
    
        
        
    def disconnect_line_add_none(self,array1,array2):
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
        
        
        
        
    def plot3d_traj(self,movie_name,fly_samples = 150,traj_samples = 20,size_x = 1/2000,size_y = 1/3000):
        """
        Plots a 3D trajectory visualization of a movie.
    
        Args:
            movie_name (str): The name of the movie.
            fly_samples (int, optional): Sampling rate for selecting center of mass (CM) of the fly. Default is 150.
            traj_samples (int, optional): Sampling rate for selecting CM of the trajectory. Default is 20.
            size_x (float, optional): Scaling factor for X body vector. Default is 1/2000.
            size_y (float, optional): Scaling factor for Y body vector. Default is 1/3000.
    
        Returns:
            None
    
        Example:
            plot3d_traj('mov_name', fly_samples=150, traj_samples=20, size_x=1/2000, size_y=1/3000)
        """

        fig = go.Figure()
  
        mov = Experiment.get_mov(self.data,[movie_name])
        center_mass = mov.iloc[::fly_samples,:] # CM of fly
        center_mass_traj = mov.iloc[::traj_samples,:] # CM of trajectory
        
        
        body_vectors = self.body_vectors[movie_name] # get body vectors
        frames_body_vectors = np.where(np.in1d(body_vectors['frames'],center_mass['frames']))
        body_x_vector = (body_vectors['X'][frames_body_vectors])*size_x + center_mass[['X','Y','Z']] # define the size of Xbody vecotr
      
        # generate two points for each "fly" and insert nan to disconnect line
        body_x_vector = self.disconnect_line_add_none(center_mass[['X','Y','Z']],body_x_vector)
        
        delta_y_on_x = (body_vectors['X'][frames_body_vectors])*size_y + center_mass[['X','Y','Z']] # define the size of Ybody vecotr
        body_y_vector = body_vectors['Y'][frames_body_vectors]*size_y
        
        body_y_vector = self.disconnect_line_add_none(-body_y_vector + delta_y_on_x,body_y_vector+delta_y_on_x)

        zero_time = np.where(mov['time [ms]'] == 0)
        
        # plot------------------------------------------------------------
        
        scene=dict(camera=dict(eye=dict(x=1.15, y=1.15, z=0.8)), #the default values are 1.25, 1.25, 1.25
        xaxis=dict(nticks=10),
        yaxis=dict(nticks=20),
        zaxis=dict(nticks=10),
        

        aspectmode='data' #this string can be 'data', 'cube', 'auto', 'manual'
        )
        fig.update_scenes(scene)
        fig.add_trace(go.Scatter3d(
        x=center_mass_traj['X'],
        y=center_mass_traj['Y'],
        z=center_mass_traj['Z'],
        customdata = center_mass_traj['time [ms]'].astype(int),
        hovertemplate='<b>time</b>: %{customdata:,.f}<br>',
        mode='markers',
        marker=dict(size=8,opacity=0.8),name = 'CM location'))
        
        
        fig.add_trace(go.Scatter3d(
        x=body_x_vector[:,0],
        y=body_x_vector[:,1],
        z=body_x_vector[:,2],
        mode='lines+markers',
        hoverinfo='skip',
        line=dict(color='black' ,width=10),marker = dict(color = 'red',size = 5),name = 'Xbody'))
    
        fig.add_trace(go.Scatter3d(
        x=body_y_vector[:,0],
        y=body_y_vector[:,1],
        z=body_y_vector[:,2],
        hoverinfo='skip',
        mode='lines+markers',
        line=dict(color='grey' ,width=10),marker = dict(color = 'red',size = 5),name = 'Ybody'))
    
        fig.add_trace(go.Scatter3d(
        x=[center_mass_traj['X'].iloc[0]],
        y=[center_mass_traj['Y'].iloc[0]],
        z=[center_mass_traj['Z'].iloc[0]],
        customdata =[ center_mass_traj['time [ms]'].iloc[0].astype(int)],
        hovertemplate='<b>time</b>: %{customdata:,.f}<br>',
        mode='markers',
        marker=dict(color = 'pink',size =13),name = 'begining of video'))
        
        fig.add_trace(go.Scatter3d(
        x=mov['X'].iloc[zero_time],
        y=mov['Y'].iloc[zero_time],
        z=mov['Z'].iloc[zero_time],
        customdata = mov['time [ms]'].iloc[zero_time].astype(int),
        hovertemplate='<b>time</b>: %{customdata:,.f}<br>',
        mode='markers',
        marker=dict(color = 'green',size =13),name = 'dark start time'))
        
        if self.pert_time[1] is not None:
            pert_end = np.where(mov['time [ms]'] == self.pert_time[1])
            fig.add_trace(go.Scatter3d(
            x=mov['X'].iloc[pert_end],
            y=mov['Y'].iloc[pert_end],
            z=mov['Z'].iloc[pert_end],
            customdata = mov['time [ms]'].iloc[pert_end].astype(int),
            hovertemplate='<b>time</b>: %{customdata:,.f}<br>',
            mode='markers',
            marker=dict(color = 'red',size =13),name = 'dark end time'))

        pio.show(fig)

        
    def plot_props_per_movie(self,propy,propx = 'time [ms]',movie_name = 'all',alpha= 0.5,size = 5, mean_data = False,return_fig = None,scatter_kwargs=None, layout_kwargs=None):
        """
        Plot properties per movie. Each movie is a different color

        Parameters
        ----------
        propy (str): y axis property to plot
        propx (str,optional): x axis property to plot. The default is 'time [ms]'.
        movie_name (str,optional): movies to plot, if 'all' - plot all The default is 'all'.
        marker_mode (str,optional): marker\marker + line\... The default is 'markers'.
        alpha (float,optional): opacity. The default is 0.5.
        size (int,optional): size of markers. The default is 5.
        mean_data (bool, optional): Plot mean_data. The default is False - dont plot.
        return_fig (fig, optional): return the object fig, used to add traces. The default is None.
        layout_kwargs (dictionary, optional): pass variables for 'update_layout'
        scatter_kwargs (dictionary, optional): pass variables for 'Scattergl'

        Returns
        -------
        fig : TYPE
            DESCRIPTION.

        """
        movie_data = self.data if not mean_data else self.mean_data
        scatter_kwargs = scatter_kwargs or {}
        layout_kwargs = layout_kwargs or {}
            
        if movie_name != 'all': movie_data = Experiment.get_mov(movie_data,movie_name)
            
        colorcd = matplotlib.cm.datad["tab10"]['listed']
        fig = go.Figure() 
        grouped_mov = movie_data.groupby('mov')
        # plot movs, different color for each movie
     
            
        traces = [(go.Scattergl(
            x=movie_data[propx],
            y=movie_data[propy],
            name = movie_name,
            hovertemplate = '<b>' + 'x' + '</b>' +   ' = %{x:,.2f} ' + '<b>'+ 'y' +'</b>' +  ' = %{y:,.2f} ' + '<b>'+ movie_name +'</b>' + '<extra></extra>',
            legendgroup = movie_name,
            marker=dict(color = f'rgba({", ".join(str(c * 255) for c in colorcd[i%len(colorcd)])}, {alpha})',size = size),**scatter_kwargs))
            for i,(movie_name, movie_data) in enumerate(grouped_mov)]
        fig.add_traces(traces)
        fig.update_layout( xaxis_title = propx, yaxis_title = propy,**layout_kwargs) 
    
        return fig if return_fig is not None else fig.show() 
 
       
        
    def plot_props_per_experiment(self,propy,propx = 'time [ms]',movie_name = 'all',alpha= 0.5,size = 5, return_fig = None, mean_data = False, fig = None,**kwargs):
        """
        Plot properties per experiment.

        Parameters
        ----------
        propy (str): y axis property to plot
        propx (str,optional): x axis property to plot. The default is 'time [ms]'.
        movie_name (str,optional): movies to plot, if 'all' - plot all The default is 'all'.
        alpha (float,optional): opacity. The default is 0.5.
        size (int,optional): size of markers. The default is 5.
        mean_data (bool, optional): Plot mean_data. The default is False - dont plot.
        return_fig (fig, optional): return the object fig, used to add traces. The default is None.

        
        color (str, optional): Color name, example: 'black' The default is None.
        fig (fig, optional): use the figure the user provides, used for add trace. The default is None.

        Returns
        -------
        fig : TYPE
            DESCRIPTION.

        """
        
        data = self.data if not mean_data else self.mean_data
        data = data.groupby('mov').apply(lambda x: pd.concat([pd.DataFrame([np.full(len(x.columns), None)], columns=x.columns), x], ignore_index=True)).reset_index(drop = True)

            
        if movie_name != 'all': data = Experiment.get_mov(data,movie_name)

        fig = go.Figure() if fig is None else fig
        
    
        grp = data
        fig.add_trace(go.Scattergl(
        x=grp[propx],
        y=grp[propy],
        name = self.exp_name,
        hovertemplate = '<b>' + 'x' + '</b>' +   ' = %{x:,.2f} ' + '<b>'+ 'y' +'</b>' +  ' = %{y:,.2f} ' + '<b>'+ self.exp_name +'</b>' + '<extra></extra>',
        showlegend = True,
        legendgroup = self.exp_name,
        **kwargs))
        fig.update_layout( xaxis_title = propx,
         yaxis_title = propy)
        return fig if return_fig is not None else fig.show()  
        
        
       
        
    def plot_min_on_movie(self,prop,idx_prop,plot_min = False,**kwargs):
        """
        Mark the minimum value or the values that cross zero on each movie with a black marker

        Parameters
        ----------
        plot_min (bool, optional). If True plot minimum, if False plot zero crossing The default is False.

        Returns
        -------
        None.

        """

        fig = go.Figure()
        
        # plot movs, different color for each movie
        
        fig = self.plot_props_per_movie(prop,propx = 'time [ms]',movie_name = 'all',alpha= 0.5,return_fig = True)
        idx_to_plot = self.interest_points[idx_prop] 
        # plot minimum on the graph
        
        traces = [go.Scattergl(
        x=self.data['time [ms]'][[indices]],
        y=self.data[prop][[indices]],
        mode='markers',
        hoverinfo='skip',
        legendgroup = movie_name,
        marker = dict(color = 'black',size = 8),
        showlegend=False,
        ) for movie_name,indices  in idx_to_plot.items() if indices is not None ]
       
        fig.add_traces(traces)    
        fig.update_layout( **kwargs)
        fig.show()
        
    def idx_histogram(self,prop,x_prop = 'time [ms]', fig = None, return_fig = True,nbinsx = 10,trace_kwargs = None,layout_kwargs = None):
        """
        plot histogram of time of minimum \ zero crossing

        Parameters
        ----------
        plot_min (bool, optional). If True plot minimum, if False plot zero crossing The default is False.
        fig (fig, optional): use the figure the user provides, used for add trace. The default is None.
        return_fig (fig, optional): return the object fig, used to add traces. The default is None.
        nbinsx (int, optional): number of bins. The default is 10.
        trace_kwargs (dictionary, optional): **kwargs to add to trace, such as marker = {opacity: 0.4} 
        layout_kwargs (dictionary, optional): **kwargs to add to layout, such as marker = {title: 'foo'} 

        Returns
        -------
        fig : TYPE
            DESCRIPTION.

        """
        trace_kwargs = trace_kwargs or {}
        layout_kwargs = layout_kwargs or {}
        fig = go.Figure() if fig is None else fig 
            
        idx_to_plot = self.interest_points[prop].dropna()
        cross_zero_idx = np.hstack(idx_to_plot)
        time = self.data[x_prop].iloc[cross_zero_idx]
        
        fig.add_trace(go.Histogram(x=time, name=self.exp_name + ' ' + prop, nbinsx=nbinsx))
        fig.update_layout( **layout_kwargs)
        fig.update_traces(**trace_kwargs)
        return fig if return_fig is not None else fig.show()
    
    
    def histogram(self,prop, fig = None, return_fig = True,nbinsx = 10,**kwargs):
        """
        plot histogram of time of minimum \ zero crossing

        Parameters
        ----------
        plot_min (bool, optional). If True plot minimum, if False plot zero crossing The default is False.
        fig (fig, optional): use the figure the user provides, used for add trace. The default is None.
        return_fig (fig, optional): return the object fig, used to add traces. The default is None.
        nbinsx (int, optional): number of bins. The default is 10.

        Returns
        -------
        fig : TYPE
            DESCRIPTION.

        """
        fig = go.Figure() if fig is None else fig 
            
        data_to_plot = self.interest_points[prop].dropna()
        
        fig.add_trace(go.Histogram(x=np.hstack(data_to_plot), name=self.exp_name + ' ' + prop, nbinsx=nbinsx))
        fig.update_layout( **kwargs)
        return fig if return_fig is not None else fig.show()
    
    def plot_multi_y(self,mov_name,props):
        """
        plot a multi Y axis graph of the same movie 

        Parameters
        ----------
        mov_name (str): name of movie
        props (str): name of properties to plot

        Returns
        -------
        None.

        """
        
        fig = go.Figure()
        data = Experiment.get_mov(self.data,[mov_name])
        colorcd = matplotlib.cm.datad["tab10"]['listed']
        yax_dict = {}
        side = ['left','right','left']
        position = [0,0.05,1]
        
        for i,prop in enumerate(props):
            color = f'rgb({", ".join(str(c * 255) for c in colorcd[i%len(colorcd)])})'
            fig.add_trace(go.Scattergl(
                        x=data['time [ms]'],
                        y=data[prop],
                        name = prop,
                        marker = dict(color = f'rgb({", ".join(str(c * 255) for c in colorcd[i%len(colorcd)])})',size = 8),
                        hovertemplate='<b>' + 'x' + '</b>' + ' = %{x:,.2f} ' + '<b>' + 'y' + '</b>' + ' = %{y:,.2f} ' + '<b>' + prop + '</b>' + '<extra></extra>',
                        yaxis=f"y{i + 1}"
                    ))  
            yax_dict[f"yaxis{i + 1}"] = (dict(title=prop,titlefont=dict(color=color),tickfont=dict(color=color),side = side[i],position=position[i]))
            

        fig.update_layout(title = f"<b> {mov_name} - {self.exp_name} <b>",xaxis = dict(domain = [0.06,1]),**yax_dict)
        fig.show()
                        
        
        
        
    
    
    
    