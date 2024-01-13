import numpy as np
import h5py
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as colormap
import plotly.graph_objects as go
from Experiment import Experiment
from ManipulatedMovie import ManipulatedMovie
from scipy.spatial.transform import Rotation as R
from Plotters import Plotters


if __name__ == '__main__':  

    loadir = 'H:/My Drive/dark 2022/csv_dark/' # dir to save the data frame in # H:\My Drive\Ronis Exp\sagiv\data_analysis
    exp_name = ['manipulated_2023_08_10_100ms','manipulated_2022_02_03_dark']#['2022_02_03','2022_03_10','2022_05_19','2022_03_03'] # 22_11_28
    ref_frame = 0

    experiment = Experiment(loadir,exp_name[0])
    experiment2 = Experiment(loadir,exp_name[1])


    fig = experiment.plot_movies('time','CM_xy_dot', case = 'plot_exp', fig = False, mov = False, color = 2)
    fig = experiment2.plot_movies('time','CM_xy_dot', case = 'plot_exp', fig = fig, mov = False, color = 3)

    fig = experiment.plot_movies('time','pitch_y_frame', case = 'plot_mean', fig = fig, mov = False, prop = 'mean_body')

    fig.show()

#%%
    fig = experiment.plot_movies('time','pitch_y_frame', case = 'plot_mean', fig = fig, mov = False, prop = 'mean_body')

# %%
