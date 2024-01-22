import numpy as np
import h5py
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as colormap
import plotly.graph_objects as go
from Experiment import Experiment
from data_analysis.Movie import ManipulatedMovie
from scipy.spatial.transform import Rotation as R
from Plotters import Plotters


def plot_and_save_movie_prop(exp,yprop,wing_body,add_horizontal_line = False):
    
    fig = go.Figure()
    exp.plot_movie_prop('time',yprop, case = 'plot_mov',prop = wing_body, fig = fig, mov = False,add_horizontal_line = add_horizontal_line) 
    fig.show()
    plotly.offline.plot(fig, filename=f'{exp.figures_path}/{yprop}_plot.html',auto_open=False)



if __name__ == '__main__':  


    from Experiment import Experiment
    import matplotlib.pyplot as plt
    from Plotters import Plotters
    import plotly.graph_objects as go
    import plotly
    import numpy as np
    import time

    loadir = 'H:/My Drive/dark 2022/csv_dark' # dir to save the data frame in # H:\My Drive\Ronis Exp\sagiv\data_analysis
    exp_name = ['manipulated_2023_08_06_40ms','manipulated_2023_08_09_60ms']
    save_fig = False
    #            'manipulated_2022_01_31_dark','manipulated_2023_08_07_5ms','manipulated_2023_08_07_10ms']#,'manipulated_2022_02_03_dark','manipulated_2023_08_06_40ms','manipulated_2023_08_09_60ms','manipulated_2023_08_07_5ms','manipulated_2023_08_07_10ms']

    ref_frame = 0
    experiments = [Experiment(loadir,exp) for exp in exp_name]
    min_max_v = [exp.min_max_point_movies('v_xy') for exp in experiments]
    min_max_pitch =[exp.min_max_point_movies('pitch_y_frame') for exp in experiments]
    zero_velocity = [exp.zero_velocity_movies('v_xy') for exp in experiments]
    [exp.pqr_movies() for exp in experiments]
    [exp.mean_mean_props_movies('phi_rw','phi_lw','mean_wing','mean_mean_phi') for exp in experiments]


    fig = go.Figure()

    [exp.plot_prop_movies('v_xy','body',color_map[0],fig,mov = False,case = 'plot_exp',mode = 'markers') for idx,exp in enumerate(experiments)]






# save figures
    if save_fig == True:
        for exp in experiments:

            exp.plot_3d_traj_movies('CM_xy_dot',save_plot = True, mov = False)

            plot_and_save_movie_prop(exp,yprop = 'CM_xy_dot',wing_body = 'body',add_horizontal_line = 0)
            plot_and_save_movie_prop(exp,yprop = 'pitch_y_frame',wing_body = 'body',add_horizontal_line = 0)

            exp.mean_props_movies('phi_rw_min_val','phi_lw_min_val','mean_wing','mean_front_phi')
            exp.mean_props_movies('phi_rw','phi_lw','mean_wing','mean_mean_phi')

            plot_and_save_movie_prop(exp,yprop = 'mean_front_phi',wing_body = 'mean_wing')
            plot_and_save_movie_prop(exp,yprop = 'mean_mean_phi',wing_body = 'mean_wing')
