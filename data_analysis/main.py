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


def plot_and_save_movie_prop(exp,yprop,wing_body,add_horizontal_line = False):
    
    fig = go.Figure()
    exp.plot_movie_prop('time',yprop, case = 'plot_mov',prop = wing_body, fig = fig, mov = False,add_horizontal_line = add_horizontal_line) 
    fig.show()
    plotly.offline.plot(fig, filename=f'{exp.figures_path}/{yprop}_plot.html',auto_open=False)



if __name__ == '__main__':  


    loadir = 'H:/My Drive/dark 2022/csv_dark' # dir to save the data frame in # H:\My Drive\Ronis Exp\sagiv\data_analysis
    exp_name = ['manipulated_2022_02_03_dark','manipulated_2023_08_06_40ms','manipulated_2023_08_09_60ms','manipulated_2023_08_07_5ms','manipulated_2023_08_07_10ms']

    # ['manipulated_2023_08_10_100ms','manipulated_2022_02_03_dark','manipulated_2023_08_06_40ms',
                    # 'manipulated_2023_08_09_60ms','manipulated_2023_08_07_5ms','manipulated_2023_08_07_10ms'
                    #['2022_02_03','2022_03_10','2022_05_19','2022_03_03'] # 22_11_28

    # exp_name = ['manipulated_2022_05_19_40ms']
    ref_frame = 0

    experiments = [Experiment(loadir,exp) for exp in exp_name]

    min_max_v = [exp.min_max_point('CM_xy_dot') for exp in experiments]
    min_max_pitch =[exp.min_max_point('pitch_y_frame') for exp in experiments]
    zero_velocity = [exp.zero_velocity('CM_xy_dot') for exp in experiments]



    def mean_props(mov,prop1,prop2,wing_body,header_name):
        mean_prop = (mov.get_prop(prop1,wing_body)  + mov.get_prop(prop2,wing_body) )/2
        mov.data[wing_body] = np.hstack((mov.data[wing_body], mean_prop[np.newaxis,:].T))
        mov.add_to_header([header_name],wing_body)

    [mean_props(experiments[0].get_mov(mov_name),'phi_rw_min_val','phi_lw_min_val','mean_wing','mean_front_phi') for mov_name in experiments[0].mov_names]




# save figures
    for exp in experiments:

        exp.plot_3d_traj_movies('CM_xy_dot',save_plot = True, mov = False)

        plot_and_save_movie_prop(exp,yprop = 'CM_xy_dot',wing_body = 'body',add_horizontal_line = 0)
        plot_and_save_movie_prop(exp,yprop = 'pitch_y_frame',wing_body = 'body',add_horizontal_line = 0)

        exp.mean_props_movies('phi_rw_min_val','phi_lw_min_val','mean_wing','mean_front_phi')
        exp.mean_props_movies('phi_rw','phi_lw','mean_wing','mean_mean_phi')

        plot_and_save_movie_prop(exp,yprop = 'mean_front_phi',wing_body = 'mean_wing')
        plot_and_save_movie_prop(exp,yprop = 'mean_mean_phi',wing_body = 'mean_wing')
