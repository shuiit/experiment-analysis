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
    exp_name = ['manipulated_2022_05_19_40ms'
                ]#['2022_02_03','2022_03_10','2022_05_19','2022_03_03'] # 22_11_28
    exp_name = ['manipulated_2023_08_10_100ms','manipulated_2022_02_03_dark','manipulated_2023_08_06_40ms',
                'manipulated_2023_08_09_60ms','manipulated_2023_08_07_5ms','manipulated_2023_08_07_10ms'
                ]#['2022_02_03','2022_03_10','2022_05_19','2022_03_03'] # 22_11_28
    experiments = [Experiment(loadir,exp) for exp in exp_name]

    min_max_v = [exp.min_max_point('CM_xy_dot') for exp in experiments]
    min_max_pitch =[exp.min_max_point('pitch_y_frame') for exp in experiments]
    zero_velocity = [exp.zero_velocity('CM_xy_dot') for exp in experiments]


    from Plotters import Plotters
    exp = experiments[0]
    
    mov = exp.get_mov('mov6')
    data = mov.data
    header = mov.header
    int(exp.pertubation_name.split('ms')[0])
    data,plot_cofnig = mov.calculation_for_3d_traj()
    ploter = Plotters(data,False)
    ploter.plot_3d_traj(data,plot_cofnig)



    def mean_props(mov,prop1,prop2,wing_body,header_name):
        mean_prop = (mov.get_prop(prop1,wing_body)  + mov.get_prop(prop2,wing_body) )/2
        mov.data[wing_body] = np.hstack((mov.data[wing_body], mean_prop[np.newaxis,:].T))
        mov.add_to_header([header_name],wing_body)

    [mean_props(experiments[0].get_mov(mov_name),'phi_rw_min_val','phi_lw_min_val','mean_wing','mean_front_phi') for mov_name in experiments[0].mov_names]


    fig = go.Figure()
    [exp.plot_movies('time','mean_front_phi', case = 'plot_mov',prop = 'mean_wing', fig = fig, mov = False,color = idx + 2) for idx,exp in enumerate(experiments)]
    [exp.plot_movies('time','pitch_y_frame', case = 'plot_mov',prop = 'body', fig = fig, mov = False,color = idx + 2) for idx,exp in enumerate(experiments)]

    fig.show()


    
    fig = go.Figure()    
    [fig.add_trace(go.Histogram(x=np.hstack(min[0]), name=exp.experiment_name, nbinsx=10)) for min,exp in zip(min_max_v,experiments)]
    fig.update_layout(title = 'time of minimum velocity')
    fig.update_traces(opacity = 0.9)
    fig.show()


    fig = go.Figure() 
    [fig.add_trace(go.Histogram(x=np.hstack(min[1]), name=exp.experiment_name, nbinsx=10)) for min,exp in zip(min_max_pitch,experiments)]
    fig.update_layout(title = 'time of maximum pitch' )
    fig.show()


    fig = go.Figure() 
    [fig.add_trace(go.Histogram(x=np.hstack(min), name=exp.experiment_name, nbinsx=10)) for min,exp in zip(zero_velocity,experiments)]
    fig.update_layout(title = 'time of first zero velocity' )
    fig.show()




    fig = go.Figure()
    [exp.plot_movies('time','CM_xy_dot', case = 'plot_exp', fig = fig, mov = False,color = idx + 2) for idx,exp in enumerate(experiments)]
    fig.show()

    
    fig = go.Figure()
    [exp.plot_movies('time','pitch_y_frame', case = 'plot_exp', fig = fig, mov = False,color = idx + 2) for idx,exp in enumerate(experiments)]
    fig.show()

