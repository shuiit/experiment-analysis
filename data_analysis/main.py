import numpy as np
import h5py
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as colormap
import plotly.graph_objects as go
from Experiment import Experiment
from Movie import Movie
from scipy.spatial.transform import Rotation as R

def rotation_matrix(yaw,pitch,roll):
    roll_mat = np.array([[1,0,0],[0 ,np.cos(roll),-np.sin(roll)],[0, np.sin(roll), np.cos(roll)]])
    pitch_mat = np.array([[np.cos(pitch),0,np.sin(pitch)],[0, 1,0],[-np.sin(pitch), 0, np.cos(pitch)]])
    yaw_mat = np.array([[np.cos(yaw),-np.sin(yaw),0],[np.sin(yaw),np.cos(yaw),0],[0, 0, 1]])
    return yaw_mat,pitch_mat,roll_mat

def fly_axes(movie101):
    
    xbody = np.array([movie101.vectors['X_x_body'],movie101.vectors['X_y_body'],movie101.vectors['X_z_body']]).T[0]
    ybody = np.array([movie101.vectors['Y_x_body'],movie101.vectors['Y_y_body'],movie101.vectors['Y_z_body']]).T[0]
    zbody = np.array([movie101.vectors['Z_x_body'],movie101.vectors['Z_y_body'],movie101.vectors['Z_z_body']]).T[0]
    axes_mat =np.vstack([xbody,ybody,zbody])
    return axes_mat

def plot_axes(ax,axes_mat,color = ['r','b','g']):
    
    ax.plot3D([0,axes_mat[0][0] + 0.1],[0,axes_mat[0][1] + 0.1],[0,axes_mat[0][2] + 0.1],color = color[0])
    ax.plot3D([0,axes_mat[1][0] + 0.1],[0,axes_mat[1][1] + 0.1],[0,axes_mat[1][2] + 0.1],color = color[1])
    ax.plot3D([0,axes_mat[2][0] + 0.1],[0,axes_mat[2][1] + 0.1],[0,axes_mat[2][2] + 0.1],color = color[2])
    return ax


if __name__ == '__main__':  

    loadir = 'H:/My Drive/dark 2022/csv_dark/' # dir to save the data frame in # H:\My Drive\Ronis Exp\sagiv\data_analysis
    exp_name = ['2022_05_19_40ms']#['2022_02_03','2022_03_10','2022_05_19','2022_03_03'] # 22_11_28
    
    experiment = Experiment(loadir,exp_name[0])
    movie101 = Movie(experiment.experiment,'mov2')
    movie101.get_strokes()
    movie101.project_body_props_on_xy(['CM_real_x_body','CM_real_y_body','CM_real_z_body'],'CM_xy')
    
    angles = np.array([movie101.body['yaw_body'],-movie101.body['pitch_body'],movie101.body['roll_body']]).T[0]
    yaw_mat,pitch_mat,roll_mat = rotation_matrix(angles[0]*np.pi/180,angles[1]*np.pi/180,angles[2]*np.pi/180)
    axes_mat = fly_axes(movie101)
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_aspect('equal')
    ax = plot_axes(ax,axes_mat)
    
    
    axes_rotated = (yaw_mat @ pitch_mat @ roll_mat  ).T @ axes_mat.T
    axes_rotated = axes_rotated/np.linalg.norm(axes_rotated,axis = 1)
    color = ['magenta','cyan','lime']
    ax = plot_axes(ax,axes_rotated,color)
    ax.set_xlim(-1,1), ax.set_xlabel('x')
    ax.set_ylim(-1,1), ax.set_ylabel('y')
    ax.set_zlim(-1,1), ax.set_zlabel('z')



    group_name = 'body_angles'
    color_map = colormap.datad["tab10"]['listed']
    size = 5
    mov_name = 'mov101'

    fig = go.Figure() 

    movs = ['mov101', 'mov103']
    prop_xy = ['time','pitch_body']
    mov_name = movs[0]
    legend_name = experiment.experiment_name
    color = color_map[0]
    experiment.plot_mov(fig,mov_name,group_name,prop_xy,color,legend_name)


    x,y = experiment.get_data(mov_name,group_name,prop_xy)
    experiment.plot_prop_movie(fig,x,y,
                        color_map[0],mov_name, size = 5,showlegend = True)
