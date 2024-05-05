import pandas as pd
import Utils
import plotly


if __name__ == '__main__': 

    # the header of the CSV file is assumed to have the following header: 
    # CoM_x,CoM_y,CoM_z, x_body_x, x_body_y, x_body_z, y_body_x, y_body_y, y_body_z, time_ms, "color prop"
    # color prop is dwtermined by the user and should be a column in the CSV file


    csv_path = 'C:/Users/Roni/Documents/experiment-analysis/3d_trajectory/tst_mov2.csv'
    input_data = pd.read_csv(csv_path)


    data = {}
    color_prop = 'x_vel'
    mov_name = 'mov1'
    exp_name = 'dark mosquito'
    save_plot = True
    figures_path = 'I:/My Drive/Research/Dark/article/figures/Noams'
    plot_cofnig = {'fly_samples':150,'traj_samples':20,'size_x':1,'size_y':1/3,'delta_y_on_x':3/4}
    pertubation = 60 # if False: there is no pertubation

    data = Utils.prepare_data(input_data,color_prop,plot_cofnig,pertubation)
    fig = Utils.plot_3d_traj(data,plot_cofnig,mov_name,exp_name,color_prop)
    plotly.offline.plot(fig, filename=f'{figures_path}/traj_3d_{mov_name}.html',auto_open=False) if save_plot == True else fig.show()



