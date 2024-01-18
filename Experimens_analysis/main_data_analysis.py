import matplotlib.pyplot as plt
from ExperimentCollection import ExperimentCollection
import numpy as np
from Plotters import Plotter
import utilities
from Experiment import Experiment
import h5py





if __name__ == '__main__':  
    
    # Define the necessary variables
    
    loadir = 'H:\\My Drive\\dark 2022\\' # dir to save the data frame in # H:\My Drive\Ronis Exp\sagiv\data_analysis
    exp_date = ['2022_02_03']#['2022_02_03','2022_03_10','2022_05_19','2022_03_03'] # 22_11_28
    exp_name = ['Led up']#['Led up','20 ms pert','40 ms pert','220Hz']
    pert_time = [[0, None],[0,20],[0,40],[0,None]]
    
    features_to_keep = ['frames', 'time [ms]', 'pitch', 'yaw', 'yaw_zer', 'roll', 'X', 'Y', 'Z', 'mov',
            'pitch smth', 'pitch smth_dot', 'pitch smth_dot_dot',
            'roll smth', 'roll smth_dot', 'roll smth_dot_dot',
            'yaw smth', 'yaw smth_dot', 'yaw smth_dot_dot',
            'Z smth','Z smth_dot_dot', 'phi_rw', 'phi_lw', 'theta_rw', 'psi_rw','phi_rw smth', 'theta_rw smth', 
            'psi_rw smth', 'phi_lw smth',   'theta_lw smth','psi_lw smth', 
            'p','q', 'r', 'p_dot', 'q_dot', 'r_dot', 'projected_x', 'projected_x_dot',
            'projected_x_dot_dot', 'projected_y', 'projected_y_dot','projected_y_dot_dot']
    
    features_to_keep_mean = features_to_keep + ['max_phi_lw','min_phi_lw','max_phi_rw','min_phi_rw']
    
    manual_del = {'2022_02_03':['mov6','mov7','mov65','mov59','mov57','mov36'],'2022_03_10':['mov37','mov32']}
    input_del_short_mov = [1500, 0]
    input_filter_manouvering = {'manouver_prop':'pitch smth_dot','time_of_maneuvers':[-5,-2,0],'prop_threshold':5*10**5,'plot_delted':False}
    
    # Create an instance of ExperimentCollection
    collection = ExperimentCollection()
#%%
     # Add experiments to the collection and loads it
    utilities.load_experiments(collection,exp_date,exp_name,pert_time,loadir)
    # Clean data for each experiment in the collection
    utilities.clean_data(collection,input_del_short_mov ,manual_del ,input_filter_manouvering)
    # Manipulate data for each experiment in the collection
    utilities.manipulate_data(collection,features_to_keep,features_to_keep_mean)
    # Calculate interest points for each experiment in the collection
    utilities.calculate_interest_points(collection)

#%%  
# =============================================================================
# plotting
# =============================================================================
  
    # Plot histograms of all experiments
    figures = [None]*8
    nbinsx = 20
    plot_flies_list = [Plotter(experiment) for experiment_id, experiment in collection.iterate_experiments()]
    hist_to_plot = ['zero_velocity','minimal_velocity','maximal_velocity','minimal_acc','maximal_acc','maximal_pitch']
    
    #%%
    
    
    prop = 'projected_x_dot'
    led_up = collection.get_experiment('Led up')
    led_up.mean_data['mean_front_phi'] = (led_up.mean_data['min_phi_lw'] + led_up.mean_data['min_phi_rw'])/2
    led_up.mean_data['mean_phi'] = (led_up.mean_data['phi_lw'] + led_up.mean_data['phi_rw'])/2

    led_up.interest_points['maximal_phi'] = led_up.mean_data.groupby('mov')['mean_phi'].idxmax() # find min to use when there are no sign changes. No time constrain


    plot_led_up = Plotter(led_up)
    plot_led_up.plot_props_per_movie(prop, propx='time [ms]', alpha=0.5, mean_data=True, scatter_kwargs={'mode': 'markers+lines'}, layout_kwargs={'title': '<b>' + plot_flies_list[0].exp_name + '<b>'})
    plot_led_up.plot_props_per_movie('mean_front_phi', propx='time [ms]', alpha=0.5, mean_data=True, scatter_kwargs={'mode': 'markers'}, layout_kwargs={'title': '<b>' + plot_flies_list[0].exp_name + '<b>'})
    plot_led_up.plot_props_per_movie('mean_phi', propx='time [ms]', alpha=0.5, mean_data=True, scatter_kwargs={'mode': 'markers+lines'}, layout_kwargs={'title': '<b>' + plot_flies_list[0].exp_name + '<b>'})

    
    #%%
    
    for plot_flies in plot_flies_list:     
        for i,prop in enumerate(hist_to_plot):
            figures[i] = plot_flies.idx_histogram( prop,fig = figures[i],layout_kwargs = {'title' : f'<b>Time of {prop} [ms]<b>'} ,nbinsx = nbinsx)
        
        figures[i+1] = plot_flies.histogram( 'delta_v',fig = figures[6],title = '<b>' +'Delta velocity [m/s]' + '<b>',nbinsx = nbinsx)
        figures[i+2] = plot_flies.histogram( 'delta_pitch',fig = figures[7],title = '<b>' +'Delta pitch [deg]' + '<b>',nbinsx = nbinsx)
        

    for fig in figures:
        fig.show()
#%%
    fig = None
    for plot_flies in plot_flies_list:
        fig = plot_flies.plot_props_per_experiment('projected_x_dot_dot',propx = 'time [ms]',movie_name = 'all',
                                  alpha= 0.5,size = 5, return_fig = True, mean_data = False, fig = fig)
    fig.show()   
#%%     
    # Plot min/max points per movie for one experiment
    for plot_flies in plot_flies_list[2:4]:
        plot_flies.plot_min_on_movie('pitch smth', 'maximal_pitch', title='<b>' + plot_flies.exp_name + '<b>')
        plot_flies.plot_min_on_movie('projected_x_dot', 'minimal_velocity', title='<b>' + plot_flies.exp_name + '<b>')
        plot_flies.plot_min_on_movie('projected_x_dot_dot', 'minimal_acc', title='<b>' + plot_flies.exp_name + '<b>')
    
    # Plot all movies of an experiment according to user-defined property
    prop = 'max_phi_lw'
    plot_flies_list[0].plot_props_per_movie(prop, propx='time [ms]', alpha=0.5, mean_data=True, scatter_kwargs={'mode': 'lines+markers'}, layout_kwargs={'title': '<b>' + plot_flies_list[0].exp_name + '<b>'})

#%%
    import scipy.signal as signal
    
    mov5 = Experiment.get_mov(led_up.data, ['mov15'])
    dt = (mov5['time [ms]'].iloc[1] - mov5['time [ms]'].iloc[0])/1000
    pitch_dot = signal.savgol_filter(mov5['pitch']/dt, 600, 3,deriv = 1)
    pitch_dot_dot = signal.savgol_filter(mov5['pitch']/dt/dt, 300, 3,deriv = 2)

    fig, axs = plt.subplots(3, 1,sharex=(True))
    axs[0].plot(mov5['time [ms]'],mov5['pitch'])
    axs[1].plot(mov5['time [ms]'],pitch_dot)
    axs[2].plot(mov5['time [ms]'],pitch_dot_dot)

#%%     
    # Plot a multi Y-axis for a specific movie and then plot its trajectory
    movie_name = 'mov11' #'mov5
    props = ['pitch smth_dot', 'pitch smth_dot_dot','pitch smth']
    plot_flies_list[0].plot_multi_y(movie_name, props)
    plot_flies_list[0].plot3d_traj(movie_name)
    
#%%
    # Plot all histograms of the same experiment together
    figures = [None]*1
    nbinsx = 40
    trace_kwargs = {'opacity': 0.9}
    title = '<b>' + plot_flies.exp_name + '<b>'
    for plot_flies in plot_flies_list[0:1]:
        figures[0] = plot_flies.idx_histogram('zero_velocity', fig=figures[0], layout_kwargs={'title': title}, trace_kwargs=trace_kwargs, nbinsx=nbinsx)
        figures[0] = plot_flies.idx_histogram('minimal_velocity', fig=figures[0], layout_kwargs={'title': title}, trace_kwargs=trace_kwargs, nbinsx=nbinsx)
        figures[0] = plot_flies.idx_histogram('minimal_acc', fig=figures[0], layout_kwargs={'title': title}, trace_kwargs=trace_kwargs, nbinsx=nbinsx)
        figures[0] = plot_flies.idx_histogram('maximal_acc', fig=figures[0], layout_kwargs={'title': title}, trace_kwargs=trace_kwargs, nbinsx=nbinsx)      
        figures[0] = plot_flies.idx_histogram('maximal_pitch', fig=figures[0], layout_kwargs={'title': title}, trace_kwargs=trace_kwargs, nbinsx=nbinsx)
        
    for fig in figures:
        fig.show()


    wakk = 2
