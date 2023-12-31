from IPython import get_ipython
get_ipython().magic('reset -sf')
import utilities
import loaders
import os
import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import time
from tqdm import tqdm
import pickle
import matplotlib as mpl
from sklearn.decomposition import PCA
import re
from pathlib import Path

# =============================================================================
# User input
# =============================================================================
expname = '22_12_05' # name of experiment 22_12_05
# =============================================================================
# expname = '22_11_28_1630-2359'

# =============================================================================
expdir = 'H:\\My Drive\\Ronis Exp\\sagiv\\' # dir to save the data frame in
pth = expdir + expname + '\\hull\\hull_Reorder\\' # dir of the specific experiment
savedir = 'H:\\My Drive\\Ronis Exp\\sagiv\\allFrames2\\'
# =============================================================================
# savedir = 'H:\\My Drive\\Ronis Exp\\sagiv_allframes\\' # dir to save the data frame in
# =============================================================================






savedir = 'H:\\My Drive\\Ronis Exp\\sagiv\\allFrames2\\'




# =============================================================================
# expname = '2022_02_03' # name of experiment 
# expdir = 'H:\\My Drive\\dark 2022\\' # dir to save the data frame in
# pth = expdir + '\\2022_02_03\\light021220_RenameV2\\light021220_RenameV2_Reorder\\' # dir of the specific experiment
# =============================================================================

# =============================================================================
# expname = '2022_01_31' # name of experiment 
# expdir = 'H:\\My Drive\\dark 2022\\' # dir to save the data frame in
# pth = expdir + '\\2022_01_31\\hullAn\\hullAn_Reorder\\' # dir of the specific experiment
# 
# =============================================================================
# =============================================================================
# expname = '2022_03_10' # name of experiment 
# expdir = 'H:\\My Drive\\dark 2022\\' # dir to save the data frame in
# pth = expdir + '\\2022_03_10\\4cams\\hull\\hull_Reorder\\' # dir of the specific experiment
# 
# =============================================================================
expname = '22_11_28' # name of experiment 
expdir = 'H:\\My Drive\\Ronis Exp\\sagiv\\' # dir to save the data frame in
pth = expdir + expname + '\\hull\\hull_Reorder\\' # dir of the specific experiment





allframes = 1
smoothwinbody = 211

loaders.generateCSV(pth) # generate the CSV (the movies and their location)
annotations_file = os.path.join(pth , "datasetFiles.csv") # location of csv file

file_loc = pd.read_csv(annotations_file)
movies = np.unique(file_loc['mov'].values)

if allframes == 1:
    nameOfMov = 'all'
    nameOfMov_all = ['mov' + str(i) for i in range(0,len(movies)) ]
    mvchnk = np.hstack((np.arange(1,len(movies),100),len(movies)))

else:
    frames_file = os.path.join(pth , "Exp_sam.txt") # open the detailed txt file with the initial and end frames of the experiment
    
    goodframes = pd.read_csv(frames_file, sep="|",comment=('#')).dropna().reset_index(drop=True)
    nameOfMov = []                                                # specific: ['mov1','mov2'...] or [nameOfMov[0],nameOfMov[1]...]
    # =============================================================================
    # for i in goodframes.index:
    # # =============================================================================
    # #     if abs(int(goodframes['frms0'][i]) - int(goodframes['frms1'][i])) < 600:
    # #         continue
    # # =============================================================================
    #     nameOfMov.append('mov' + str(goodframes['mov'][i]))
    # =============================================================================
    nameOfMov_all = ['mov' + str(i) for i in goodframes['mov'] ] # movies to open. all: ['mov' + str(i) for i in goodframes['mov'] ], 


    mvchnk = np.hstack((np.arange(1,len(nameOfMov_all),100),len(nameOfMov_all)))

for i in range(0,len(mvchnk)-1):
    nameOfMov = nameOfMov_all[mvchnk[i]:mvchnk[i + 1]]                                            
    expname_sv =  expname +'_' +  nameOfMov_all[mvchnk[i]] +'_' + nameOfMov_all[mvchnk[i + 1]-1]                                                  
    
    
    # =============================================================================
    # create basic dataframe - body/wing angles, [X Y Z]_lab vector. smooth and diff all variables
    # =============================================================================
    angdict,body_vectors = loaders.mat_of_movies(annotations_file,pth,movnames = nameOfMov)
    angdict['body'] = angdict['body'].dropna(subset = ['pitch', 'roll', 'yaw','X','Y','Z']).reset_index(drop=True)
    angdict['rightwing'] = angdict['rightwing'].dropna(subset = ['phi_rw','theta_rw', 'psi_rw']).reset_index(drop=True)
    angdict['leftwing'] = angdict['leftwing'].dropna(subset =  ['phi_lw', 'theta_lw', 'psi_lw']).reset_index(drop=True)
    
        
    for nms in ['pitch', 'roll', 'yaw','X','Y','Z']:
        angdict['body'].groupby('mov').apply(lambda x: utilities.interpOLAndSmoothV2(angdict['body'], nms,smthwin=smoothwinbody, smthpol=4)) 
    for nms in ['phi', 'theta', 'psi']:    
        angdict['rightwing'].groupby('mov').apply(lambda x: utilities.interpOLAndSmoothV2(angdict['rightwing'], nms + '_rw',smthwin=7, smthpol=2)) 
        angdict['leftwing'].groupby('mov').apply(lambda x: utilities.interpOLAndSmoothV2(angdict['leftwing'], nms + '_lw',smthwin=7, smthpol=2)) 
    
       
    # =============================================================================
    # trim frames dfrm_min and dfrm_max from the beggining and end of the movie
    # =============================================================================
    
    
    #!! CAN MAKE IT BETTER - IM IN A RUSH... makes sure that the sa,me frames are in both wings and body dictionary
    indx = []
    pdrw = []
    pdlw = []
    pdbod = []
    u,exp2calc = np.unique(angdict['body']['mov'], return_index=True)
    mvorder = u[np.argsort(exp2calc)]
    utilities.addStrks2DF(angdict['rightwing'],debugFBstrk = 0)
    utilities.addStrks2DF(angdict['leftwing'],debugFBstrk = 0)
    

    # =============================================================================
    # Calculate velocity, acceleration and mark strokes. Concatenate angdict to get one dataframe - trans_df
    # =============================================================================
    # add forward and backward strokes, mark forward as 0 and backward as 1
   
    
    
    
    # devide by dt to get velocity and acceleration. calculate body velocity and acceleration

    trans_df = pd.merge(angdict['body'], angdict['rightwing'], on=['mov','frames'], how="left", indicator=True)
    trans_df = trans_df.groupby('_merge').get_group('both')
    trans_df = trans_df.drop(columns=['_merge'])
    trans_df = pd.merge(trans_df, angdict['leftwing'], on=['mov','frames'], how="left", indicator=True)
    trans_df = trans_df.groupby('_merge').get_group('both')
    trans_df = trans_df.drop(columns=['_merge'])
    trans_df = trans_df.drop(columns=['time [ms]_x'])

    
    
    # angular accelerations and angular velocity
    dt = np.diff(angdict['body']['time [ms]'])[0]/1000
    trans_df[['yaw smth_dot_dot','pitch smth_dot_dot','roll smth_dot_dot']] = trans_df[['yaw smth_dot_dot','pitch smth_dot_dot','roll smth_dot_dot']]/dt**2
    trans_df[['yaw smth_dot','pitch smth_dot','roll smth_dot']] = trans_df[['yaw smth_dot','pitch smth_dot','roll smth_dot']]/dt
    trans_df['expname'] = expname
    
    
    
    # calculate the body angular acceleratrion and velocity: pqr and pqr_dot
    utilities.pqr_pqr_dot(trans_df)
    
    Path(savedir + 'data_analysis').mkdir(parents=True, exist_ok=True) # create dir to save the output
    with open(savedir + 'data_analysis\\' + expname_sv + '_angles_dict.pkl', 'wb') as f:
        pickle.dump(angdict, f)
    with open(savedir + 'data_analysis\\' + expname_sv + '_body_vec.pkl', 'wb') as f:
        pickle.dump(body_vectors, f)
    with open(savedir + 'data_analysis\\' + expname_sv + '_DF_all.pkl', 'wb') as f:
        pickle.dump(trans_df, f)
        
        
    # =============================================================================
    fig,axes =plt.subplots(3,1, figsize=(12, 9)) # 3 columns each containing 10 figures, total 30 features
    
    axes[0].plot(trans_df['roll']),axes[2].set_title('roll')
    axes[1].plot(trans_df['yaw']),axes[2].set_title('yaw')
    axes[2].plot(trans_df['pitch']),axes[2].set_title('pitch')
    
    
    fig,axes =plt.subplots(3,1, figsize=(12, 9),sharex=True) # 3 columns each containing 10 figures, total 30 features
    
    axes[1].plot(trans_df['frames'],trans_df['phi_rw smth_dot'])
    axes[1].plot(trans_df['frames'],trans_df['phi_lw smth_dot'])
    
    axes[0].plot(trans_df['frames'],trans_df['phi_rw smth'])
    axes[0].plot(trans_df['frames'],trans_df['phi_lw smth'])
    
    
    axes[2].plot(trans_df['frames'],trans_df['phi_rw smth_dot_dot'])
    axes[2].plot(trans_df['frames'],trans_df['phi_lw smth_dot_dot'])


