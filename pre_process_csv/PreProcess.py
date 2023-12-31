import pandas as pd
import matplotlib.pyplot as plt



class PreProcess:

    def __init__(self,exp_dir,exp_name,mov):
        self.mov_path = f'{exp_dir}/{exp_name}/{exp_name}_mov_{mov}'
        self.load_csv()

    def load_csv(self):
        self.angles = pd.read_csv(f'{self.mov_path}_angles.csv')
        self.vectors = pd.read_csv(f'{self.mov_path}_vectors.csv')
        self.ew_to_lab_rotmat = pd.read_csv(f'{self.mov_path}_ew_to_lab_rotmat.csv', header = None)

    def manual_clip_frames(self):
        fig, ax = plt.subplots()
    
        phi,  = ax.plot(self.angles['frames_frame'],self.angles['phi_rw'])
        ax_twin = ax.twinx()

        pitch,  = ax_twin.plot(self.angles['frames_frame'],self.angles['pitch_body'],color = 'red')
        zero_frame = self.angles['frames_frame'][self.angles['timeframe_time']==0]
        if len(zero_frame)>0: ax_twin.axvline(x = zero_frame.iloc[0], color = 'red',linewidth = 5)
        ax_twin.autoscale()
        ax_twin.relim()
       
        x = plt.ginput(2)
        self.angles = self.angles[(self.angles['frames_frame']>=int(x[0][0])) & (self.angles['frames_frame']<=int(x[1][0]))]
        self.vectors = self.vectors[(self.vectors['frames_frame']>=int(x[0][0])) & (self.vectors['frames_frame']<=int(x[1][0]))]
