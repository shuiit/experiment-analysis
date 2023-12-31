import pandas as pd
import matplotlib.pyplot as plt



class PreProcess:

    def __init__(self,exp_dir,exp_name,mov):
        self.mov_path = f'{exp_dir}/{exp_name}/{exp_name}_mov_{mov}'
        self.load_csv()

    def load_csv(self):
        self.angles = pd.read_csv(f'{self.mov_path}_angles.csv')
        self.vectors = pd.read_csv(f'{self.mov_path}_vectors.csv')
        self.ew_to_lab_rotmat = pd.read_csv(f'{self.mov_path}_ew_to_lab_rotmat.csv')

    def manual_clip_frames(self):
        fig, ax = plt.subplots()

        ax.plot(self.angles['frames_frame'],self.angles['phi_rw'])
        ax.plot(self.angles['frames_frame'],self.angles['pitch_body'])
        x = plt.ginput(2)
        self.angles = self.angles[(self.angles['frames_frame']>=int(x[0][0])) & (self.angles['frames_frame']<=int(x[1][0]))]
        self.vectors = self.vectors[(self.vectors['frames_frame']>=int(x[0][0])) & (self.vectors['frames_frame']<=int(x[1][0]))]
        
