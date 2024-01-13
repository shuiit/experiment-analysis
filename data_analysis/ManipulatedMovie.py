import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px
import matplotlib.cm
from numpy import linalg 

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from plotly.subplots import make_subplots
import h5py
from scipy.signal import argrelextrema, savgol_filter,find_peaks
from scipy.spatial.transform import Rotation as R



pio.renderers.default='browser'

class ManipulatedMovie():
    def __init__(self,experiment,mov_name):    
        self.mov = {}
        self.data = {dict_name:np.array(experiment[mov_name][dict_name]) for dict_name in experiment[mov_name].keys()}
        self.header = {dict_name:self.get_header(experiment[mov_name][dict_name]) for dict_name in experiment[mov_name].keys()}
        self.name = mov_name
    

    def get_header(self,dataset):
        return { header: idx for idx,header in enumerate(dataset.attrs['header'])}
    
    
    def add_to_header(self, string_to_add,dict_name):
        [self.header[dict_name].update({name:len(self.header[dict_name])}) for name in string_to_add]
