import os
from glob import glob
import csv
from scipy import io
import pandas as pd
from functools import reduce
import operator
import numpy as np
import matplotlib.pyplot as plt

def generatefile_loc(pth):
  # pthfiles = [y for x in os.walk(pth) for y in glob(os.path.join(x[0], '*.mat').replace(r'\\',r'/'))]
    pthfiles = []
    for x in os.walk(pth):
      # for y in glob(os.path.join(x[0], '*.mat')):
      #     y.replace(r'\\', r'/')
      if x[0].split('\\')[-1] == 'hull_op':
          for files in os.walk(x[0]):
                for y in glob(os.path.join(files[0], '*.mat')):
                    pthfiles.append(y.replace(r'\\', r'/'))

    return pthfiles

def fileLoc2csv(pth,pthfiles,header):
  f = open(pth + "datasetFiles.csv", "w")
  writer = csv.DictWriter(f, fieldnames=header)
  writer.writeheader()
  for path in pthfiles:
    writer = csv.writer(f)
    writer.writerow([path + os.linesep,path[path.find('mov'):-1].split('\\')[0]])

def generateCSV(pth):
  pthfiles = generatefile_loc(pth)
  fileLoc2csv(pth,pthfiles,['location','mov'])


class fly():
    def __init__(self, annotations_file, img_dir, transform=None, target_transform=None, np_image=True):
        self.file_loc = pd.read_csv(annotations_file)
        self.path = img_dir
        self.movies = np.unique(self.file_loc['mov'].values)

    def __len__(self):
        return len(self.img_labels)

    def getitem(self, mov, filenm_part):
        # filenm_part: part of the file name (exp Smov for hull_Smov1)
        movfls = self.file_loc['location'][self.file_loc['mov'] == mov]
        movfile = movfls[movfls.str.contains(filenm_part)]
        indmat = movfile.iloc[0][0:-1].find('mat')
        mat = io.loadmat(movfile.iloc[0][0:indmat+3])
        self.movie = mov
        return mat

# arange matlab structure in dictionary
def struct2dict(dic, nms_out, nms_in, datMat):
    for nmstm in nms_in:
            if type(datMat.item()[nmstm].item().dtype.names) is tuple:

                names = datMat.item()[nmstm].item().dtype.names
                dic[nms_out].update({nmstm: {names[0]: np.real(datMat.item()[nmstm].item()[names[0]][0, 0])}})
                for nm in names[1:]:
                    dic[nms_out][nmstm].update({nm: datMat.item()[nmstm].item()[nm][0, 0]})
            else:
                dic[nms_out].update({nmstm: datMat.item()[nmstm].item()[0, 0]})
    return dic


# create dataframe from dictionary
def createDataFrame(hull, prop, prop2, timeV, frmsV):
    dfdict = {}
    for bodywing in prop:
        if 'body' not in bodywing:
            nms = list(hull[bodywing]['angles'].keys())
            for nm in nms:
                if 'sec' in nm:
                    hull[bodywing]['angles'][nm + 'LE'] = hull[bodywing]['angles'][nm]['LE'][0][0]
                    hull[bodywing]['angles'][nm + 'TE'] = hull[bodywing]['angles'][nm]['TE'][0][0]
                    del hull[bodywing]['angles'][nm]
                    
    
    for bodywing in prop:

        datf = pd.DataFrame.from_dict(hull[bodywing])
        datf = datf[prop2].dropna().to_frame()  
        df = pd.DataFrame({'time [ms]': timeV, 'frames': frmsV})  # add the time vector ato the data frame
        for nm in datf[prop2].keys():
            if nm in ['phi','theta','psi','pitch','yaw','roll'] or 'sec' in nm:
                val = reduce(operator.iconcat, datf[prop2][nm], [])
                
                if nm == 'yaw':
                    arrval = np.array(val)
                    arrval[~np.isnan(arrval)] = np.unwrap(arrval[~np.isnan(arrval)])
                    val = list(arrval)
                
                    
                if '_old' not in nm:
                    if len(val) > len(timeV) or len(val) > len(df):
                        val = val[0:len(df)]
                    if len(val) < len(timeV) or len(val) < len(df):
                        df = df[0:len(val)]
                if bodywing == 'rightwing':
                    nm = nm + '_rw'
                if bodywing == 'leftwing':
                    nm = nm + '_lw'
                
                
                df[nm] = val
        if bodywing == 'body':
            df['X'] = hull[bodywing]['coords']['CM_real'][:,0]
            df['Y'] = hull[bodywing]['coords']['CM_real'][:, 1]
            df['Z'] = hull[bodywing]['coords']['CM_real'][:, 2]

        del hull[bodywing][prop2]
        dfdict[bodywing] = df
    return hull, dfdict

def mat_of_movies(annotations_file,pth,st = 0,en = -1,movnames = 'all',delmov = 0):
    # generate a dictionary of body, rightwing and leftwing of multiple movies
    hull = {}
    body = []
    rwng = []
    lwng = []
    dict_angles =dict.fromkeys(['body','leftwing','rightwing'])
    count = 1
    datafly = fly(annotations_file, pth)  # initilize dataset class
    body_vectors = {}
    if movnames == 'all':
        movs = datafly.movies[st:en]
        if delmov != 0:
            for mv in delmov:
                movs = np.delete(movs,movs == mv)
    else:
        movs = movnames
    for nms in movs: # for loop on all movies
       
        try:
          movfls = datafly.file_loc['location'][datafly.file_loc['mov'] == nms]
          movfile = movfls[movfls.str.contains('Shull')]
          if len(movfile) == 0:
              continue
          mat  = datafly.getitem(nms,'Shull')
          ## create a dictionary for the 3D and general fields------------
          nms_out = mat['Shull'].dtype.names

          for nmefield in nms_out:
              if nmefield in ['sprs', 'wingana', 'bodyana']:
                  continue
              else:
                datMat = mat['Shull'][nmefield]
                nms_in = datMat.item().dtype.names
                if nms_in:
                    nms_in = list(nms_in)
                    if 'hull3d' in nms_in:
                        nms_in.remove('hull3d')
                if nms_in is None:
                    if nmefield == 'real_coord':
                        hull[nmefield] = mat['Shull'][nmefield].item()[0,:]
                    else:
                        hull[nmefield] = mat['Shull'][nmefield].item()
                elif nmefield == 'video':
                    hull.update({nmefield: {nms_in[0]: datMat.item()[nms_in[0]].item()}}, )
                    hull['video'][nms_in[1]] = datMat.item()[nms_in[1]].item()
                else:
                    hull.update({nmefield: {nms_in[0]:datMat.item()[nms_in[0]].item()}})
                    hull = struct2dict(hull,nmefield,nms_in,datMat)
        #----------------------------------------------------------------

        # calculate the time and frame vectors and create a pandas data frame for the angles---
          timeVec,frms = calcTimeVec(hull)
          timeV = pd.Series(timeVec[0])
          frmsV = pd.Series(frms)
          prop = dict_angles.keys()
          prop2 = 'angles'
          if count == 1: # initialize dataframe for body and wings
            bodmov = pd.DataFrame({'time [ms]': timeV, 'frames': frmsV})  # add the time vector ato the data frame
            Rwing = pd.DataFrame({'time [ms]': timeV, 'frames': frmsV})  # add the time vector ato the data frame
            Lwing = pd.DataFrame({'time [ms]': timeV, 'frames': frmsV})  # add the time vector ato the data frame
            count = 0
          hull,angles = createDataFrame(hull,prop,prop2,timeV,frmsV) # hull - 3d data and general field. angles - the angles data frame

          #-------------------------------------------------------------------------------------
          try:
# =============================================================================
#             angles['leftwing'] = angles['leftwing'].dropna()
#             angles['rightwing'] = angles['rightwing'].dropna()
#             angles['body'] = angles['body'].dropna()
# =============================================================================
            maxfrm = min([max(angles['leftwing']['frames']), max(angles['rightwing']['frames']), max(angles['body']['frames'])]) 
            minfrm = min(angles['leftwing']['frames']) 
            body_vectors[nms] = hull['body']['vectors']
            body_vectors[nms]['frames'] = hull['frames'][0]

            for nms_prop in prop:
                angles[nms_prop]['mov'] = nms
            #angles[nms_prop] = angles[nms_prop].set_index('frames', drop = False).dropna()
                angles[nms_prop] = angles[nms_prop].set_index('frames', drop = False)
            body.append(angles['body'].loc[minfrm :maxfrm ])
            rwng.append(angles['rightwing'].loc[minfrm :maxfrm ])
            lwng.append(angles['leftwing'].loc[minfrm :maxfrm ])
# =============================================================================
#                 dict_angles[nms_prop] = pd.concat([dict_angles[nms_prop],angles[nms_prop].loc[minfrm :maxfrm ]],ignore_index = True)# generate body angles dataframe
# =============================================================================
    
          except:
            continue
        except:
          continue
    dict_angles['body'] = pd.concat(body)
    dict_angles['rightwing'] = pd.concat(rwng)
    dict_angles['leftwing'] = pd.concat(lwng)

    return dict_angles,body_vectors


def calcTimeVec(hull):
    # calculate time vector

    Xvec = hull['frames'][0]
    XvecOP = Xvec;
    XvecTime = hull['video']['timeframe']
    return XvecTime, Xvec


def getDF_mov(angdict,nameOfMov):
    if nameOfMov != '':
        rows = []
        for nmmov in nameOfMov:
            rowsOfMov = np.where(angdict['body']['mov'] == nmmov)[0]
            rows = np.append(rows,rowsOfMov)
        body = angdict['body'].iloc[rows]
        rw = angdict['rightwing'].iloc[rows]
        lw = angdict['leftwing'].iloc[rows]
    else:
        body = angdict['body']
        rw = angdict['rightwing']
        lw = angdict['leftwing']
    return body,rw,lw