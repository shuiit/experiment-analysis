from scipy.signal import savgol_filter
import numpy as np
import pandas as pd
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from scipy.signal import savgol_coeffs


def zscore(x, window):
    # calculate rolling zscore
    # inp: x = data
    #      window = size of window for rolling
    r = x.rolling(window=window)
    m = r.mean().shift(1)
    s = r.std(ddof=0).shift(1)
    z = (x - m) / s
    return z

def interpOLAndSmoothV2(dat, ColRep, window=200, TH=3, mth='cubic', smthwin=15, smthpol=2):
    # find outliers using Zscore, interpulate the OL
    # use savgol filter to smooth
    # inp: dat = data from data frame
    #      ColRep = the name of the relavent column (the data column)
    #      window = window size to calculate zscore (7)
    #      TH = threashold of zscore (10)
    #      mth = method used for interpulation (cubic)
    #      smthwin = window used for smothing (15)
    #      smthpol = polynom degree for savgol filter (2)
    #      returns: data and new columns names
    if 'phi' in ColRep:
        dat[ColRep][dat[ColRep]>180] = 360-dat[ColRep]
    zscr = zscore(dat[ColRep], window)
    newCol_interp = ColRep + ' inter nan'
    newCol_smooth = ColRep + ' smth'
    newCol_dot = ColRep + ' smth' + '_dot'
    newCol_dot_dot = ColRep + ' smth' + '_dot'+ '_dot'

    dat[newCol_interp] = np.where(abs(zscr) > TH, float('NaN'), dat[ColRep])
    dat[newCol_interp] = dat[newCol_interp].interpolate(method=mth)
    if dat[newCol_interp].isna().eq(True).any().any():
      extrapidx = np.where(dat[newCol_interp].isna().eq(False).any())
      dat[newCol_interp] = dat[newCol_interp].interpolate(method='nearest', axis=0).ffill().bfill()  

    dat[newCol_smooth] = savgol_filter(dat[newCol_interp], smthwin, smthpol)
    dat[newCol_dot] = savgol_filter(dat[newCol_interp], smthwin, smthpol,deriv = 1)
    dat[newCol_dot_dot] = savgol_filter(dat[newCol_interp], smthwin, smthpol, deriv = 2)
    return dat

def interpOLAndSmooth(dat, ColRep, window=7, TH=10, mth='cubic', smthwin=15, smthpol=2):
    # find outliers using Zscore, interpulate the OL
    # use savgol filter to smooth
    # inp: dat = data from data frame
    #      ColRep = the name of the relavent column (the data column)
    #      window = window size to calculate zscore (7)
    #      TH = threashold of zscore (10)
    #      mth = method used for interpulation (cubic)
    #      smthwin = window used for smothing (15)
    #      smthpol = polynom degree for savgol filter (2)
    #      returns: data and new columns names

    zscr = zscore(dat[ColRep], window)
    newCol_interp = ColRep + ' inter nan'
    newCol_smooth = ColRep + ' smth'
    newCol_dot = ColRep + ' smth' + '_dot'
    newCol_dot_dot = ColRep + ' smth' + '_dot'+ '_dot'

    dat[newCol_interp] = np.where(abs(zscr) > TH, float('NaN'), dat[ColRep])
    dat[newCol_interp] = dat[newCol_interp].interpolate(method=mth)

    dat[newCol_smooth] = savgol_filter(dat[newCol_interp], smthwin, smthpol)
    dat[newCol_dot] = savgol_filter(dat[newCol_interp], smthwin, smthpol,deriv = 1)
    dat[newCol_dot_dot] = savgol_filter(dat[newCol_interp], smthwin, smthpol, deriv = 2)
    return dat, newCol_interp, newCol_smooth

def bin_prop(prop,bodmov,bins):
    for kprp in range(0, len(prop)):
        count_pitch, division_pitch = np.histogram(bodmov[prop[kprp]], bins=bins[kprp])
        bined_angles = pd.DataFrame()  # add the time vector ato the data frame
        bined_angles[prop[kprp]] = division_pitch
        pitchind = np.digitize(bodmov[prop[kprp]], division_pitch, right=True)
        bodmov['ind' + prop[kprp]] = pitchind
    return bodmov


def create_state_dictionary(statesName,bodmov):
    states = bodmov[statesName].drop_duplicates().to_numpy()
    indrepdict = {}
    npinter = bodmov[statesName].to_numpy()

    for rowtmp in range(0, np.size(states, 0)):  # len(indrep)
        a = np.where((npinter == states[rowtmp, :]).all(axis=1))
        if str(states[rowtmp]) not in indrepdict:
            indrepdict[str(states[rowtmp])] = a
    return indrepdict, states


def create_transitionMat(a1,states, indrepdict):
    transitionMat = np.zeros([len(states), len(states)])
    statesBatch = states[a1[0]:a1[1]]
    for kstate in range(0, np.size(statesBatch, 0)):
        inds1 = np.array(indrepdict[str(statesBatch[kstate])]).astype(int)
        for inkstate in range(0, np.size(states, 0)):
            inds2 = np.array(indrepdict[str(states[inkstate])]).astype(int)

            a = np.hstack([inds1, inds2 - 1])
            unique, counts = np.unique(a, return_counts=True)
            value_cou_eq_win = unique[counts == 2]
            transitionMat[kstate + a1[0], inkstate] = len(value_cou_eq_win)
    return transitionMat

def create_transitionMatV2(a1, indrepdict,allstates):
    transitionMat = np.zeros([len(indrepdict.keys()), len(indrepdict.keys())])
    allstates_batch = allstates[a1[0]:a1[1]]
    for k in range(0,len(allstates_batch)-1):
        key1 = str(allstates_batch[k])
        key2 = str(allstates_batch[k+1])
        inds1 = np.array(indrepdict[key1]).astype(int)
        inds2 = np.array(indrepdict[key2]).astype(int)
        value_cou_eq_win = findrep(inds1,inds2 - 1)
        transitionMat[list(indrepdict.keys()).index(key2), list(indrepdict.keys()).index(key1)] = len(value_cou_eq_win)
    return transitionMat

def findrep(inds1,inds2):
    a = np.hstack([inds1, inds2])
    unique, counts = np.unique(a, return_counts=True)
    value_cou_eq_win = unique[counts == 2]
    return value_cou_eq_win

def find_FBstrk_updateDF(wing,wn,plot = 0):
    # find forward and backward stroke. update in data frame and return indices
    # 0 - forward stroke
    # 1 - backward stroke
    
    phismth4maxmin = savgol_filter(wing['phi' + wn], 30, 2)
    indofmax = argrelextrema(phismth4maxmin, np.greater, order=30)
    indofmin = argrelextrema(phismth4maxmin, np.less_equal, order=30)
    for ind in range(0, min(len(indofmax[0]), len(indofmin[0]))):
        if indofmax[0][0] < indofmin[0][0]:
            wing.at[indofmax[0][ind]:indofmin[0][ind], 'FW_Bck_strk'] = 0
            wing.at[indofmin[0][ind]:indofmax[0][ind + 1], 'FW_Bck_strk'] = 1
           
        else:
            wing.at[indofmin[0][ind]:indofmax[0][ind], 'FW_Bck_strk'] = 0
            wing.at[indofmax[0][ind]:indofmin[0][ind + 1], 'FW_Bck_strk'] = 1
        wing.at[indofmax[0][ind]:indofmax[0][ind + 1], 'FullStrk'] = ind

    if plot == 1:
        fig = plt.figure()
        plt.plot(wing['FW_Bck_strk'] * 180)
        plt.plot(wing.index, wing['phi' + wn], '-*')
        plt.plot(wing.index[indofmax], wing['phi' + wn].iloc[indofmax], 'b*')
        plt.plot(wing.index[indofmin], wing['phi' + wn].iloc[indofmin], 'r*')
        plt.grid()
        plt.xlabel('index')
        plt.ylabel('phi [deg]')
        fig.legend(['0 - FWD ,1 - BCK','phi','maximum','minimum'])

    return wing,indofmin,indofmax


def find_FBstrk(wing,plot = 0):
    # find forward and backward stroke. update in data frame and return indices
    if 'rw' in str(wing.keys()):
        wng = 'rw'
    else:
        wng = 'lw'
    phival = wing['phi_'+ wng + ' smth'].values
    indofmax = argrelextrema(phival, np.greater, order=50,mode = 'wrap')
    indofmin = argrelextrema(phival, np.less, order=50,mode = 'wrap')

    inds = [wing.index[0],wing.index[-1]]
    FB_strk = np.zeros((inds[1] - inds[0] + 1,9))

    ini_inds = [max(min(indofmax[0]),min(indofmin[0]))]
    if (ini_inds[0] == indofmin).any():
        ini_inds = indofmax[0][1]
    
    indofmax = indofmax[0][indofmax[0] >= ini_inds]
    indofmin = indofmin[0][indofmin[0] >= ini_inds]

    cnt = 5
    for ind in range(0, min(len(indofmax), len(indofmin)) - 1):
# =============================================================================
#         if indofmax[0] < indofmin[0]:
# =============================================================================
        frphimx = wing['frames'][indofmax[ind] + wing.index[0] - 1]
        frphimin = wing['frames'][indofmin[ind] + wing.index[0] - 1]
        tmphimx = wing['time [ms]'][indofmax[ind] + wing.index[0] - 1]
        tmphimin = wing['time [ms]'][indofmin[ind] + wing.index[0] - 1]
        
        
        FB_strk[indofmax[ind]:indofmin[ind]] = [cnt,0,phival[indofmax[ind]],phival[indofmin[ind]],cnt,
                                                frphimx,frphimin,tmphimx,tmphimin]
        FB_strk[indofmin[ind]:indofmax[ind + 1]] = [cnt + 1,1,phival[indofmin[ind]],phival[indofmax[ind]],cnt,
                                                    frphimx,frphimin,tmphimx,tmphimin]

# =============================================================================
#         else:
#             FB_strk[indofmin[ind]:indofmax[ind]] = [cnt,0,phival[indofmin[ind]],phival[indofmax[ind]]]
#             FB_strk[indofmax[ind]:indofmin[ind + 1]] = [cnt + 1,1,phival[indofmax[ind]],phival[indofmin[ind]]]
# =============================================================================
        cnt = cnt + 2

    if plot == 1:
        fig = plt.figure()
        plt.plot(range(0,len(FB_strk)),FB_strk[:,1] * 180)
        plt.plot( wing['phi_'+ wng + ' smth'].values, '-*')
        plt.plot( indofmax,wing['phi_'+ wng + ' smth'].values[indofmax], 'b*')
        plt.plot(indofmin,wing['phi_'+ wng + ' smth'].values[indofmin], 'r*')
        plt.grid()
        plt.xlabel('index')
        plt.ylabel('phi [deg]')
        fig.legend(['0 - FWD ,1 - BCK','phi','maximum','minimum'])
    return FB_strk,inds,indofmin,indofmax

def addStrks2DF(rw,debugFBstrk = 0):
    strks = rw.groupby('mov').apply(lambda x: find_FBstrk(x, plot=0))
    if len(strks) == 0:
        return 0,0
    if debugFBstrk == 1:
        find_FBstrk(rw, plot=1)

    addst = []
    mrkFB = []
    addinds_max = []
    addinds_min = []
    fullST = []
    maxfr = []
    minfr = []
    maxtm = []
    mintm = []
    addinds_max = []
    addinds_min = []
    if 'rw' in str(rw.keys()):
        wng = 'rw'
    else:
        wng = 'lw'
    for nm in rw['mov'].unique()[::-1]:
        addst = np.append(strks[nm][0][:,0],addst)
        mrkFB = np.append(strks[nm][0][:,1],mrkFB)
        fullST = np.append(strks[nm][0][:,4],fullST)
        maxfr = np.append(strks[nm][0][:,5],maxfr)
        minfr = np.append(strks[nm][0][:,6],minfr)
        maxtm = np.append(strks[nm][0][:,7],maxtm)
        mintm = np.append(strks[nm][0][:,8],mintm)
        addinds_max = np.append(strks[nm][0][:,2],addinds_max)
        addinds_min = np.append(strks[nm][0][:,3],addinds_min)
    rw['strks_' + wng] = addst
    rw['BCK FWD_' + wng] = mrkFB
    rw['Full_strk' + wng] = fullST
    rw['maxmin_val1' + wng] = addinds_max
    rw['maxmin_val2' + wng] = addinds_min
    rw['maxphi_fr' + wng] = maxfr
    rw['minphi_fr' + wng] = minfr
    rw['maxphi_tm' + wng] = maxtm
    rw['minphi_tm' + wng] = mintm
    rw.loc[rw['strks_' + wng] == 0,'strks_' + wng] = None
    return rw,strks


def lstsqr(x,y):
    missing = np.isnan(x)
    ynew = y[~missing]
    xnew = x[~missing]

    A1 = np.vstack([xnew, np.ones(len(xnew))]).T

    m1, c1 = np.linalg.lstsq(A1, ynew, rcond=None)[0]
    return m1,c1,ynew,xnew



def pqr_pqr_dot(df):
    # calculate the body angular acceleratrion and velocity: pqr and pqr_dot
    psi=df['roll smth'];theta=df['pitch smth'];phi = df['yaw smth']
    psi_dot=df['roll smth_dot'];theta_dot = df['pitch smth_dot'];phi_dot = df['yaw smth_dot'];
    psi_dot_dot=df['roll smth_dot_dot'];theta_dot_dot = df['pitch smth_dot_dot'];phi_dot_dot = df['yaw smth_dot_dot'];
    
    df['p'] = psi_dot-phi_dot*np.sin(np.radians(theta));
    df['q'] = theta_dot*np.cos(np.radians(psi))+phi_dot*np.sin(np.radians(psi))*np.cos(np.radians(theta));
    df['r'] = -theta_dot*np.sin(np.radians(psi))+phi_dot*np.cos(np.radians(psi))*np.cos(np.radians(theta));
    
    
    df['p_dot'] = psi_dot_dot-phi_dot_dot*np.sin(np.radians(theta))-phi_dot*theta_dot*np.cos(np.radians(theta));
    df['q_dot'] = (theta_dot_dot*np.cos(np.radians(psi))-theta_dot*np.sin(np.radians(psi))*psi_dot
                    +phi_dot_dot*np.sin(np.radians(psi)*np.cos(theta)+
                    +phi_dot*np.cos(np.radians(psi))*np.cos(np.radians(theta))*psi_dot
                    -phi_dot*np.sin(np.radians(psi))*np.sin(np.radians(theta))*theta_dot))
                
                
    df['r_dot'] = (-theta_dot_dot*np.sin(np.radians(psi))-theta_dot*np.cos(np.radians(psi))*psi_dot
                    +phi_dot_dot*np.cos(np.radians(psi))*np.cos(np.radians(theta))
                    -phi_dot*np.sin(np.radians(psi))*np.cos(np.radians(theta))*psi_dot
                    -phi_dot*np.cos(np.radians(psi))*np.sin(np.radians(theta))*theta_dot)
    return df



def calculateSP_V(angdict,body_vectors):
    ind = 0
    dt = np.diff(angdict['body']['time [ms]'])[0] / 1000

    for mvnm in np.unique(angdict['body']['mov']):
        mvdf = angdict['body'].loc[angdict['body']['mov'] == mvnm]
        velocity_deltaCM = mvdf[['X smth_dot', 'Y smth_dot', 'Z smth_dot']].to_numpy()
        VzSP = np.sum(velocity_deltaCM * body_vectors[mvnm]['strkPlan'], axis=1).T/dt

        XonSP = body_vectors[mvnm]['X'] - body_vectors[mvnm]['strkPlan']
        YonSP = body_vectors[mvnm]['Y'] - body_vectors[mvnm]['strkPlan']
        XonSP = XonSP / np.sqrt(np.sum(XonSP ** 2, axis=1))[:, np.newaxis]
        YonSP = YonSP / np.sqrt(np.sum(YonSP ** 2, axis=1))[:, np.newaxis]

        VxSP = np.sum(XonSP * velocity_deltaCM, axis=1).T/dt
        VySP = np.sum(YonSP * velocity_deltaCM, axis=1).T/dt
        if ind == 0:
            Vx = VxSP
            Vy = VxSP
            Vz = VxSP
        else:
            Vx = np.hstack((Vx,VxSP))
            Vy = np.hstack((Vy, VySP))
            Vz = np.hstack((Vz, VzSP))
        ind += 1
    angdict['body'][['Vxsp','Vysp','Vzsp']] = pd.DataFrame(np.hstack((Vx[:, np.newaxis],Vy[:, np.newaxis],Vz[:, np.newaxis])))
    angdict['body'][['X smth_dot','Y smth_dot','Z smth_dot']] = angdict['body'][['X smth_dot', 'Y smth_dot', 'Z smth_dot']].to_numpy()/dt
    angdict['body'][['X smth_dot_dot','Y smth_dot_dot','Z smth_dot_dot']] = angdict['body'][['X smth_dot_dot', 'Y smth_dot_dot', 'Z smth_dot_dot']].to_numpy()/dt/dt


    return angdict


