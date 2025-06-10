import pickle
import scipy as sp
import numpy as np
import cmath as cm
import scipy.stats as stats 
import scipy.signal as signal
import scipy.fftpack as fftpack
import pylab as pl
import dataanalysis as da

def arangevar(traj_collection,var,nperseg,step):
# split trajectories in windows for spectral analysis and arange all windows in one array
# nperseg: window length for spectral analysis
# step: determine the overlag (if same as nperseg no overlap will be made)
    wl=[]
    for tjn, traj in traj_collection.iteritems():
        w_traj = sp.array(traj[var])
        #print tjn, w_traj.shape
        # to do: make sure all data is evenly spaced 
        if tjn == 96: # this buoy has 15min time steps 
            continue#w_traj=w_traj[::2]
        for i in range(0,w_traj.shape[0]-nperseg,step):
#        w=w+list(w_traj[i:nperseg])
            wl.append(w_traj[i:i+nperseg])

    w = sp.zeros((len(wl),nperseg))
    for i in range(len(wl)):
        w[i] = wl[i]
    w[sp.isfinite(w)==False]=0
    return w


def rotaryspectra(x, dt, rescale=1.):
    #window = signal.blackman(nperseg)    
    nperseg = x.shape[1]
    M = x.shape[0]
    N = M*nperseg
    Yw = fftpack.fft(x, n=nperseg) * dt * rescale # discrete fourier transform, rescale fft by constant if a window was applied.
    #Gww = 2/(N*dt) * (Yw.conj()*Yw).sum(axis=0) # spectral density function as defined in Emery and Thomson 1997, pp. 427 in 
    Gww = 1/(N*dt) * (Yw.conj()*Yw).sum(axis=0) # same, but for two-sided spectra
    #Gww = 1. / dt * (Yw.conj()*Yw).mean(axis=0) 
    f = np.fft.fftfreq(nperseg,dt) # equiv. to f = sp.arange(nperseg/2.)/(nperseg * dt) but with shifted frequencies

    Scc = Gww[:nperseg/2] #counterclockwise spectrum
    Scl = Gww[nperseg/2:] #clockwise spectrum
    fpos = f[:nperseg/2]
    fneg = -f[nperseg/2:]

    return Scc, Scl, fpos, fneg


def rotary_crossspectra(x1, x2, dt, rescale=1., averagesegments=True ):
    #window = signal.blackman(nperseg)    
    nperseg = x1.shape[1]
    M = x1.shape[0]
    N = M*nperseg
    Y1 = fftpack.fft(x1, n=nperseg) * dt * rescale
    Y2 = fftpack.fft(x2, n=nperseg) * dt * rescale
    f = np.fft.fftfreq(nperseg,dt) #f = sp.arange(nperseg/2.)/(nperseg * dt)

    if averagesegments == True:
        S = 1/(N*dt) * sp.sum(Y1 * Y2.conj(), axis=0)
    else:
        S = Y1 * Y2.conj() / (nperseg*dt)

    S = fftpack.fftshift(S)
    f = fftpack.fftshift(f)

    return S, f


def demean(x,axis=1):
    average = x.mean(axis=axis)
    average = np.tile(average,(x.shape[axis],axis)).transpose()
    return x-average

detrend = demean

def applywindow(x,window):
    W = np.tile(window, (x.shape[0], 1))
    return x*W

def pad(x,n):
    z=sp.zeros((x.shape[0],n))
    return np.concatenate((z,x,z), axis=1)


