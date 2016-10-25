#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 11:25:38 2016

@author: Sivanand Achanta
"""

import numpy as np


def ola(y, w, frSize, frShift, tsmf):

    len_synsig = int(np.floor(tsmf*len(y[0])))

    ati = np.floor(np.arange(0, len(y[0]) - frSize, frShift)).astype('int')  # analysis time instants
    sti = np.floor(np.arange(0, len_synsig - frSize, int(tsmf*frShift))).astype('int')  # synthesis time instants

    if len(sti) > len(ati):
        sti = sti[:-1]

    v = w**2

    ys_ola = np.zeros((len_synsig, 1))
    norm_sig = np.zeros((len_synsig, 1)) + 1e-6

    for ti in range(len(sti)):

        sp = int(sti[ti])
        ep = int(sti[ti] + frSize)

        x_w = w*y[0][int(ati[ti]):int(ati[ti] + frSize)]

        # OLA signal
        ys_ola[sp:ep] = ys_ola[sp:ep] + ((w.reshape(1, len(w)))*x_w).T

        # Normalization signal
        norm_sig[sp:ep] = norm_sig[sp:ep] + v.reshape(len(v), 1)

    ysn_ola = ys_ola/norm_sig

    return(ysn_ola)


def sola(y, w, frSize, frShift, tsmf, max_delay):

    len_synsig = int(np.ceil(tsmf*len(y[0])))
    ati = np.floor(np.arange(0, len(y[0]) - frSize, frShift)).astype('int')  # analysis time instants
    sti = np.floor(np.arange(0, len_synsig - frSize, tsmf*frShift)).astype('int')  # synthesis time instants

    if len(sti) > len(ati):
        sti = sti[:-1]

    v = w**2
    delay = 0

    ys_sola = np.zeros((len_synsig, 1))
    norm_sig = np.zeros((len_synsig, 1)) + 1e-6

    for ti in range(len(sti)):
        sp_nodelay = int(sti[ti])
        ep_nodelay = int(sti[ti] + frSize)

        x_w = w*y[0][int(ati[ti]):int(ati[ti] + frSize)]

        if ti > 0:
            r1 = np.correlate(x_w.reshape(frSize,), ys_sola[sp_nodelay:ep_nodelay].reshape(frSize,), 'full')
            r2 = np.correlate(ys_sola[sp_nodelay:ep_nodelay].reshape(frSize,), x_w.reshape(frSize,), 'full')

            r1_max = np.max(r1[frSize:frSize+max_delay])
            r2_max = np.max(r2[frSize:frSize+max_delay])

            d1 = np.argmax(r1[frSize:frSize+max_delay])
            d2 = np.argmax(r2[frSize:frSize+max_delay])

            if r1_max > r2_max:
                delay = -d1
            else:
                delay = d2

        sp = int(sti[ti] + delay)
        ep = int(sti[ti] + frSize + delay)

        # OLA signal
        ys_sola[sp:ep] = ys_sola[sp:ep] + ((w.reshape(1, len(w)))*x_w).T

        # Normalization signal
        norm_sig[sp:ep] = norm_sig[sp:ep] + v.reshape(len(v), 1)

    ysn_sola = ys_sola/norm_sig

    return(ysn_sola)


'''
# Debugging OLA methods
fs = 16000
frSizems = 30
frShiftms = 10
frSize = int((frSizems/1000)*fs)
frShift = int((frShiftms/1000)*fs)
frOvlap = frSize - frShift
nfft = 512
nfftby2 = np.round(nfft/2 + 1)
hfhz = np.linspace(0, fs/2, nfftby2)
w = np.hamming(frSize)
max_delay = 80  # SOLA parmaeter

# Plot flags
pflag = 1
plot_crosscorr = 0

# Matlab plot parameters
font_size = 14

# Time scale modification factor
tsmf = 1

# Generate a synthetic signal (pulse train or a rectangular wave)
t_on = 10
t_off = 160 - t_on

pulse = np.concatenate((np.ones((1, t_on)), np.zeros((1, t_off))), axis=1)
y = np.concatenate((np.zeros((1, 40)), np.tile(pulse, 16)), axis=1)


ysn_ola = ola(y, w, frSize, frShift, tsmf)
ysn_sola = sola(y, w, frSize, frShift, tsmf, max_delay)

if pflag:
    import matplotlib.pyplot as plt

    plt.subplot(311)
    plt.plot(y[0])
    plt.xlim(0, len(y[0]))

    plt.subplot(312)
    plt.plot(ysn_ola, 'r-')
    plt.xlim(0, len(ysn_ola))

    plt.subplot(313)
    plt.plot(ysn_sola, 'm-')
    plt.xlim(0, len(ysn_sola))

'''
