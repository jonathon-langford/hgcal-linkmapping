#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os

def plot(variable,savename="hist.png",binwidth=1,xtitle='Number of words on a single lpGBT',outdir="."):

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    fig = plt.figure()
    binwidth=binwidth
    plt.hist(variable, bins=np.arange(min(variable), max(variable) + binwidth, binwidth))
    plt.ylabel('Number of Entries')
    plt.xlabel(xtitle)
    plt.savefig(os.path.join(outdir, savename))
    

def plot2D(variable_x,variable_y,savename="hist2D.png",binwidthx=1,binwidthy=1,xtitle='Number of words on a single lpGBT',outdir="."):

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    fig = plt.figure()
    binwidthx=binwidthx
    binwidthy=binwidthy
    plt.hist2d(variable_x,variable_y,bins=[np.arange(min(variable_x), max(variable_x) + binwidthx, binwidthx),np.arange(min(variable_y), max(variable_y) + binwidthy, binwidthy)])
#    plt.hist2d(variable_x,variable_y,bins=[np.arange(0.9, max(variable_x) + binwidthx, binwidthx),np.arange(min(variable_y), max(variable_y) + binwidthy, binwidthy)])
    plt.colorbar()
    plt.ylabel('Layer')
    plt.xlabel(xtitle)
    plt.savefig(os.path.join(outdir, savename))

