import matplotlib
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
import os
from pathlib import Path

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Some general Settings
SMALL_SIZE = 10
MEDIUM_SIZE = 11
BIGGER_SIZE = 12
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

myred1 = (175/255, 22/255, 31/255)
myred2 = (224/255, 49/255, 47/255)
myblue1 = (65/255, 67/255, 158/255)

RES_ID = "fullrun2/all_constaints"


res_path = Path.cwd()

# Load Data
r_series = pd.read_csv(res_path.joinpath("rand_max200_min30_n10000.csv"), float_precision='round_trip', header=2)
nyc_series = pd.read_csv(res_path.joinpath("nyiso_nyc_2018-10-08_to_2018-10-10.csv"), float_precision='round_trip',header=2)

fig = plt.figure()
fig.set_size_inches(6,3)
gs = GridSpec(1,2)

ax1 = plt.subplot(gs[0,0])
tr_series = r_series['price'].values[0:500]
ax1.plot(tr_series)

ax2 = plt.subplot(gs[0,1])
tr_series = nyc_series['price'].values[0:500]
ax2.plot(tr_series)
ax2.set_ylim(-10,60)

for ax in [ax1,ax2]:
    ax.set_ylabel("Price at substation in $/MWh")
    ax.set_xlabel("Timestep")
    ax.grid('on')
    ax.set_xticks([0,100,200,300,400,500])
    
plt.tight_layout()


