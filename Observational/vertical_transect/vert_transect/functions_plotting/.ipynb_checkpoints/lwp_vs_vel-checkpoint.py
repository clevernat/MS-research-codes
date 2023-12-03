import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from matplotlib.transforms import Bbox
from matplotlib.lines import Line2D

def lwp_vs_vel(lwp, vel, id, st):
    fig, ax = plt.subplots(1,1, figsize=(5,5))
    ax.scatter(lwp,vel, c = 'k', s = 1)
    os.makedirs("lwp_vel", exist_ok = True)
    plt.savefig("lwp_vel/lwp_vel_"+st+"_"+id+".png", dpi = 300)
    plt.close()
