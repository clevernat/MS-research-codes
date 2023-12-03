import numpy as np
import sys
sys.path.append('../')
from functions_calc import geo_idx

def aeri_reg(aeri):
    aeri_time_reg = np.arange(aeri['time'][0],aeri['time'][-1]+1,300)
    aeri_t_reg = np.full((len(aeri_time_reg),len(aeri['h'])),np.nan)
    aeri_rh_reg = np.full((len(aeri_time_reg),len(aeri['h'])),np.nan)
    aeri_time_good = []
    aeri_t_good = []
    aeri_rh_good = []
    for j in range(0, len(aeri['time'])):
        good1 = (aeri['gamma'][j] <= 3 and aeri['rmsr'][j] < 6)
        if good1 == True:
            aeri_time_good.append(aeri['time'][j])
            aeri_t_good.append(aeri['t'][j])
            aeri_rh_good.append(aeri['RH'][j])

    for j in range(0, len(aeri_time_good)):
        aeri_time_idx = geo_idx.geo_idx(aeri_time_good[j], aeri_time_reg)
        aeri_t_reg[aeri_time_idx] = aeri_t_good[j]
        aeri_rh_reg[aeri_time_idx] = aeri_rh_good[j]
        if j < len(aeri_time_good)-1 and j > 0:
            if aeri_time_good[j]-aeri_time_good[j-1] <= 1500:
                aeri_rh_reg[aeri_time_idx+1] = aeri_rh_good[j]
            if aeri_time_good[j+1]-aeri_time_good[j] == 600:
                aeri_t_reg[aeri_time_idx+1] = np.mean([aeri_t_good[j],aeri_t_good[j+1]],axis=0)
                aeri_rh_reg[aeri_time_idx+1] = np.mean([aeri_rh_good[j],aeri_rh_good[j+1]],axis=0)
            if aeri_time_good[j+1]-aeri_time_good[j] == 900:
                aeri_t_reg[aeri_time_idx+1] = np.mean([aeri_t_good[j],aeri_t_good[j],aeri_t_good[j+1]],axis=0)
                aeri_rh_reg[aeri_time_idx+1] = np.mean([aeri_rh_good[j],aeri_rh_good[j],aeri_rh_good[j+1]],axis=0)
                aeri_t_reg[aeri_time_idx+2] = np.mean([aeri_t_good[j],aeri_t_good[j+1],aeri_t_good[j+1]],axis=0)
                aeri_rh_reg[aeri_time_idx+2] = np.mean([aeri_rh_good[j],aeri_rh_good[j+1],aeri_rh_good[j+1]],axis=0)
            if aeri_time_good[j+1]-aeri_time_good[j] == 1200:
                aeri_t_reg[aeri_time_idx+1] = np.mean([aeri_t_good[j],aeri_t_good[j],aeri_t_good[j],aeri_t_good[j+1]],axis=0)
                aeri_rh_reg[aeri_time_idx+1] = np.mean([aeri_rh_good[j],aeri_rh_good[j],aeri_rh_good[j],aeri_rh_good[j+1]],axis=0)
                aeri_t_reg[aeri_time_idx+2] = np.mean([aeri_t_good[j],aeri_t_good[j+1]],axis=0)
                aeri_rh_reg[aeri_time_idx+2] = np.mean([aeri_rh_good[j],aeri_rh_good[j+1]],axis=0)
                aeri_t_reg[aeri_time_idx+3] = np.mean([aeri_t_good[j],aeri_t_good[j+1],aeri_t_good[j+1],aeri_t_good[j+1]],axis=0)
                aeri_rh_reg[aeri_time_idx+3] = np.mean([aeri_rh_good[j],aeri_rh_good[j+1],aeri_rh_good[j+1],aeri_rh_good[j+1]],axis=0)
            if aeri_time_good[j+1]-aeri_time_good[j] == 1500:
                aeri_t_reg[aeri_time_idx+1] = np.mean([aeri_t_good[j],aeri_t_good[j],aeri_t_good[j],aeri_t_good[j],aeri_t_good[j+1]],axis=0)
                aeri_rh_reg[aeri_time_idx+1] = np.mean([aeri_rh_good[j],aeri_rh_good[j],aeri_rh_good[j],aeri_rh_good[j],aeri_rh_good[j+1]],axis=0)
                aeri_t_reg[aeri_time_idx+2] = np.mean([aeri_t_good[j],aeri_t_good[j],aeri_t_good[j],aeri_t_good[j+1],aeri_t_good[j+1]],axis=0)
                aeri_rh_reg[aeri_time_idx+2] = np.mean([aeri_rh_good[j],aeri_rh_good[j],aeri_rh_good[j],aeri_rh_good[j+1],aeri_rh_good[j+1]],axis=0)
                aeri_t_reg[aeri_time_idx+3] = np.mean([aeri_t_good[j],aeri_t_good[j],aeri_t_good[j+1],aeri_t_good[j+1],aeri_t_good[j+1]],axis=0)
                aeri_rh_reg[aeri_time_idx+3] = np.mean([aeri_rh_good[j],aeri_rh_good[j],aeri_rh_good[j+1],aeri_rh_good[j+1],aeri_rh_good[j+1]],axis=0)
                aeri_t_reg[aeri_time_idx+4] = np.mean([aeri_t_good[j],aeri_t_good[j+1],aeri_t_good[j+1],aeri_t_good[j+1],aeri_t_good[j+1]],axis=0)
                aeri_rh_reg[aeri_time_idx+4] = np.mean([aeri_rh_good[j],aeri_rh_good[j+1],aeri_rh_good[j+1],aeri_rh_good[j+1],aeri_rh_good[j+1]],axis=0)
    return aeri_t_reg,aeri_rh_reg,aeri_time_reg
