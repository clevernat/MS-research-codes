import numpy as np
import metpy as mp
import netCDF4
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import itertools
import copy
import metpy.calc as mpcalc
from metpy.units import units
from datetime import datetime
from scipy import ndimage
### functions from other files that are required
from functions_plotting import cfad, timeseries,interp_plot,lwp_vs_vel,cfad_fillplot
from functions_calc import geo_idx, date_finder, aeri_reg
from functions_dics import vert_transect_dics, predefined
from functools import partial
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

print("Use predefined or interactive environment? Enter p for predefined or i for interactive:")
pre_or_int = str(input())
while pre_or_int != "i" and pre_or_int != "p":
    print("Enter correct letter name: p or i")
    pre_or_int = str(input())

### Interactive environment ###
if pre_or_int == "p":
    data_set, panel_name, cao_times, folder_save = predefined.predefined()

### Interactive environment ###
if pre_or_int == "i":
    print("----------------------------------------------------------------------------------------------")
    print("Enter ge to use general mode or arscl to use processed radar data:")
    data_set = str(input())
    while data_set != "ge" and data_set != "arscl":
        print("Enter correct dataset name: ge or arscl")
        data_set = str(input())
    print("----------------------------------------------------------------------------------------------")
    print("What to call the folder where figures are saved:")
    folder_save = str(input())

    cao_times = [[0,0,0,0,0]]
    cao_times_idx = ["start day between 0 and 182", "start day minute between 0 and 1439", "end day between 0 and 182", "end day minute between 0 and 1439"]
    for i in range(0,4):
        print("----------------------------------------------------------------------------------------------")
        print("Enter "+ cao_times_idx[i])
        cao_times[0][i] = int(input())
    cao_times[0][4] = int((cao_times[0][2]-cao_times[0][0])*1440+(cao_times[0][3]-cao_times[0][1]))
    print(cao_times)
    ### Plotting Positions ###
    panel_name = vert_transect_dics.panel_name_dic()
    for var in panel_name.keys():
        print("----------------------------------------------------------------------------------------------")
        print("Plot "+var+" at which position? Enter -1 if variable is not supposed to be plotted.")
        panel_name[var] = int(input())
        if panel_name[var] >= 1:
            panel_name[var] -= 1

### define basic variables of datasets ###
if data_set == "ge":
    folder, nc_file, ref_name, h_name, noise_name  = "anxkazrcfrgeM1", "anxkazrcfrgeM1.a1.20", "reflectivity", "range", "signal_to_noise_ratio_copolar_h"
if data_set == "arscl":
    folder, nc_file, ref_name, h_name, noise_name  = "anxarsclkazr1kolliasM1", "anxarsclkazr1kolliasM1.c1.20", "reflectivity_best_estimate", "height", "signal_to_noise_ratio"
### Function to take required files ###
def fi(files, s, e):
    files = files[s:e]
    return files

### Function to flatten arrays ###
def f2d(x):
    y = np.asarray(list(itertools.chain(*x)))
    return y


#path = "/netdata/virga/lackner_comble/"
path = "/glade/scratch/noteng/ftp.archive.arm.gov/otengn1/240475/"
### basic directions ###
shift_0, len_cases = 0, 48
cfad_plot = False

### main file paths with dictornaries ###
files, files_path = vert_transect_dics.main_dics(folder, nc_file)

#for i in range(103, 105):
for i in range(0, 183):
    day, month = date_finder.date_finder(i)
    for var in files.keys():
        files[var].append(sorted(glob.glob(path+files_path[var]+month+day+"*")))
### Load data ###
ranges_all,reflectivity_all,velocity_all,sw_all = [], [], [], []
#aeri_days = [8,65,102,103,104,110,117,118,119,120]
for i in range(shift_0, len(cao_times)):
    print("Running case "+str(i))

    ### Base Time of case ###
    case = cao_times[i]
    base_time_start = (365*38+12*366-31)*86400 # Base time on 1 Dec 2019 0:00:00
    bts, bte = case[0]*86400+case[1]*60+base_time_start, case[2]*86400+case[3]*60+base_time_start # Base time Start and endof case
    #start_day, end_day, start_min, end_min = case[0]-103, case[2]-103, case[1], case[3]
    start_day, end_day, start_min, end_min = case[0], case[2], case[1], case[3]

    #ldpr, backscatter, ranges_mpl, time_mpl = [], [], [], []
    ldpr, backscatter = [], []
    u_vwp, v_vwp, ranges_vwp, u_vwp1, v_vwp1, ranges_vwp1 = [], [], [], [], [], []


    vars, files_case = vert_transect_dics.vars_dics()
   

    for var in files_case.keys():
        #if var == "mpl" and case[0] >= 73:
            #files_case[var] = f2d(fi(files[var], start_day, end_day+1))
        if var == "cld" and case[0] >= 73:
            files_case[var] = f2d(fi(files[var], start_day, end_day+1))
        elif var == "vwp" and case[0] >= 91:
            files_case[var] = f2d(fi(files[var], start_day, end_day+1))
       # elif var == "aeri" and case[0] in aeri_days and case[2] in aeri_days:
        #    files_case[var] = f2d(fi(files[var], start_day, end_day+1))
        else:
            files_case[var] = f2d(fi(files[var], start_day, end_day+1))
    print(files_case)

    for j in range(0,len(files_case['kazr'])):
        f = netCDF4.Dataset(files_case['kazr'][j])
        base_time_kazr = np.array(f.variables['base_time'][:])
        vars['kazr']['time'].append(np.array(f.variables['time'][:])+base_time_kazr)
        vars['kazr']['ref'].append(np.array(f.variables[ref_name][:]))
        vars['kazr']['h'].append(np.array(f.variables[h_name][:]))
        vars['kazr']['noise'].append(np.array(f.variables[noise_name][:]))
        vars['kazr']['vel'].append(np.array(f.variables['mean_doppler_velocity'][:]))
        vars['kazr']['sw'].append(np.array(f.variables['spectral_width'][:]))
        f.close()

    for j in range(0,len(files_case['micro'])):
        f = netCDF4.Dataset(files_case['micro'][j])
        base_time_micro = np.array(f.variables['base_time'][:])
        vars['micro']['time'].append(np.array(f.variables['time'][:])+base_time_micro)
        vars['micro']['lwc'].append(np.array(f.variables['liquid_water_content'][:]))
        vars['micro']['lwc_u'].append(np.array(f.variables['liquid_water_content_uncertainty_random'][:]))
        vars['micro']['iwc'].append(np.array(f.variables['ice_water_content'][:]))
        vars['micro']['iwc_u'].append(np.array(f.variables['ice_water_content_uncertainty_random'][:]))
        vars['micro']['h'].append(np.array(f.variables['height'][:]))
        f.close()

    for j in range(0, len(files_case['vpt'])):
        f = netCDF4.Dataset(files_case['vpt'][j])
        base_time_vpt = np.array(f.variables['base_time'][:])
        vars['vpt']['time'].append(np.array(f.variables['time'][:])+base_time_vpt)
        vars['vpt']['ref'].append(np.array(f.variables['reflectivity'][:]))
        vars['vpt']['h'].append(np.array(f.variables["range"][:]))
        vars['vpt']['noise'].append(np.array(f.variables['signal_to_noise_ratio_copolar_h'][:]))
        vars['vpt']['vel'].append(np.array(f.variables['mean_doppler_velocity'][:]))
        vars['vpt']['sw'].append(np.array(f.variables['spectral_width'][:]))
        for var in ['ref','noise','vel','sw']:
            vars['vpt'][var][j][-1,:] = np.nan
        f.close()

    for j in range(0, len(files_case['mwr'])):
        f = netCDF4.Dataset(files_case['mwr'][j])
        base_time_mwr = np.array(f.variables['base_time'][:])
        vars['mwr']['time'].append(np.array(f.variables['time'][:])+base_time_mwr)
        vars['mwr']['lwp'].append(np.array(f.variables['be_lwp'][:]))
        vars['mwr']['pwv'].append(np.array(f.variables['be_pwv'][:]))
        f.close()

        f = netCDF4.Dataset(files_case['cloud'][j])
        base_time_cloud = np.array(f.variables['base_time'][:])
        vars['cloud']['time'].append(np.array(f.variables['time'][:])+base_time_cloud)
        vars['cloud']['fcth'].append(np.array(f.variables['radar_first_top'][:]))
        vars['cloud']['cb'].append(np.array(f.variables['cloud_base_best_estimate'][:]))
        f.close()
        cdt, cdb = copy.deepcopy(vars['cloud']['fcth'][j]), copy.deepcopy(vars['cloud']['cb'][j])
        invalidt, invalidb = np.logical_or(cdt < 0, cdt > 15000), np.logical_or(cdb < 0, cdb > 15000)
        cdt[invalidt], cdb[invalidb] = np.nan, np.nan
        vars['cloud']['cd'].append(cdt - cdb)

        #f = netCDF4.Dataset(files_case['iwp'][j])
        #vars['iwp']['iwp'].append(np.array(f.variables['iwp'][:]))
        #f.close()

        f = netCDF4.Dataset(files_case['is'][j])
        base_time_is = np.array(f.variables['base_time'][:])
        vars['is']['time'].append(np.array(f.variables['time'][:])+base_time_is)
        vars['is']['t'].append(np.array(f.variables['temp'][:]))
        vars['is']['dp'].append(np.array(f.variables['dp'][:]))
        vars['is']['pres'].append(np.array(f.variables['bar_pres'][:]*10))
        vars['is']['wspd'].append(np.array(f.variables['wspd'][:]))
        vars['is']['h'].append(np.array(f.variables['height'][:]*1000))
        vars['is']['u'].append(np.array(f.variables['u_wind'][:]))
        vars['is']['v'].append(np.array(f.variables['v_wind'][:]))
        rh_scaled = np.array(f.variables['rh_scaled'][:])/100.0
        f.close()
        dp_scaled = np.asarray(mpcalc.dewpoint_from_relative_humidity(vars['is']['t'][-1]*units.degC,rh_scaled))
        vars['is']['ept'].append(np.asarray(mpcalc.equivalent_potential_temperature(vars['is']['pres'][-1]*units.hPa, vars['is']['t'][-1]*units.degC, dp_scaled*units.degC)))

        #if case[0] in aeri_days and case[2] in aeri_days:
        #    f = netCDF4.Dataset(files_case['aeri'][j])
        #    vars['aeri']['time'].append(np.array(f.variables['time'][:])+base_time_is)
        #    vars['aeri']['t'].append(np.array(f.variables['temperature'][:]))
        #    vars['aeri']['RH'].append(np.array(f.variables['rh'][:]))
        #    vars['aeri']['gamma'].append(np.array(f.variables['gamma'][:]))
        #    vars['aeri']['rmsr'].append(np.array(f.variables['rmsr'][:]))
        #    vars['aeri']['h'].append(np.array(f.variables['height'][:])*1000)
        #    f.close()

        #if case[0] >= 73:
            #f = netCDF4.Dataset(files_case['mpl'][j])
            #ldpr.append(np.array(f.variables['linear_depol_ratio'][:]))
            #backscatter.append(np.array(f.variables['backscatter'][:]))
            #time_mpl.append(np.array(f.variables['time'][:])+base_time_cloud)
            #ranges_mpl.append(np.array(f.variables['height'][:]*1000))
            #f.close()

        f = netCDF4.Dataset(files_case['cld'][j])
        base_time_phase = np.array(f.variables['base_time'][:])
        vars['cld']['time'].append(np.array(f.variables['time'][:])+base_time_phase)
        vars['cld']['phase'].append(np.array(f.variables['cloud_phase_mplgr'][:]))
        vars['cld']['h'].append(np.array(f.variables['height'][:])*1000)
        f.close()

        if case[0] >= 91:
            f = netCDF4.Dataset(files_case['vwp'][j])
            ranges_vwp.append(np.array(f.variables['height_p'][:])[:,0])
            u_vwp.append(np.array(f.variables['u_wind'][:])[:,:,0])
            v_vwp.append(np.array(f.variables['v_wind'][:])[:,:,0])
            ranges_vwp1.append(np.array(f.variables['height_p'][:])[:,1])
            u_vwp1.append(np.array(f.variables['u_wind'][:])[:,:,1])
            v_vwp1.append(np.array(f.variables['v_wind'][:])[:,:,1])
            f.close()

    for j in range(0, len(files_case['sur'])):
        f = netCDF4.Dataset(files_case['sur'][j])
        base_time_sur = np.array(f.variables['base_time'][:])
        vars['sur']['time'].append(np.array(f.variables['time'][:])+base_time_sur)
        vars['sur']['t'].append(np.array(f.variables['temperature'][:]))
        vars['sur']['wspd'].append(np.array(f.variables['wind_speed'][:]))
        vars['sur']['pres'].append(np.array(f.variables['pressure'][:]))
        vars['sur']['prec'].append(np.array(f.variables['precipitation_pwd'][:]))
        vars['sur']['dirc'].append(np.array(f.variables['wind_direction'][:]))
        f.close()
        with np.errstate(invalid='ignore'):
            vars['sur']['t'][j][np.less(vars['sur']['t'][j],-40)] = np.nan
            vars['sur']['pres'][j][np.less(vars['sur']['pres'][j],900)] = np.nan
            vars['sur']['prec'][j][np.less(vars['sur']['prec'][j],0)] = np.nan

    ### Flatten the dimension of day out ###
    #ranges_idx, ranges_is_idx, vars['kazr']['h'], vars['is']['h'], vars['vpt']['h'], ranges_v_idx, ranges_aeri_idx, vars['aeri']['h'], vars['micro']['h'] = f2d(vars['kazr']['h']),  f2d(vars['is']['h']),        f2d(vars['kazr']['h']), f2d(vars['is']['h']), f2d(vars['vpt']['h']), f2d(vars['vpt']['h']), f2d(vars['aeri']['h']), f2d(vars['aeri']['h']), f2d(vars['micro']['h'])

    
    ranges_idx, ranges_is_idx, vars['kazr']['h'], vars['is']['h'], vars['vpt']['h'], ranges_v_idx, vars['micro']['h'] = f2d(vars['kazr']['h']), f2d(vars['is']['h']), f2d(vars['kazr']['h']), f2d(vars['is']['h']), f2d(vars['vpt']['h']), f2d(vars['vpt']['h']), f2d(vars['micro']['h'])
    
    
    for var in list(vars['vpt'].keys())[1:]:
        vars['vpt'][var] = f2d(vars['vpt'][var])
    for var in list(vars['kazr'].keys())[1:]:
        vars['kazr'][var] = f2d(vars['kazr'][var])
    for var in list(vars['micro'].keys())[1:]:
        vars['micro'][var] = f2d(vars['micro'][var])
    for var in list(vars['is'].keys())[3:]:
        vars['is'][var] = f2d(vars['is'][var])
    for var in vars['cloud'].keys():
        vars['cloud'][var] = f2d(vars['cloud'][var])
    for var in vars['mwr'].keys():
        vars['mwr'][var] = f2d(vars['mwr'][var])
    for var in vars['sur'].keys():
        vars['sur'][var] = f2d(vars['sur'][var])
    #for var in vars['iwp'].keys():
     #   vars['iwp'][var] = f2d(vars['iwp'][var])
   # if case[0] in aeri_days and case[2] in aeri_days:
   #     for var in list(vars['aeri'].keys())[1:]:
   #         vars['aeri'][var] = f2d(vars['aeri'][var])
    if case[0] >= 73:
        #ranges_mpl, ldpr, ranges_mpl_idx, backscatter, time_mpl = f2d(ranges_mpl), f2d(ldpr), f2d(ranges_mpl), f2d(backscatter), f2d(time_mpl)
        vars['cld']['h'], ranges_cld_idx = f2d(vars['cld']['h']), f2d(vars['cld']['h'])
        for var in list(vars['cld'].keys())[1:]:
            vars['cld'][var] = f2d(vars['cld'][var])
        if case[0] >= 91:
            u_vwp, v_vwp, ranges_vwp, ranges_vwp_idx = f2d(u_vwp), f2d(v_vwp), f2d(ranges_vwp), f2d(ranges_vwp)
            u_vwp1, v_vwp1, ranges_vwp1, ranges_vwp_idx1 = f2d(u_vwp1), f2d(v_vwp1), f2d(ranges_vwp1), f2d(ranges_vwp1)

    ### Calc max possible PWV ###
    for j in range(0, len(vars['is']['t'])):
        #vars['is']['q'].append(np.asarray(mpcalc.specific_humidity_from_dewpoint(vars['is']['t'][j]*units.degreeC, vars['is']['pres'][j]*units.hPa)))
        vars['is']['q'].append(np.asarray(mpcalc.specific_humidity_from_dewpoint(vars['is']['pres'][j]*units.hPa, vars['is']['t'][j]*units.degreeC)))
    vars['is']['q'] = np.asarray(vars['is']['q'])

    for j in range(0,len(vars['is']['q'])):
        pwv_max = 0
        for k in range(0, len(vars['is']['q'][j])-1):
            if vars['is']['pres'][j,k+1] > 0 and vars['is']['pres'][j,k] > 0:
                pwv_max += np.mean([vars['is']['q'][j,k],vars['is']['q'][j,k+1]])*(vars['is']['pres'][j,k]-vars['is']['pres'][j,k+1])*100
        vars['is']['pwv'].append(1/9.81*pwv_max)
    vars['is']['pwv'] = np.asarray(vars['is']['pwv'])

    ### Only take data thats within the CAO times ###
    k_i = np.logical_and(vars['kazr']['time']>=bts,vars['kazr']['time']<=bte) # kazr_invalid time
    v_i = np.logical_and(vars['vpt']['time']>=bts,vars['vpt']['time']<=bte) # kasacr vertical invalid time
    s_i = np.logical_and(vars['sur']['time']>=bts,vars['sur']['time']<=bte) # surface_invalid time
    i_i = np.logical_and(vars['is']['time']>=bts,vars['is']['time']<=bte) # interpolated sonde invalid time
    c_i = np.logical_and(vars['cloud']['time']>=bts,vars['cloud']['time']<=bte) # cloud boundaries invalid time
    m_i = np.logical_and(vars['mwr']['time']>=bts,vars['mwr']['time']<=bte) # mwr invalid time
    b_i = np.logical_and(vars['micro']['time']>=bts,vars['micro']['time']<=bte) # microbase invalid time
    for var in list(vars['kazr'].keys())[1:]:
        vars['kazr'][var] = vars['kazr'][var][k_i]
    for var in list(vars['micro'].keys())[1:]:
        vars['micro'][var] = vars['micro'][var][b_i]
    for var in list(vars['vpt'].keys())[1:]:
        vars['vpt'][var] = vars['vpt'][var][v_i]
    for var in list(vars['is'].keys())[1:]:
        vars['is'][var] = vars['is'][var][i_i]
    for var in vars['cloud'].keys():
        vars['cloud'][var] = vars['cloud'][var][c_i]
    #for var in list(vars['iwp'].keys()):
    #    vars['iwp'][var] = vars['iwp'][var][c_i]
    for var in vars['sur'].keys():
        vars['sur'][var] = vars['sur'][var][s_i]
        if var == 'time':
            print(vars['sur'][var])
            vars['sur'][var][-1] = vars['sur'][var][-2] + 60
            
    #if case[0] in aeri_days and case[2] in aeri_days:
    #    a_i = np.logical_and(vars['aeri']['time']>=bts,vars['aeri']['time']<=bte) # aeri invalid time
    #    for var in list(vars['aeri'].keys())[1:]:
    #        vars['aeri'][var] = vars['aeri'][var][a_i]
    if case[0] >= 73:
        p_i = np.logical_and(vars['cld']['time']>=bts,vars['cld']['time']<=bte) ### phase
        for var in list(vars['cld'].keys())[1:]:
            vars['cld'][var] = vars['cld'][var][p_i]
    for var in vars['mwr'].keys():
        vars['mwr'][var] = vars['mwr'][var][m_i]
        if var == 'lwp':
            lwp_i = np.greater(vars['mwr'][var],0)
            vars['mwr'][var] = vars['mwr'][var][lwp_i]
        if var == 'pwv':
            pwv_i = np.greater(vars['mwr'][var],0)
            vars['mwr'][var] = vars['mwr'][var][pwv_i]
        if var == 'time':
            lwp_time = vars['mwr'][var][lwp_i]
            pwv_time = vars['mwr'][var][pwv_i]

   # if case[0]>=73:
       # mpl_i = np.logical_and(time_mpl>=bts,time_mpl<=bte)
        #ldpr, backscatter = ldpr[mpl_i], backscatter[mpl_i]
    u_vwp = u_vwp[int(start_min/60):int((1440*(end_day-start_day)+(end_min+1))/60)]
    v_vwp = v_vwp[int(start_min/60):int((1440*(end_day-start_day)+(end_min+1))/60)]
    u_vwp1 = u_vwp1[int(start_min/60):int((1440*(end_day-start_day)+(end_min+1))/60)]
    v_vwp1 = v_vwp1[int(start_min/60):int((1440*(end_day-start_day)+(end_min+1))/60)]
    time_vwp = []
    for j in range(0,len(u_vwp)):
        time_vwp.append(int(start_min/60)*60-start_min+j*60)
    time_vwp = np.asarray(time_vwp)+60

    ### Remove areas with low signal to noise ratio ###
    noise_invalid = np.less_equal(vars['kazr']['noise'], -15)
    vars['kazr']['ref'][noise_invalid], vars['kazr']['vel'][noise_invalid], vars['kazr']['sw'][noise_invalid] = np.nan, np.nan, np.nan
    noise_invalid_vpt = np.less_equal(vars['vpt']['noise'], -15)
    vars['vpt']['ref'][noise_invalid_vpt], vars['vpt']['vel'][noise_invalid_vpt], vars['vpt']['sw'][noise_invalid_vpt] = np.nan, np.nan, np.nan

    ### Convert time to string for labels ###
    time_surface = []
    for j in range(0,len(vars['sur']['time'])):
        dt_object = datetime.utcfromtimestamp(vars['sur']['time'][j])
        time_surface.append(dt_object.strftime('%H:%M UTC\n%d %b %Y'))
        if j == 0:
            time_stamp = dt_object.strftime('%y%m%d_%H%M')
        if j == len(vars['sur']['time'])-1:
            time_stamp = time_stamp + "_" + dt_object.strftime('%y%m%d_%H%M')
    time_surface = np.asarray(time_surface)

    ### limit data to below certain level ###
    bot, top = 150, 7000
   
    bot_i, top_i = geo_idx.geo_idx(bot, ranges_idx), geo_idx.geo_idx(top, ranges_idx)
    bot_i_v, top_i_v = geo_idx.geo_idx(bot, ranges_v_idx), geo_idx.geo_idx(top, ranges_v_idx)
    bot_i_is, top_i_is  = geo_idx.geo_idx(bot, ranges_is_idx), geo_idx.geo_idx(top, ranges_is_idx)
    #if case[0] in aeri_days and case[2] in aeri_days:
    #    bot_i_aeri, top_i_aeri  = geo_idx.geo_idx(0, ranges_aeri_idx), geo_idx.geo_idx(3300, ranges_aeri_idx)
    if case[0] >= 73:
        #bot_i_MPL, top_i_MPL = geo_idx.geo_idx(bot, ranges_mpl_idx), geo_idx.geo_idx(top, ranges_mpl_idx)
        bot_i_cld, top_i_cld = geo_idx.geo_idx(bot, ranges_cld_idx), geo_idx.geo_idx(top, ranges_cld_idx)
        if case[0] >= 91:
            bot_i_vwp, top_i_vwp = geo_idx.geo_idx(bot, ranges_vwp_idx*1000), geo_idx.geo_idx(top, ranges_vwp_idx*1000)
            bot_i_vwp1, top_i_vwp1  = geo_idx.geo_idx(bot, ranges_vwp_idx1*1000), geo_idx.geo_idx(top, ranges_vwp_idx1*1000)

    for inst in ['kazr','vpt']:
        bot_loop, top_loop = bot_i, top_i
        if inst == 'vpt':
            bot_loop, top_loop = bot_i_v, top_i_v
        for var in ['h','noise','ref','vel','sw']:
            if var == 'h':
                vars[inst][var] = vars[inst][var][bot_loop:top_loop]
            else:
                vars[inst][var] = vars[inst][var][:,bot_loop:top_loop]

    for var in list(vars['is'].keys())[:-1]:
        if var == 'pwv':
            continue
        if var == 'h':
            vars['is'][var] = vars['is'][var][bot_i_is:top_i_is]
        else:
            vars['is'][var] = vars['is'][var][:,bot_i_is:top_i_is]

    reflectivity_cfad = copy.deepcopy(vars['kazr']['ref'])
    velocity_cfad = copy.deepcopy(vars['kazr']['vel'])
    sw_cfad = copy.deepcopy(vars['kazr']['sw'])
    #if case[0] in aeri_days and case[2] in aeri_days:
    #    for var in list(vars['aeri'].keys())[:-3]:
    #        if var == 'h':
    #            vars['aeri'][var] = vars['aeri'][var][bot_i_aeri:top_i_aeri]
    #        else:
    #            vars['aeri'][var] = vars['aeri'][var][:,bot_i_aeri:top_i_aeri]
    if case[0] >= 73:
        #ranges_mpl = ranges_mpl[bot_i_MPL:top_i_MPL]
        #ldpr = ldpr[:,bot_i_MPL:top_i_MPL]
       # backscatter = backscatter[:,bot_i_MPL:top_i_MPL]
        vars['cld']['phase'] = vars['cld']['phase'][:,bot_i_cld:top_i_cld]
        vars['cld']['h'] = vars['cld']['h'][bot_i_cld:top_i_cld]
        if case[0] >= 91:
            ranges_vwp = ranges_vwp[bot_i_vwp:top_i_vwp]
            v_vwp = v_vwp[:,bot_i_vwp:top_i_vwp]
            u_vwp = u_vwp[:,bot_i_vwp:top_i_vwp]
            ranges_vwp1 = ranges_vwp1[bot_i_vwp1:top_i_vwp1]
            v_vwp1 = v_vwp1[:,bot_i_vwp1:top_i_vwp1]
            u_vwp1 = u_vwp1[:,bot_i_vwp1:top_i_vwp1]
            with np.errstate(invalid='ignore'):
                v_vwp[np.logical_or(np.absolute(v_vwp)>100,np.absolute(u_vwp)>100)] = np.nan
                u_vwp[np.logical_or(np.absolute(v_vwp)>100,np.absolute(u_vwp)>100)] = np.nan
                v_vwp1[np.logical_or(np.absolute(v_vwp1)>100,np.absolute(u_vwp1)>100)] = np.nan
                u_vwp1[np.logical_or(np.absolute(v_vwp1)>100,np.absolute(u_vwp1)>100)] = np.nan

    ### aeri temperature perturbation ###
    ### AERI regular grid ###
    #if case[0] in aeri_days and case[2] in aeri_days:
    #    aeri_t_reg,aeri_rh_reg,aeri_time_reg = aeri_reg.aeri_reg(vars['aeri'])
    #    aeri_tbar = np.zeros(np.shape(aeri_t_reg))
    #    for j in range(0,np.shape(aeri_tbar)[0]):
    #        for k in range(0,np.shape(aeri_tbar)[1]):
    #            aeri_tbar[j,k] = aeri_t_reg[j,k]-np.nanmean(aeri_t_reg[:,k])
    t_is_level_mean = np.zeros(np.shape(vars['is']['t'])[1])
    t_is_anomaly = np.zeros(np.shape(vars['is']['t']))
    for j in range(0,np.shape(vars['is']['t'])[1]):
        t_is_level_mean[j]=np.nanmean(vars['is']['t'][:,j])
        t_is_anomaly[:,j] = vars['is']['t'][:,j]-t_is_level_mean[j]

    ### Calculate equivalent potential temperature ###
    eth = mpcalc.equivalent_potential_temperature(vars['is']['pres']*units.hPa,vars['is']['t']*units.degC,vars['is']['dp']*units.degC)
    eth = np.asarray(eth/units.kelvin)
    th = mpcalc.potential_temperature(vars['is']['pres']*units.hPa,vars['is']['t']*units.degC)
    th = np.asarray(th/units.kelvin)
    #q = mpcalc.specific_humidity_from_dewpoint(vars['is']['dp']*units.degC, vars['is']['pres']*units.hPa)
    #q = mpcalc.specific_humidity_from_dewpoint(vars['is']['pres']*units.hPa, vars['is']['dp']*units.degCvars['is']['pres']*units.hPa)
    #q = np.asarray(q*1000/units.dimensionless)
    interp_list = [eth[::-1],t_is_anomaly[::-1], vars['is']['v'][::-1]]

    ### define length of bins for cfad and max,min values for cfad ###
    num_r, num_t = len(vars['kazr']['h']), len(vars['kazr']['ref'])
    maxref_cfad, minref_cfad, maxvel_cfad, minvel_cfad, maxsw_cfad, minsw_cfad = 30.0, -30.0, 1.5, -2.5, 1.0, 0.0

    ### Eliminate unrealistic data and limit data range for cfad ###
    with np.errstate(invalid='ignore'):
        vars['kazr']['ref'][np.logical_or(vars['kazr']['ref']<-40.0,vars['kazr']['ref']>40.0)]=np.nan
        vars['kazr']['vel'][np.logical_or(vars['kazr']['vel']<-6.5,vars['kazr']['vel']>6.5)]=np.nan
        vars['kazr']['sw'][np.logical_or(vars['kazr']['sw']<0.0,vars['kazr']['sw']>11.0)]=np.nan
        vars['vpt']['ref'][np.logical_or(vars['vpt']['ref']<-40.0,vars['vpt']['ref']>40.0)]=np.nan
        vars['vpt']['vel'][np.logical_or(vars['vpt']['vel']<-200,vars['vpt']['vel']>6.5)]=np.nan
        vars['vpt']['sw'][np.logical_or(vars['vpt']['sw']<0.0,vars['vpt']['sw']>11.0)]=np.nan

    ### Cumulatative arrays of all conditions ###
    with np.errstate(invalid='ignore'):
        if i <= len_cases and cfad_plot == True:
            reflectivity_cfad[np.logical_or(reflectivity_cfad<minref_cfad,reflectivity_cfad>maxref_cfad)]=np.nan
            velocity_cfad[np.logical_or(velocity_cfad<minvel_cfad,velocity_cfad>maxvel_cfad)]=np.nan
            sw_cfad[np.logical_or(sw_cfad<minsw_cfad,sw_cfad>maxsw_cfad)]=np.nan
            for j in range(0,num_t):
                invalid_ranges = np.greater(vars['kazr']['h'],vars['cloud']['fcth'][geo_idx.geo_idx(vars['kazr']['time'][j],vars['cloud']['time'])])
                reflectivity_cfad[j,invalid_ranges] = np.nan
                velocity_cfad[j,invalid_ranges] = np.nan
                sw_cfad[j,invalid_ranges] = np.nan
            cfad.cfad(reflectivity_cfad[i], bin_r, ranges, num_r, "ref", time_stamp)
            cfad.cfad(velocity_cfad[i], bin_v, ranges, num_r, "vel", time_stamp)
            cfad.cfad(sw_cfad[i], bin_sw, ranges, num_r, "sw", time_stamp)
            reflectivity_all.append(reflectivity_cfad)
            velocity_all.append(velocity_cfad)
            sw_all.append(sw_cfad)
            if i == len_cases:
                reflectivity_all = f2d(reflectivity_all)
                velocity_all = f2d(velocity_all)
                sw_all = f2d(sw_all)
                bin_r = np.linspace(minref_cfad,maxref_cfad,num=61)
                bin_v = np.linspace(minvel_cfad,maxvel_cfad,num=41)
                bin_sw = np.linspace(minsw_cfad,maxsw_cfad,num=51)
                cfad.cfad(reflectivity_all, bin_r, ranges, num_r, "ref", "_closed_cell")
                cfad.cfad(velocity_all, bin_v, ranges, num_r, "vel", "_cao_all")
                cfad.cfad(sw_all, bin_sw, ranges, num_r, "sw", "_cao_all")

        #cfad_fillplot(reflectivity_s, reflectivity_w, velocity_s, velocity_w, sw_s, sw_w,bin_r, ranges, num_r, "ref")
    ### Run cfad and vertical time series plots for individual cases ###
    # timeseries.timeseries(vars, panel_name, bts, bte,  lwp_time, pwv_time, [aeri_tbar,aeri_rh_reg,aeri_time_reg-bts], time_surface, time_stamp,"ANX_"+folder_save+"_"+data_set)
    timeseries.timeseries(vars, panel_name, bts, bte, lwp_time, pwv_time, [], time_surface, time_stamp,"ANX_"+folder_save+"_"+data_set)
    
    """
    u_500, v_500, v_925, u_925, h_500, h_925, u_mean, v_mean = [], [], [] , [], [], [], [], []
    for j in range(0,len(vars['is']['pres'])):
        is_idx_500 = geo_idx.geo_idx(500.0,vars['is']['pres'][j])
        is_idx_925 = geo_idx.geo_idx(850.0,vars['is']['pres'][j])
        u_500.append(vars['is']['u'][j,is_idx_500])
        u_925.append(vars['is']['u'][j,is_idx_925])
        u_mean.append(np.nanmean(vars['is']['u'][j,is_idx_925:is_idx_500]))
        v_mean.append(np.nanmean(vars['is']['v'][j,is_idx_925:is_idx_500]))
        v_500.append(vars['is']['v'][j,is_idx_500])
        v_925.append(vars['is']['v'][j,is_idx_925])
        h_500.append(vars['is']['h'][is_idx_500])
        h_925.append(vars['is']['h'][is_idx_925])
    du = (np.nanmean(u_500) - np.nanmean(u_925))
    dv = (np.nanmean(v_500) - np.nanmean(v_925))
    du_dz = (du**2+dv**2)**(0.5)/(np.nanmean(h_500)-(np.nanmean(h_925)))
    shear_d = ((-90.0+180.0/np.pi*np.arctan2(du,dv))-(-90.0+180.0/np.pi*np.arctan2(np.nanmean(u_mean),np.nanmean(v_mean))))
    print(du_dz)
    print(shear_d)
    print(np.nanmax(vars['cloud']['fcth']))
    print(np.nanmin(vars['sur']['pres']))
    print(np.nanpercentile(vars['sur']['wspd'], 90))
    print(np.nanmax(vars['sur']['wspd']))
    """
