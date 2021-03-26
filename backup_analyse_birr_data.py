import pandas as pd 
import matplotlib.pyplot as plt 
from pathlib import Path
from sunpy import timeseries as ts
from sunpy.time import parse_time
from matplotlib import dates
import glob
import numpy as np 
import sys
sys.path.append("..")
from goes_event_list import get_goes_event_list 
from datetime import datetime, timedelta
from scipy import stats
from scipy.optimize import leastsq 
import os
import matplotlib.dates as mdates
from sklearn.metrics import r2_score


def get_flarelist(goes_class_filter, filename):
    """
    getting flare list with flares greater than specifed flare class.
    Saved as a CSV file filename
    """ 
    t_start = "2012-08-22 00:00"
    t_end = "2018-04-20 00:00"
    get_goes_event_list(t_start, t_end, filename=Path.cwd().joinpath(filename), goes_class_filter=goes_class_filter)


# get_flarelist('C1', filename='goes_c_flares_birr_dates.csv')
# get_flarelist('M1', filename='goes_m_flares_birr_dates.csv')
# get_flarelist('C5', filename='goes_c5_flares_birr_dates.csv')

goes_data_dir = "C:/Users/oscar/Desktop/SS_Project/Statistical Study/sid_all/goes_txt_files/"
vlf_data_dir = "C:/Users/oscar/Desktop/SS_Project/Statistical Study/sid_all/birr_sid_archive/"
save_dir = "C:/Users/oscar/Desktop/SS_Project/Statistical Study/sid_all/vlf_plots_birr/"


goes_flares = pd.read_csv("goes_c_flares_birr_dates.csv")
# goes_flares = pd.read_csv("goes_m_flares_birr_dates.csv")
# goes_flares = pd.read_csv("goes_c5_flares_birr_dates.csv")
goes_flares = goes_flares.drop_duplicates(subset="start_time") 

goes_flares["peak_times_hours"] = [x.hour for x in pd.to_datetime(goes_flares["peak_time"])]
daytime_flares = goes_flares[(goes_flares["peak_times_hours"]>10) & (goes_flares["peak_times_hours"]<17)]

days_to_plot = daytime_flares["event_date"].unique()


from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import TimeRange


def goes_data_fetch(filename):
    
    t_start = "2012-08-22 00:00"
    t_end = "2018-04-20 00:00"
    tr = TimeRange([t_start, t_end])
    results = Fido.search(a.Time(tr), a.Instrument.xrs)
    files = Fido.fetch(results, path=goes_data_dir)

# goes_data_fetch(goes_data_dir)

def read_files(files):
    """
    Reads the file(s) in the list `files` and returns
    a pandas dataframe. If the list is > 1 then this
    function concatenates the list data into one dataframe.

    Parameters
    ----------
    files : `list`
        a list contains the filepath to the file(s)

    Returns
    -------
    pd.DataFrame

    """
    if len(files) == 1:
        return pd.read_csv(files[0], comment='#', names=["time", "volts"])

    elif len(files)>1:
        df = []
        for f in files:
            data = pd.read_csv(f, comment='#', names=["time", "volts"])
            df.append(data)
        new_df = pd.concat(df)
        new_df = new_df.drop_duplicates(subset='time')
        new_df.reset_index(drop=True, inplace=True)
        return new_df


def make_vlf_flare_list():
    vlf_days = []
    for i in range(len(days_to_plot)):
        tt = parse_time(days_to_plot[i]).strftime("%Y%m%d")
        files_vlf = glob.glob(vlf_data_dir + tt + '*.csv')
        if len(files_vlf) != 0:
            vlf_days.append(days_to_plot[i])


make_vlf_flare_list()


def plot(i, subtraction=None):

    tt = parse_time(days_to_plot[i]).strftime("%Y%m%d")

    files_vlf = glob.glob(vlf_data_dir + tt + '*.csv')
    if len(files_vlf) == 0:
        print("No VLF data")
        return

    goes_file = goes_data_dir + "go15" + tt + ".fits"
    if not Path(goes_file).exists():
        print("No goes data")
        return

    data_vlf = read_files(files_vlf)
    goes_data = ts.TimeSeries(goes_file).to_dataframe()


    flares_ind = np.where(daytime_flares["event_date"].isin([days_to_plot[i]])==True)[0]
    flares = daytime_flares.iloc[flares_ind]


    fig, ax = plt.subplots(2, figsize=(8,6), sharex=True)

    ax[0].plot(goes_data['xrsb'], color='r', label='GOES 1-8 $\mathrm{\AA}$')
    ax[0].plot(goes_data['xrsa'], color='b', label='GOES 0.5-4 $\mathrm{\AA}$')
    ax[0].legend()
    ax[0].set_yscale('log')
    ax[0].set_xlim(days_to_plot[i] + " 00:00", days_to_plot[i] + " 23:59")
    
    
    if subtraction == None:
        ax[1].plot(pd.to_datetime(data_vlf['time']), data_vlf['volts'], color='grey')
    else:
        baseline = read_files(subtraction)
        ax[1].plot(pd.to_datetime(data_vlf['time']), data_vlf['volts']-baseline['volts'], color='grey')
    
    for f in flares["peak_time"]:
        ax[0].axvline(parse_time(f).datetime, color="k", ls="dashed")
        ax[1].axvline(parse_time(f).datetime, color="k", ls="dashed")

    ax[1].xaxis.set_major_locator(dates.HourLocator(interval=3))
    ax[1].xaxis.set_minor_locator(dates.HourLocator(interval=1))
    ax[1].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
    ax[0].set_ylim(1e-9, 1e-3)
    ax[1].set_ylim(-5.5, 5.5)

    for a in ax:
        a.tick_params(which='both', direction='in')

    ax[0].set_ylabel('Flux Wm$^{-2}$')
    ax[1].set_ylabel('Volts')
    ax[1].set_xlabel('Time ' + days_to_plot[i] + ' UT')
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.05)
    plt.savefig(save_dir + 'birr_vlf_' +  days_to_plot[i] + '.png', dpi=100)
    plt.close()
    

save_dir_flares = "C:/Users/oscar/Desktop/SS_Project/Statistical Study/sid_all/vlf_plots_flares/flares/"

goes_flares = pd.read_csv("goes_c_flares_birr_dates.csv")
goes_flares = goes_flares.drop_duplicates(subset="start_time") 
goes_flares["peak_times_hours"] = [x.hour for x in pd.to_datetime(goes_flares["peak_time"])]
daytime_c_flares = goes_flares[(goes_flares["peak_times_hours"]>10) & (goes_flares["peak_times_hours"]<17)]
c_flare_dates = daytime_c_flares["event_date"].unique()

"""
for i in range(len(days_to_plot)):
    plot(i)
"""

usable = "C:/Users/oscar/Desktop/SS_Project/Statistical Study/sid_all/vlf_plots_birr_usable/"

flare_data_dates = [date[9:-4] for date in os.listdir(usable)]


def start_sub_method():
    vlf_abs_amp = []
    goes_peak_amp = []
    peak_hour = []
    for i in range(len(days_to_plot)):
        if days_to_plot[i] in  flare_data_dates:
                tt = parse_time(days_to_plot[i]).strftime("%Y%m%d")
    
                files_vlf = glob.glob(vlf_data_dir + tt + '*.csv')
                if len(files_vlf) == 0:
                    print("No VLF data")
                    return
            
                goes_file = goes_data_dir + "go15" + tt + ".fits"
                if not Path(goes_file).exists():
                    print("No goes data")
                    return
    
                data_vlf = read_files(files_vlf)
                goes_data = ts.TimeSeries(goes_file).to_dataframe()
                data_vlf['time'] = pd.to_datetime(data_vlf['time'])
                data_vlf = data_vlf.set_index('time')
                data_vlf, goes_data = data_vlf.resample('5S').mean(), goes_data.resample('5S').mean()
                goes_data = goes_data.drop(index = goes_data.index[0])
    
                flares_ind = np.where(daytime_flares["event_date"].isin([days_to_plot[i]])==True)[0]
                flares = daytime_flares.iloc[flares_ind]
                            
                for s, p, e in zip(flares['start_time'].values, flares['peak_time'].values, flares['end_time']):
                    peak_index = data_vlf.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
                    start_index = data_vlf.index.get_loc(datetime.strptime(s, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
                    end_index = data_vlf.index.get_loc(datetime.strptime(e, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
                    flare_peak = data_vlf.iloc[peak_index]
                    flare_start = data_vlf.iloc[start_index]
                    
                    vlf_peak_index = data_vlf['volts'][start_index-300:end_index+300].idxmax()
                    vlf_flare_peak = data_vlf['volts'][vlf_peak_index]
                        
                    vlf_abs_amp.append(abs(flare_peak['volts'] - flare_start['volts']))
                    # vlf_abs_amp.append(abs(flare_peak['volts']))
                    goes_peak = goes_data.iloc[goes_data.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f'), method='nearest')]
                    goes_peak_amp.append(goes_peak.values[1])
                    
                    peak_hour.append(pd.Timestamp(goes_data.index.values[peak_index]).to_pydatetime().hour)
                    
    return goes_peak_amp, vlf_abs_amp, 0 , peak_hour


def inter_sub_method():
    vlf_abs_amp = []
    goes_peak_amp = []
    peak_hour = []
    for i in range(len(days_to_plot)):
        if days_to_plot[i] in  flare_data_dates:
                tt = parse_time(days_to_plot[i]).strftime("%Y%m%d")

                files_vlf = glob.glob(vlf_data_dir + tt + '*.csv')
                if len(files_vlf) == 0:
                    print("No VLF data")
                    return
            
                goes_file = goes_data_dir + "go15" + tt + ".fits"
                if not Path(goes_file).exists():
                    print("No goes data")
                    return

                data_vlf = read_files(files_vlf)
                goes_data = ts.TimeSeries(goes_file).to_dataframe()
                data_vlf['time'] = pd.to_datetime(data_vlf['time'])
                data_vlf = data_vlf.set_index('time')
                data_vlf, goes_data = data_vlf.resample('5S').mean(), goes_data.resample('5S').mean()
                goes_data = goes_data.drop(index = goes_data.index[0])

                flares_ind = np.where(daytime_flares["event_date"].isin([days_to_plot[i]])==True)[0]
                flares = daytime_flares.iloc[flares_ind]
                
            
                for s, p, e in zip(flares['start_time'].values, flares['peak_time'].values, flares['end_time']):
                    peak_index = data_vlf.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
                    start_index = data_vlf.index.get_loc(datetime.strptime(s, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
                    end_index = data_vlf.index.get_loc(datetime.strptime(e, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
                    flare_start = data_vlf.iloc[start_index]
                    flare_end = data_vlf.iloc[end_index+600]
                    flare_peak = data_vlf.iloc[peak_index]
                    interval = end_index - start_index
                    
                    vlf_peak_index = data_vlf['volts'][start_index-300:end_index+300].idxmax()
                    vlf_flare_peak = data_vlf['volts'][vlf_peak_index]
                    
                    baseline = np.linspace(flare_start['volts'], flare_end['volts'], interval)
                    base_amp = baseline[start_index - peak_index]
                    
                    # vlf_abs_amp.append(abs(vlf_flare_peak - base_amp))
                    vlf_abs_amp.append(abs(flare_peak['volts'] - base_amp))
                    
                    goes_peak = goes_data.iloc[goes_data.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f'), method='nearest')]
                    goes_peak_amp.append(goes_peak.values[1])
                    
                    peak_hour.append(pd.Timestamp(goes_data.index.values[peak_index]).to_pydatetime().hour)
                    
    return goes_peak_amp, vlf_abs_amp, 1, peak_hour


def neigh_sub_method():
    vlf_abs_amp = []
    goes_peak_amp = []
    peak_hour = []
    flare_list = pd.DataFrame()
    for i in range(len(days_to_plot)):
        if days_to_plot[i] in  flare_data_dates:
                tt = parse_time(days_to_plot[i]).strftime("%Y%m%d")

                files_vlf = glob.glob(vlf_data_dir + tt + '*.csv')
                if len(files_vlf) == 0:
                    print("No VLF data")
                    return
            
                goes_file = goes_data_dir + "go15" + tt + ".fits"
                if not Path(goes_file).exists():
                    print("No goes data")
                    return

                data_vlf = read_files(files_vlf)
                goes_data = ts.TimeSeries(goes_file).to_dataframe()
                data_vlf['time'] = pd.to_datetime(data_vlf['time'])
                data_vlf = data_vlf.set_index('time')
                data_vlf, goes_data = data_vlf.resample('5S').mean(), goes_data.resample('5S').mean()
                goes_data = goes_data.drop(index = goes_data.index[0])

                flares_ind = np.where(daytime_flares["event_date"].isin([days_to_plot[i]])==True)[0]
                flares = daytime_flares.iloc[flares_ind]
                flare_list = flare_list.append(flares, ignore_index=True)
            
                for s, p, e in zip(flares['start_time'].values, flares['peak_time'].values, flares['end_time']):
                    peak_index = data_vlf.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
                    flare_peak = data_vlf.iloc[peak_index]
                    
                    df = pd.DataFrame(data = {'date': [datetime.strptime(item[0:8], '%Y%m%d').strftime("%Y-%m-%d") for item in os.listdir(vlf_data_dir)], 'filename': os.listdir(vlf_data_dir)}, index=[datetime.strptime(item[0:15], '%Y%m%d_%H%M%S') for item in os.listdir(vlf_data_dir)])
                    df = df[~df['date'].isin(c_flare_dates)]
                    
                    base_date = df.iloc[df.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')]
                    tt = parse_time(base_date['date']).strftime("%Y%m%d")
                    files_vlf_base = glob.glob(vlf_data_dir + tt + '*.csv')
                    baseline = read_files(files_vlf_base)
                    baseline.index = pd.to_datetime(baseline['time'])
                    baseline = baseline.sort_index()
                    p_datetime = datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f')
                    base_time_peak = datetime.strptime(base_date['date'], '%Y-%m-%d') + timedelta(hours=p_datetime.hour, minutes=p_datetime.minute, seconds=p_datetime.second)
                    base_flare_peak = baseline.iloc[baseline.index.get_loc(base_time_peak ,method='nearest')]
                    
                    vlf_abs_amp.append(abs(flare_peak['volts']-base_flare_peak['volts']))
                    
                    goes_peak = goes_data.iloc[goes_data.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f'), method='nearest')]
                    goes_peak_amp.append(goes_peak.values[1])
                    
                    peak_hour.append(pd.Timestamp(goes_data.index.values[peak_index]).to_pydatetime().hour)
                    
    return goes_peak_amp, vlf_abs_amp, 2, peak_hour
                    
def subtraction(method):
    goes_peak_amp, vlf_abs_amp, method_num, peak_hour = method
    
    xy_data = pd.DataFrame({'goes': goes_peak_amp, 'vlf': vlf_abs_amp, 'hour': peak_hour})
    xy_data = xy_data[~xy_data.isin([np.nan, np.inf, -np.inf]).any(1)]

    spear = stats.spearmanr(xy_data['goes'].values, xy_data['vlf'].values)
    fig, ax = plt.subplots(1, figsize=(8,8), sharex=True)
    
    coef = np.polyfit(np.log(xy_data['goes'].values), xy_data['vlf'].values, 1)
    # ax.plot(xy_data['goes'].values, coef[0] * np.log(xy_data['goes'].values) + coef[1], 'g', label='f(x)=' +  str(coef[0])[0:-5] + '*log(x)+' + str(coef[1])[0:-5])
    
    # cb = ax.scatter(xy_data['goes'].values, xy_data['vlf'].values, c=xy_data['hour'])
    ax.scatter(xy_data['goes'].values, xy_data['vlf'].values, color='black', s=10)
    ax.set_ylabel('Absolute VLF Signal (Volts)', fontsize=15)
    ax.set_xlabel('Goes Long x-ray Flux (Wm$^{-2}$)', fontsize=15)
    ax.set_xlim(5e-7, 5e-4)
    ax.xaxis.set_tick_params(labelsize=13)
    ax.tick_params(axis='x', labelsize=13)
    ax.yaxis.set_tick_params(labelsize=13)
    ax.tick_params(axis='y', labelsize=13)
    ax.set_xscale('log')
    ax.legend()
    plt.title("Spearman R Result (correlation = " + str(spear[0])[0:5] + ", p-value = " + str(spear[1])[0:5] + str(spear[1])[-4:] + ")", fontsize=15)
    plt.tight_layout()
    method_name = ['start_subtraction', 'interpolate_subtraction', 'neighbour_subtraction']
    plt.savefig(save_dir_flares + method_name[method_num] + '.png', dpi=100)
    plt.close()
   

#subtraction(method=start_sub_method())
#subtraction(method=neigh_sub_method())
#subtraction(method=inter_sub_method())


def seasonal_plots():
    vlf_abs_amp = []
    vlf_abs_month = []
    goes_peak_amp = []
    for i in range(len(days_to_plot)):
        if days_to_plot[i] in  flare_data_dates:
                tt = parse_time(days_to_plot[i]).strftime("%Y%m%d")

                files_vlf = glob.glob(vlf_data_dir + tt + '*.csv')
                if len(files_vlf) == 0:
                    print("No VLF data")
                    return
            
                goes_file = goes_data_dir + "go15" + tt + ".fits"
                if not Path(goes_file).exists():
                    print("No goes data")
                    return

                data_vlf = read_files(files_vlf)
                goes_data = ts.TimeSeries(goes_file).to_dataframe()

                flares_ind = np.where(daytime_flares["event_date"].isin([days_to_plot[i]])==True)[0]
                flares = daytime_flares.iloc[flares_ind]
            
                for s, p in zip(flares['start_time'].values, flares['peak_time'].values):
                    data_vlf.index = pd.to_datetime(data_vlf['time'])
                    data_vlf = data_vlf.sort_index()
                    flare_peak = data_vlf.iloc[data_vlf.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')]
                    flare_start = data_vlf.iloc[data_vlf.index.get_loc(datetime.strptime(s, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')]

                    #vlf_abs_amp.append(abs(flare_peak['volts'] - flare_start['volts']))
                    vlf_abs_amp.append(abs(flare_peak['volts']))
                    vlf_abs_month.append(datetime.strptime(flare_peak['time'], '%Y-%m-%d %H:%M:%S').replace(year=2020))
                    goes_peak = goes_data.iloc[goes_data.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f'), method='nearest')]
                    goes_peak_amp.append(goes_peak.values[1])

    
    fig, ax = plt.subplots(1, figsize=(24,8), sharex=True) 
    
    df = pd.DataFrame({'month' : vlf_abs_month, 'amp' : vlf_abs_amp})
    df.replace([np.inf, -np.inf], np.nan, inplace=True) 
    df.dropna(inplace=True)
    df = df[np.abs(df.amp-df.amp.mean()) <= (3*df.amp.std())]
    vlf_abs_month, vlf_abs_amp = df['month'].values, df['amp'].values
            
    ax.scatter(vlf_abs_month, vlf_abs_amp, color='black')
    ax.xaxis.set_major_formatter(dates.DateFormatter("%B"))
    ax.xaxis.set_major_locator(dates.MonthLocator(interval=1))
    ax.set_xlabel('Event month', fontsize=15)
    
    x = mdates.date2num(vlf_abs_month)
    xx = np.linspace(x.min(), x.max(), 1000)
    dd = mdates.num2date(xx)
    
    coef = np.polyfit(x, vlf_abs_amp,2)
    polystat = r2_score(np.array(vlf_abs_amp), np.poly1d(coef)(x))
    ax.plot(dd, np.poly1d(coef)(xx), color = 'red', label='3rd poly fit: $R^{2}$ = ' + str(polystat))
    
    guess_mean = np.mean(vlf_abs_amp)
    guess_phase = 0
    guess_freq = 0.005
    guess_amp = 5
    
    optimize_func = lambda l: l[0]*np.sin(l[1]*np.array(x)+l[2]) + l[3] - np.array(vlf_abs_amp)
    est_amp, est_freq, est_phase, est_mean = leastsq(optimize_func, [guess_amp, guess_freq, guess_phase, guess_mean])[0]
    
    y_data = est_amp*np.sin(est_freq*xx+est_phase) + est_mean
    sqstat = r2_score(np.array(vlf_abs_amp), est_amp*np.sin(est_freq*x+est_phase) + est_mean)
    ax.plot(xx, y_data, color='blue', label='Leastsq fit: $R^{2}$ = ' + str(sqstat))
    
    plt.legend(fontsize=18)
     
    ax.set_ylabel('Absolute VLF Signal (Volts)', fontsize=20)
    ax.xaxis.set_tick_params(labelsize=18)
    ax.tick_params(axis='x', labelsize=18)
    ax.yaxis.set_tick_params(labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    plt.tight_layout()
    plt.savefig(save_dir_flares + 'seasonal_dependance' + '.png', dpi=100)
    plt.close()
    

seasonal_plots()


def hour_plots():
    vlf_abs_amp = []
    vlf_abs_hour = []
    goes_peak_amp = []
    for i in range(len(days_to_plot)):
        if days_to_plot[i] in  flare_data_dates:
                tt = parse_time(days_to_plot[i]).strftime("%Y%m%d")

                files_vlf = glob.glob(vlf_data_dir + tt + '*.csv')
                if len(files_vlf) == 0:
                    print("No VLF data")
                    return
            
                goes_file = goes_data_dir + "go15" + tt + ".fits"
                if not Path(goes_file).exists():
                    print("No goes data")
                    return

                data_vlf = read_files(files_vlf)
                goes_data = ts.TimeSeries(goes_file).to_dataframe()

                flares_ind = np.where(daytime_flares["event_date"].isin([days_to_plot[i]])==True)[0]
                flares = daytime_flares.iloc[flares_ind]
            
                for s, p in zip(flares['start_time'].values, flares['peak_time'].values):
                    data_vlf.index = pd.to_datetime(data_vlf['time'])
                    data_vlf = data_vlf.sort_index()
                    flare_peak = data_vlf.iloc[data_vlf.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')]
                    flare_start = data_vlf.iloc[data_vlf.index.get_loc(datetime.strptime(s, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')]

                    # vlf_abs_amp.append(abs(flare_peak['volts'] - flare_start['volts']))
                    vlf_abs_amp.append(abs(flare_peak['volts']))
                    vlf_abs_hour.append(datetime.strptime(flare_peak['time'], '%Y-%m-%d %H:%M:%S').replace(year=2020, month=1, day=1))
                    goes_peak = goes_data.iloc[goes_data.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f'), method='nearest')]
                    goes_peak_amp.append(goes_peak.values[1])

    
    fig, ax = plt.subplots(1, figsize=(24,8), sharex=True) 
    
    df = pd.DataFrame({'hour' : vlf_abs_hour, 'amp' : vlf_abs_amp})
    df.replace([np.inf, -np.inf], np.nan, inplace=True) 
    df.dropna(inplace=True)
    df = df[np.abs(df.amp-df.amp.mean()) <= (3*df.amp.std())]
    vlf_abs_hour, vlf_abs_amp = df['hour'].values, df['amp'].values
            
    ax.scatter(vlf_abs_hour, vlf_abs_amp, color='black')
    ax.xaxis.set_major_formatter(dates.DateFormatter("%H"))
    ax.xaxis.set_major_locator(dates.HourLocator(interval=1))
    ax.set_xlabel('Event Hour')
    
    x = mdates.date2num(vlf_abs_hour)
    xx = np.linspace(x.min(), x.max(), 1000)
    dd = mdates.num2date(xx)
    
    coef = np.polyfit(x, vlf_abs_amp,2)
    polystat = r2_score(np.array(vlf_abs_amp), np.poly1d(coef)(x))
    #ax.plot(xx, np.poly1d(coef)(xx), color = 'red', label='3rd poly fit: $R^{2}$ = ' + str(polystat))
    
    guess_mean = np.mean(vlf_abs_amp)
    guess_phase = 0
    guess_freq = 0.005
    guess_amp = 5
    
    optimize_func = lambda l: l[0]*np.sin(l[1]*np.array(x)+l[2]) + l[3] - np.array(vlf_abs_amp)
    est_amp, est_freq, est_phase, est_mean = leastsq(optimize_func, [guess_amp, guess_freq, guess_phase, guess_mean])[0]
    
    y_data = est_amp*np.sin(est_freq*xx+est_phase) + est_mean
    sqstat = r2_score(np.array(vlf_abs_amp), est_amp*np.sin(est_freq*x+est_phase) + est_mean)
    #ax.plot(xx, y_data, color='blue', label='Leastsq fit: $R^{2}$ = ' + str(sqstat))
    ax.set_ylabel('log')
    
    plt.legend()
     
    ax.set_ylabel('Absolute VLF Signal (Volts)')
    plt.tight_layout()
    plt.savefig(save_dir_flares + 'time_dependance' + '.png', dpi=100)
    plt.close()
    

#hour_plots()


def plot_flare(i):

    tt = parse_time(days_to_plot[i]).strftime("%Y%m%d")

    files_vlf = glob.glob(vlf_data_dir + tt + '*.csv')
    if len(files_vlf) == 0:
        print("No VLF data")
        return

    goes_file = goes_data_dir + "go15" + tt + ".fits"
    if not Path(goes_file).exists():
        print("No goes data")
        return
    
    data_vlf = read_files(files_vlf)
    goes_data = ts.TimeSeries(goes_file).to_dataframe()


    flares_ind = np.where(daytime_flares["event_date"].isin([days_to_plot[i]])==True)[0]
    flares = daytime_flares.iloc[flares_ind]
    
    temp_vlf = data_vlf
    temp_vlf.index = pd.to_datetime(temp_vlf['time'])
    temp_vlf = temp_vlf.sort_index()
    
    
    df = pd.DataFrame(data = {'date': [datetime.strptime(item[0:8], '%Y%m%d').strftime("%Y-%m-%d") for item in os.listdir(vlf_data_dir)], 'filename': os.listdir(vlf_data_dir)}, index=[datetime.strptime(item[0:15], '%Y%m%d_%H%M%S') for item in os.listdir(vlf_data_dir)])
    df = df[~df['date'].isin(c_flare_dates)]
    
    base_date = df.iloc[df.index.get_loc(datetime.strptime(tt, '%Y%m%d') ,method='nearest')]
    tt = parse_time(base_date['date']).strftime("%Y%m%d")
    files_vlf_base = glob.glob(vlf_data_dir + tt + '*.csv')
    neigh_baseline = read_files(files_vlf_base)    
    avg_baseline = pd.read_csv(r'C:/Users/oscar/Desktop/SS_Project/Statistical Study/sid_all/vlf_plots_flares/avg_trend.csv', comment='#', names=["volts"])
    

    for s, e, p in zip(flares['start_time'].values, flares['end_time'].values, flares['peak_time'].values):
        fig, ax = plt.subplots(3, figsize=(8,9), sharex=True)
    
        ax[1].plot(goes_data['xrsb'], color='r', label='1-8$\mathrm{\AA}$')
        ax[1].plot(goes_data['xrsa'], color='b', label='0.5-4$\mathrm{\AA}$')
        ax[1].set_yscale('log')
        ax[1].set_xlim((datetime.strptime(s, '%Y-%m-%dT%H:%M:%S.%f') - timedelta(minutes=5)).strftime("%Y-%m-%d %H:%M"), (datetime.strptime(e, '%Y-%m-%dT%H:%M:%S.%f') + timedelta(minutes=40)).strftime("%Y-%m-%d %H:%M"))
        time = pd.to_datetime(data_vlf['time'])
        
        ax[0].plot(time, data_vlf['volts'], color='black', label='data')
        
        peak_index = data_vlf.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
        start_index = data_vlf.index.get_loc(datetime.strptime(s, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
        end_index = data_vlf.index.get_loc(datetime.strptime(e, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
        flare_start = data_vlf.iloc[start_index]
        flare_end = data_vlf.iloc[end_index+1200]
        flare_peak = data_vlf.iloc[peak_index]
        interval = end_index - start_index
        
        vlf_peak_index = data_vlf['volts'][start_index-300:end_index+300].idxmax()
        vlf_min_index = data_vlf['volts'][start_index-300:end_index+1200].idxmin()
        neigh_min_index = neigh_baseline['volts'][start_index-300:end_index+1200].idxmin()
        neigh_peak_index = neigh_baseline['volts'][start_index-300:end_index+1200].idxmax()
        vlf_flare_peak = data_vlf['volts'][vlf_peak_index]
        vlf_flare_min = data_vlf['volts'][vlf_min_index]
        vlf_flare_time = data_vlf['time'][vlf_peak_index]
        neigh_flare_min = neigh_baseline['volts'][neigh_min_index]
        neigh_flare_peak = neigh_baseline['volts'][neigh_peak_index]
        
        
        #plotted neighbour method baseline
        if len(data_vlf) <= len(neigh_baseline):
            ax[0].plot(time, neigh_baseline['volts'][0:len(data_vlf['volts'])], color = 'blue', label='neighbouring method')
            neigh_sub = data_vlf['volts'].values - neigh_baseline['volts'][0:len(data_vlf['volts'])].values
            ax[2].plot(time, neigh_sub, color='blue', label='neighbouring background')
        else:
            ax[0].plot(time[0:len(neigh_baseline['volts'])], neigh_baseline['volts'], color = 'blue', label='neighbouring background')
            neigh_sub = data_vlf['volts'][0:len(neigh_baseline['volts'])].values - neigh_baseline['volts'].values
            ax[2].plot(time[0:len(neigh_baseline['volts'])], neigh_sub, color='blue',  label='neighbouring method')
        
        #plotted start method baseline 
        start_baseline = np.linspace(flare_start['volts'], flare_start['volts'], len(data_vlf))
        if len(data_vlf) <= len(start_baseline):
            ax[0].plot(time, start_baseline[0:len(data_vlf['volts'])], color = 'red', label='start method')
            start_sub = data_vlf['volts'].values - start_baseline[0:len(data_vlf['volts'])]
            ax[2].plot(time, start_sub, color='red', label='start method')
        else:
            ax[0].plot(time[0:len(start_baseline)], start_baseline, color = 'red', label='start method')
            start_sub = data_vlf[0:len(start_baseline['volts'])] - start_baseline
            ax[2].plot(time[0:len(start_baseline)], start_sub, color='red', label='start method')
        
        #ploting interpolated method baseline
        inter_baseline = np.concatenate([np.linspace(flare_start['volts'], flare_start['volts'], start_index), np.linspace(flare_start['volts'], flare_end['volts'], interval), np.linspace(flare_end['volts'], flare_end['volts'], len(data_vlf) - end_index)])
        ax[0].plot(time, inter_baseline, color='green', label='inter method')
        inter_sub = data_vlf['volts'].values - inter_baseline
        ax[2].plot(time, inter_sub, color='green', label='inter method')
    
        for f in flares["peak_time"]:
            ax[1].axvline(parse_time(f).datetime, color="k", ls="dashed")
            ax[2].axvline(parse_time(f).datetime, color="k", ls="dashed")
            ax[0].axvline(parse_time(f).datetime, color="k", ls="dashed")
            
            ax[1].axvline(datetime.strptime(vlf_flare_time, '%Y-%m-%d %H:%M:%S'), color="grey", ls="dashed")
            ax[2].axvline(datetime.strptime(vlf_flare_time, '%Y-%m-%d %H:%M:%S'), color="grey", ls="dashed")
            ax[0].axvline(datetime.strptime(vlf_flare_time, '%Y-%m-%d %H:%M:%S'), color="grey", ls="dashed")
    
        ax[2].xaxis.set_major_locator(dates.MinuteLocator(interval=10))
        ax[2].xaxis.set_minor_locator(dates.MinuteLocator(interval=2))
        ax[2].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
        ax[1].set_ylim(1e-9, 1e-3)
        ax[0].set_ylim(min([vlf_flare_min, neigh_flare_min]) - 1, max([vlf_flare_peak, neigh_flare_peak]) + 1)
        ax[2].set_ylim(min([min(start_sub[start_index-300:end_index+1200]), min(neigh_sub[start_index-300:end_index+1200]), min(inter_sub[start_index-300:end_index+1200])]) - 1,
                       max([max(start_sub[start_index-300:end_index+1200]), max(neigh_sub[start_index-300:end_index+1200]), max(inter_sub[start_index-300:end_index+1200])]) + 1)
        #ax[1].set_ylim(-5, 5)
    
        for a in ax:
            a.tick_params(which='both', direction='in')
    
        ax[1].set_ylabel('Flux Wm$^{-2}$')
        ax[0].set_ylabel('VLF signal (Volts)')
        ax[2].set_ylabel('Absolute VLF signal after subtraction (Volts)')
        ax[2].set_xlabel('Time ' + days_to_plot[i] + ' UT')
        ax[0].legend()
        ax[2].legend()
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.05)
        plt.savefig(save_dir_flares + 'birr_vlf_flare_' +  parse_time(p).datetime.strftime("%m%d%Y-%H%M%S") + '.png', dpi=100)
        plt.close()


usable = "C:/Users/oscar/Desktop/SS_Project/Statistical Study/sid_all/vlf_plots_birr_usable/"
flare_data_dates = [date[9:-4] for date in os.listdir(usable)]

"""
for i in range(len(days_to_plot)):    
    if days_to_plot[i] in flare_data_dates:
        try:
            plot_flare(i)
        except:
            pass
"""


def goes_vlf_single_flare(i):

    tt = parse_time(days_to_plot[i]).strftime("%Y%m%d")

    files_vlf = glob.glob(vlf_data_dir + tt + '*.csv')
    if len(files_vlf) == 0:
        print("No VLF data")
        return

    goes_file = goes_data_dir + "go15" + tt + ".fits"
    if not Path(goes_file).exists():
        print("No goes data")
        return
    
    data_vlf = read_files(files_vlf)
    goes_data = ts.TimeSeries(goes_file).to_dataframe()
    data_vlf['time'] = pd.to_datetime(data_vlf['time'])
    data_vlf = data_vlf.set_index('time')
    data_vlf, goes_data = data_vlf.resample('5S').mean(), goes_data.resample('5S').mean()
    goes_data = goes_data.drop(index = goes_data.index[0])

    flares_ind = np.where(daytime_flares["event_date"].isin([days_to_plot[i]])==True)[0]
    flares = daytime_flares.iloc[flares_ind]
    
    xy_data = pd.DataFrame({'xrsa': goes_data['xrsa'].values, 'xrsb': goes_data['xrsb'].values, 'vlf': data_vlf['volts'], 'time': data_vlf.index.values})
    xy_data = xy_data[~xy_data.isin([np.nan, np.inf, -np.inf]).any(1)]   
    xy_data = xy_data.set_index('time')
    
    date_to_val = xy_data.index.map(pd.Series(data=np.arange(len(xy_data)), index=xy_data.index.values).to_dict())
        
    for s, e, p in zip(flares['start_time'].values, flares['end_time'].values, flares['peak_time'].values):
   
        start_index = xy_data.index.get_loc(datetime.strptime(s, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
        end_index = xy_data.index.get_loc(datetime.strptime(e, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
        peak_index = xy_data.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
    
        spear = stats.spearmanr(xy_data['xrsa'].values[start_index:end_index], xy_data['vlf'].values[start_index:end_index])
        fig, ax = plt.subplots(1, figsize=(8, 8), sharex=True)
        
        cb = ax.scatter(xy_data['xrsb'].values[start_index-60:peak_index+480], xy_data['vlf'].values[start_index-60:peak_index+480], c=date_to_val[start_index-60:peak_index+480], cmap='brg')
        ax.set_ylabel('VLF Signal Strength (Volts)')
        ax.set_xlabel('Goes Long x-ray Flux (Wm$^{-2}$)')
        plt.title(spear, fontsize=11)
        plt.tight_layout()
        method_name = ['start_subtraction', 'interpolate_subtraction', 'neighbour_subtraction']
        plt.savefig(s[0:10] + '.png', dpi=100)
        plt.close()
    
    
    return


def colour_plot_flare(i):

    tt = parse_time(days_to_plot[i]).strftime("%Y%m%d")

    files_vlf = glob.glob(vlf_data_dir + tt + '*.csv')
    if len(files_vlf) == 0:
        print("No VLF data")
        return

    goes_file = goes_data_dir + "go15" + tt + ".fits"
    if not Path(goes_file).exists():
        print("No goes data")
        return
    
    data_vlf = read_files(files_vlf)
    goes_data = ts.TimeSeries(goes_file).to_dataframe()


    flares_ind = np.where(daytime_flares["event_date"].isin([days_to_plot[i]])==True)[0]
    flares = daytime_flares.iloc[flares_ind]
    
    temp_vlf = data_vlf
    temp_vlf.index = pd.to_datetime(temp_vlf['time'])
    temp_vlf = temp_vlf.sort_index()
    
    
    df = pd.DataFrame(data = {'date': [datetime.strptime(item[0:8], '%Y%m%d').strftime("%Y-%m-%d") for item in os.listdir(vlf_data_dir)], 'filename': os.listdir(vlf_data_dir)}, index=[datetime.strptime(item[0:15], '%Y%m%d_%H%M%S') for item in os.listdir(vlf_data_dir)])
    df = df[~df['date'].isin(c_flare_dates)]
    
    base_date = df.iloc[df.index.get_loc(datetime.strptime(tt, '%Y%m%d') ,method='nearest')]
    tt = parse_time(base_date['date']).strftime("%Y%m%d")
    files_vlf_base = glob.glob(vlf_data_dir + tt + '*.csv')
    neigh_baseline = read_files(files_vlf_base)    
    avg_baseline = pd.read_csv(r'C:/Users/oscar/Desktop/SS_Project/Statistical Study/sid_all/vlf_plots_flares/avg_trend.csv', comment='#', names=["volts"])
    

    for s, e, p in zip(flares['start_time'].values, flares['end_time'].values, flares['peak_time'].values):
        fig, ax = plt.subplots(2, figsize=(8,8), sharex=True)
    
        ax[1].plot(goes_data['xrsb'], color='r', label='1-8$\mathrm{\AA}$')
        ax[1].plot(goes_data['xrsa'], color='b', label='0.5-4$\mathrm{\AA}$')
        ax[1].set_yscale('log')
        ax[1].set_xlim((datetime.strptime(s, '%Y-%m-%dT%H:%M:%S.%f') - timedelta(minutes=5)).strftime("%Y-%m-%d %H:%M"), (datetime.strptime(e, '%Y-%m-%dT%H:%M:%S.%f') + timedelta(minutes=40)).strftime("%Y-%m-%d %H:%M"))
        time = pd.to_datetime(data_vlf['time'])
        
        peak_index = data_vlf.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
        start_index = data_vlf.index.get_loc(datetime.strptime(s, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
        end_index = data_vlf.index.get_loc(datetime.strptime(e, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
        flare_start = data_vlf.iloc[start_index]
        flare_end = data_vlf.iloc[end_index+480]
        flare_peak = data_vlf.iloc[peak_index]
        interval = end_index - start_index
        
        vlf_peak_index = data_vlf['volts'][start_index-60:end_index+300].idxmax()
        vlf_min_index = data_vlf['volts'][start_index-600:end_index+480].idxmin()
        neigh_min_index = neigh_baseline['volts'][start_index-60:end_index+480].idxmin()
        neigh_peak_index = neigh_baseline['volts'][start_index-60:end_index+480].idxmax()
        vlf_flare_peak = data_vlf['volts'][vlf_peak_index]
        vlf_flare_min = data_vlf['volts'][vlf_min_index]
        vlf_flare_time = data_vlf['time'][vlf_peak_index]
        neigh_flare_min = neigh_baseline['volts'][neigh_min_index]
        neigh_flare_peak = neigh_baseline['volts'][neigh_peak_index]
        
        #plotted start method baseline 
        start_baseline = np.linspace(flare_start['volts'], flare_start['volts'], len(data_vlf))
        date_to_val = data_vlf['time'].map(pd.Series(data=np.arange(len(data_vlf)), index=data_vlf['time'].values).to_dict())
        
        if len(data_vlf) <= len(start_baseline):
            start_sub = data_vlf['volts'].values - start_baseline[0:len(data_vlf['volts'])]
        else:
            start_sub = data_vlf[0:len(start_baseline['volts'])] - start_baseline

        cb = ax[0].scatter(time[start_index-60:end_index+480], start_sub[start_index-60:end_index+480],
                               s=5, c=date_to_val[start_index-60:end_index+480].values,
                               label='start method', cmap = 'brg', edgecolor='none') 
        
        for f in flares["peak_time"]:
            ax[1].axvline(parse_time(f).datetime, color="k", ls="dashed")
            ax[0].axvline(parse_time(f).datetime, color="k", ls="dashed")
            
            ax[1].axvline(datetime.strptime(vlf_flare_time, '%Y-%m-%d %H:%M:%S'), color="grey", ls="dashed")
            ax[0].axvline(datetime.strptime(vlf_flare_time, '%Y-%m-%d %H:%M:%S'), color="grey", ls="dashed")

        ax[0].xaxis.set_major_locator(dates.MinuteLocator(interval=10))
        ax[0].xaxis.set_minor_locator(dates.MinuteLocator(interval=2))
        ax[0].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
        ax[1].set_ylim(1e-9, 1e-3)
        ax[0].set_ylim(min([min(start_sub[start_index-60:end_index+480])]) - 1,
                       max([max(start_sub[start_index-60:end_index+480])]) + 1)
        #ax[1].set_ylim(-5, 5)
    
        for a in ax:
            a.tick_params(which='both', direction='in')
    
        ax[1].set_ylabel('Flux Wm$^{-2}$')
        ax[0].set_ylabel('Absolute VLF signal after subtraction (Volts)')
        ax[1].set_xlabel('Time ' + days_to_plot[i] + ' UT')
        ax[0].legend()
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.05)
        plt.savefig('birr_color_lightcurve_' +  parse_time(p).datetime.strftime("%m%d%Y-%H%M%S") + '.png', dpi=100)
        plt.close()

case = "C:/Users/oscar/Desktop/SS_Project/Statistical Study/sid_all/vlf_plots_flares/case_studies/"
flare_data_dates = [date for date in os.listdir(case)]

"""
for i in range(len(days_to_plot)):    
    if days_to_plot[i] in flare_data_dates:
        goes_vlf_single_flare(i)
        colour_plot_flare(i)
"""
      
flare_data_dates = [date[9:-4] for date in os.listdir(usable)]

def csv_flarelist(save=False):
    flares_save = pd.DataFrame()
    for i in range(len(days_to_plot)):
        if days_to_plot[i] in  flare_data_dates:
            flares_ind = np.where(daytime_flares["event_date"].isin([days_to_plot[i]])==True)[0]
            flares = daytime_flares.iloc[flares_ind]
            flares_save = flares_save.append(flares, ignore_index=True)
    
    if save == True:
        flares_save.to_csv('flare_list.csv', index=False)
    else:
        return flares_save
    
#csv_flarelist(save=True)

def flare_position():
    position_data = pd.read_csv("new_flarelist_oscar.csv", comment='#', names=['event_starttime','event_endtime','fl_goescls','SOL_standard','ar_noaanum','hgc_x','hgc_y','hgs_coord','hgs_x','hgs_y','hpc_coord','hpc_radius','hpc_x','hpc_y','hrc_r'])
    
    vlf_abs_amp = []
    hgc_x_abs = []
    flare_list = pd.DataFrame()
    for i in range(len(days_to_plot)):
        if days_to_plot[i] in  flare_data_dates:
                tt = parse_time(days_to_plot[i]).strftime("%Y%m%d")

                files_vlf = glob.glob(vlf_data_dir + tt + '*.csv')
                if len(files_vlf) == 0:
                    print("No VLF data")
                    return
            
                goes_file = goes_data_dir + "go15" + tt + ".fits"
                if not Path(goes_file).exists():
                    print("No goes data")
                    return

                data_vlf = read_files(files_vlf)
                goes_data = ts.TimeSeries(goes_file).to_dataframe()
                
                data_vlf['time'] = pd.to_datetime(data_vlf['time'])
                data_vlf = data_vlf.set_index('time')
                data_vlf, goes_data = data_vlf.resample('5S').mean(), goes_data.resample('5S').mean()
                goes_data = goes_data.drop(index = goes_data.index[0])

                flares_ind = np.where(daytime_flares["event_date"].isin([days_to_plot[i]])==True)[0]
                flares = daytime_flares.iloc[flares_ind]
                flare_list = flare_list.append(flares, ignore_index=True)
                
                for s, p, e in zip(flares['start_time'].values, flares['peak_time'].values, flares['end_time']):
                    if s[0:-4] in position_data['event_starttime'].values:
                        index = np.where(position_data['event_starttime'].isin([s[0:-4]])==True)
                        
                        peak_index = data_vlf.index.get_loc(datetime.strptime(p, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
                        start_index = data_vlf.index.get_loc(datetime.strptime(s, '%Y-%m-%dT%H:%M:%S.%f') ,method='nearest')
                        flare_peak = data_vlf.iloc[peak_index]
                        flare_start = data_vlf.iloc[start_index]

                        vlf_abs_amp.append(abs(flare_peak['volts'] - flare_start['volts']))
                                                
                        hgc_x_abs.append(position_data['hgc_x'].values[index][0])
                    
               
    fig, ax = plt.subplots(1, figsize=(8,8), sharex=True) 
    ax.scatter(hgc_x_abs, vlf_abs_amp, color='k')
    ax.set_xlabel('Heliographic Carrington x coordinate')
    ax.set_ylabel('VLF Absolute Signal Strength')
    plt.savefig('position.png', dpi=100)
    plt.close()
                
    return hgc_x_abs, vlf_abs_amp
    
#u, y = flare_position()
