#!/usr/bin/python3
from obspy.clients.fdsn import Client
import matplotlib,imp
# to edit text in Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42

import tqdm
import pandas as pd
import numpy as np
from obspy import UTCDateTime

# For pqlx
import subprocess,sys

# For hour map
from matplotlib import colors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patheffects as pe

# For main plot
import os
import datetime
import textwrap
wrapper = textwrap.TextWrapper(width=15,break_long_words=False)
# For maps
#import cartopy.crs as ccrs
#from cartopy.io.img_tiles import OSM
##from mpl_toolkits.basemap import Basemap

class PSDs(object):
    def __init__(self,
                 count={},psd={},per={},times=[],mseedids=[],
                 reloadme=None):
        if reloadme is None:
            self.count=count
            self.psd=psd
            self.per=per
            self.times=times
            self.mseedids=mseedids
        else:
            self.count=reloadme.count
            self.psd=reloadme.psd
            self.per=reloadme.per
            self.times=reloadme.times
            self.mseedids=reloadme.mseedids
            
    def add(self,time,mseedid):
        if (mseedid,time) not in self.psd:
            self.count[(mseedid,time)]=[]
            self.psd[(mseedid,time)]=[]
            self.per[(mseedid,time)]=[]
            self.times+=[(mseedid,time)]
            self.mseedids+=[mseedid]
    
    def clientpqlx(self,
                   sshuserhost='user@hostname',
                   #start = UTCDateTime()-3*24*60*60,
                   #end = UTCDateTime(),
                   **args):

        pqlx2psds(sshuserhost,self=self,**args)
    
    def load(self,
             network = 'CH',
             station = 'SGEV',
             location = '',
             channel = 'HGZ,HGE,HGN',
             start = UTCDateTime()-3*24*60*60,#"2020-03-07")
             end = UTCDateTime(),# means "now"
             freqs = [(0.1,1.0),(1.0,20.0),(4.0,14.0),(4.0,20.0)],
             save='./',
             clientpqlx=True,
             clientobspy=False,
             steps={'clientpqlx':30,'clientobspy':15},
             sshuserhost='user@hostname',
             tocsv=False,
             slow=False,
             output="DISP",
             **args):
        
        self.displacement_RMS = {}
        if clientpqlx:
            step = steps['clientpqlx']
        if clientobspy:
            step = steps['clientobspy']
        if save is not None and not os.path.isdir(save):
            os.makedirs(save)
        loadfile = '%sSeismoSocialDistancing.h5'%save
        store = pd.HDFStore(loadfile)
        for n in network.split(','):
            for s in station.split(','):
                for l in location.split(','):
                    for c in channel.split(','):
                        backfill = [[start.datetime,end.datetime]]
                        mseedid='%s.%s.%s.%s'%(n,s,l,c)
                        if '/'+mseedid.replace('.','_') in store:
                            tmp=store.select(mseedid.replace('.','_'), 
                                             columns=["%.1f-%.1f"%f for f in freqs], 
                                             where=['index>=start.datetime and index<=end.datetime'])
                            if len(tmp)>0:
                                print('Loaded',mseedid,min(tmp.index),max(tmp.index))
                                tlist = pd.date_range(start.datetime, end.datetime, freq="%dmin"%step)
                                backfill=[[UTCDateTime("1920-03-05").datetime]]
                                for i,t in enumerate(tlist):
                                    if min(abs((t-tmp.index).total_seconds())) > step*60*1.5 :
                                        if (t-backfill[-1][-1]).total_seconds() > step*60*1.5:
                                            backfill+=[[t, t]]
                                        else :
                                            backfill[-1][-1]=t
                                backfill = backfill[1:]
                                self.displacement_RMS[mseedid]=tmp
                        if len(backfill)==0:
                            continue
                        if clientpqlx:
                            for bf in backfill:
                                print('Loading',mseedid,bf)
                                tmp = pqlx2psds(sshuserhost,
                                                network = n,
                                                station = s,
                                                location = l,
                                                channel = c,
                                                start = UTCDateTime(bf[0]),
                                                end = UTCDateTime(bf[1]),
                                                **args)
                                print('Computing',mseedid,bf)
                                if slow:
                                    tmp.dRMS(freqs=freqs)
                                else:
                                    tmp.dfRMS(freqs=freqs,output=output)

                                if mseedid not in tmp.displacement_RMS:
                                    print('Missing',mseedid,bf)
                                    continue
                                print('Appending',mseedid,bf)
                                store.append(mseedid.replace('.','_'),
                                             tmp.displacement_RMS[mseedid])
                        print('Selecting',mseedid,(start.datetime,end.datetime))   
                        if '/'+mseedid.replace('.','_') in store:
                            tmp = store.select(mseedid.replace('.','_'),
                                               columns=["%.1f-%.1f"%f for f in freqs],
                                               where=['index>=start.datetime and index<=end.datetime'])
                            self.displacement_RMS[mseedid] = tmp
                            if tocsv:
                                self.displacement_RMS[mseedid].to_csv("%s%s.csv" % (save,mseedid))
                        else:
                            print('Missing',mseedid,(start.datetime,end.datetime))
        store.close()

    def plot(self,
             type='timeseries',
             **args):
        plot(self.displacement_RMS,
             type=type,
             **args)

    def sitemap(self,
                **args):
        plot(self.displacement_RMS,
             type='sitemaps',
             **args)
            
    def clockplot(self,
                  **args):
        plot(self.displacement_RMS,
             type='clockplots',
             **args)
           
    def clockmap(self,
                 **args):
        plot(self.displacement_RMS,
             type='clockmaps',
             **args)

    def gridmap(self,
                 **args):
        plot(self.displacement_RMS,
             type='gridmaps',
             **args)

    def dfRMS(self,
              freqs = [(0.1,1.0),(1.0,20.0),(4.0,14.0),(4.0,20.0)],
              output="DISP"):

        if not hasattr(self,'displacement_RMS'):
            self.displacement_RMS = {}
        displacement_RMS = {}
        mseedids = list(set([mseedid for mseedid,time in self.times])) 
        for mseedid in mseedids:
            ind_times = []
            psd_values = []
            period_bin_centers = []
            for mseedid,time in tqdm.tqdm(self.times):
                if time in self.displacement_RMS[mseedid].index:
                    continue
                period_bin_centers = np.sort(np.unique(self.per[(mseedid,time)]))
                ind_times += [time]
                psd_values += [self.psd[(mseedid,time)][:len(period_bin_centers)]]
            ind_times = pd.DatetimeIndex(ind_times)
            data = pd.DataFrame(psd_values, 
                                index=ind_times, 
                                columns=1.0/np.sort(period_bin_centers))
            data = data.sort_index(axis=1)
            data = df_rms(data, freqs, output=output)
            if mseedid in self.displacement_RMS:
                self.displacement_RMS[mseedid].append(data)
            else:
                self.displacement_RMS[mseedid] = df_rms(data, freqs, output=output)

    def dRMS(self,
             freqs=[(0.1,1.0),
                    (1.0,20.0),
                    (4.0,14.0),
                    (4.0,20.0)]):
        displacement_RMS={}
        times={}
        for mseedid in self.mseedids:
            displacement_RMS[mseedid] = []
            times[mseedid] = []
        for mseedid,time in self.times:
            # acceleration power spectrum in Hz
            f = 1.0/np.sort(self.per[(mseedid,time)])[::-1]
            spec = np.asarray(self.psd[(mseedid,time)])
            spec = spec[np.argsort(self.per[(mseedid,time)])[::-1]]
            # remove NaNs from the list
            valid = np.where(np.isfinite(spec))[0]
            spec = spec[valid]
            f = f[valid]
            w2f = (2.0 * np.pi * f)
            # The acceleration amplitude spectrum (dB to Power! = divide by 10 and not 20!)
            amp = 10.0**(spec/10.)
            # velocity spectrum (divide by omega**2)
            vamp = amp / w2f**2
            # displacement spectrum (divide by omega**2)
            damp =  vamp / w2f**2
            dRMS={}
            for fmin, fmax in freqs:
                ix = np.where((f<=fmax) & (f>=fmin))
                # Parseval: the RMS in time domain is the sqrt of the integral of the power spectrum
                rms = np.sqrt(np.trapz(damp[ix], 
                                       f[ix]))
                frange = "%.1f-%.1f"%(fmin, fmax)
                if rms>0:
                    dRMS[frange] = rms
            if len(list(dRMS.keys())):
                displacement_RMS[mseedid].append(dRMS)
                times[mseedid].append(time)
        self.displacement_RMS = {}
        for mseedid in self.mseedids:
            index = pd.DatetimeIndex(times[mseedid])
            self.displacement_RMS[mseedid] = pd.DataFrame(displacement_RMS[mseedid],index=index)

def dfrms(a):
    return np.sqrt(np.trapz(a.values, a.index))

def df_rms(d, freqs, output="VEL"):
    d = d.dropna(axis=1, how='all')
    RMS = {}
    for fmin, fmax in freqs:
        
        ix = np.where((d.columns>=fmin) & (d.columns<=fmax))[0]
        spec = d.iloc[:,ix]
        f = d.columns[ix]
        
        w2f = (2.0 * np.pi * f)

        # The acceleration power spectrum (dB to Power! = divide by 10 and not 20!)
        amp = 10.0**(spec/10.) 
        if output == "ACC":
            RMS["%.1f-%.1f"%(fmin, fmax)] = amp.apply(dfrms, axis=1)
            continue
        
        # velocity power spectrum (divide by omega**2)
        vamp = amp / w2f**2
        if output == "VEL":
            RMS["%.1f-%.1f"%(fmin, fmax)] = vamp.apply(dfrms, axis=1)
            continue
                
        # displacement power spectrum (divide by omega**2)
        damp = vamp / w2f**2
       
        RMS["%.1f-%.1f"%(fmin, fmax)] = damp.apply(dfrms, axis=1)

    return pd.DataFrame(RMS, index=d.index)#.tz_localize("UTC")#.dropna()

def pqlx2psds(sshuserhost,
              network = 'CH',
              station = 'SGEV',
              location = '',
              channel = 'HGZ,HGE,HGN',
              dbname = 'AllNetworks',
              start = UTCDateTime()-3*24*60*60,#"2020-03-07")
              end = UTCDateTime(),# means "now"
              blocksize = 31*24*2, # equivalent to 9 days 1 channel
              save='./',
              self = None):
    """
    Get PSDs from PQLX
    
    :type sshuserhost: string.
    :param sshuserhost: ssh connection string, e.g. login@hostname.
    :type network,station,location,channel: string.
    :param network,station,location,channel: the mseed codes, use ',' as separator to get several channels.
    :type start, end: `obspy.UTCDateTime``.
    :param start, end: time window.
    :type freqs: list of tuples.
    :param freqs: frequency ranges (one each tuple).
    :return: `PSDs`object.

    .. rubric:: Basic Usage

    You may omit everyhting but sshuserhost.

    >>>myPSDs = sqlx2drms('login@hostname')
    """
    rflag=False
    if self is None:
        rflag=True
        self=PSDs()
    commands = []
    files = []
    datelist = pd.date_range(start.datetime,
                             end.datetime,
                             freq="30min")
    for date1 in datelist:
        date2 = date1+pd.Timedelta(minutes=30)
        date3 = date1+pd.Timedelta(minutes=15)
        if date2 > end.datetime:
            break
        for n in network.split(','):
            for s in station.split(','):
                for l in location.split(','):
                    for c in channel.split(','):
                        savef = '%s/%s/%s/%s%s/%s/%s'%(save,
                                                       n,s,l,c,
                                                       date1.strftime("%Y-%m-%d"),
                                                       date1.strftime('%X').replace(':','.'))
                        mseedid = '.'.join([n,s,l,c])
                        command = 'exPSDhour'
                        command += ' AllNetworks'
                        command += ' %s'%mseedid.replace('..','.--.').replace('.',' ')
                        command += ' %s'%date1.strftime("%Y-%m-%d")
                        command += ' %s'%date2.strftime("%Y-%m-%d")
                        command += ' %s'%date1.strftime('%X')
                        command += ' %s'%date2.strftime('%X')
                        addons = (mseedid,date3.strftime("%Y-%m-%d"),date3.strftime("%X"))
                        command += ' P | sed "s/$/\t%s\t%s\t%s\tmyprecious/"\n'%addons
                        commands += [command]
 
    for c in range(0,len(commands),blocksize):
        ssh = subprocess.Popen(["ssh",
                                "-i .ssh/id_rsa",
                                sshuserhost],#sys.argv[1]],
                               stdin =subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True,
                               bufsize=0)
        stop = c+blocksize
        stop = min([len(commands),stop])
        for cc,command in enumerate(commands[c:stop]):
            ssh.stdin.write(command)
        ssh.stdin.close()

        # Fetch output
        for line in ssh.stdout:        
            if 'myprecious' in line:
                try:
                    data = [v for v in line.strip().split('\t')[:-1]]
                except:
                    print(line.strip(),'unexpected line')
                    continue
                mseedid = data[-3]
                time = UTCDateTime('%s %s'%(data[0],data[1])).datetime
                self.add(time,mseedid)
                self.count[(mseedid,time)] += [1]
                self.psd[(mseedid,time)] += [float(data[3])]
                self.per[(mseedid,time)] += [float(data[2])]
    
    if rflag:
        return self
    
def sitemap(mseedid,
            data_provider='ETH',
            ax=None,
            self=None):

    if self is None:
        self=PSDs()
    if not hasattr(self,'resps'):
        self.resps = {}
    c = Client(data_provider)
    if mseedid not in self.resps:
        msid = mseedid.split('.')
        self.resps[mseedid] = c.get_stations(network=msid[0], 
                                             station=msid[1], 
                                             location=msid[2],
                                             channel=msid[3], 
                                             level="response")
    if ax is None:
        ax=plt.figure(figsize=(7,9)).add_subplot(111)
    longitude = self.resps[mseedid][-1][-1][-1].longitude
    latitude = self.resps[mseedid][-1][-1][-1].latitude
    print(longitude,latitude)
    if False:
        map = Basemap(llcrnrlon=longitude-0.001,
                   llcrnrlat=latitude-0.001,
                   urcrnrlon=longitude+0.001,
                   urcrnrlat=latitude-0.001, 
                   epsg=4326,
                   projection='merc',
                   ax=ax)
        map.arcgisimage(service='ESRI_Imagery_World_2D', 
                     xpixels = 1500, 
                     verbose= True)
    return ax

def pivot_for_hourmap(data, columns="angles"):
    band = data.columns[0]
    data["day"] = [d.year * 365 + d.dayofyear for d in data.index]
    data["time"] = [d.hour + d.minute / 60. for d in data.index]

    data = data.pivot(index="day", columns="time", values=band)
    data.index -= data.index[0]
    data.index = data.index.astype(float)
    if columns == "angles":
        data.columns = 2 * np.pi * data.columns / 24.0
    return data

def hourmap(data,
            bans = {"2020-03-13":'Groups >100 banned',
                    "2020-03-20":'Groups >5 banned'},
            ax=None,
            scale = 1e9,
            unit = 'nm'):
    """
    Make a polar plot of rms

    :type data: dataframe.
    :param data: the rms.
    :type bans: dict.
    :param bans: some annotation, keys are date strings, fields are text desc strings.
    :type ax: axe.
    :param ax: use the provided exiting axe if provided.
    :type scale: float.
    :param scale: scale amplitudes (to nm by default).
    :type unit: string
    :param unit: units for amplitudes (to nm by default).
    :return: A axe with the plot.

    .. rubric:: Basic Usage

    You may omit bans, ax and scale parameters.

    >>> ax = hourmap(data[mseedid])
    """
    origin_time = data.index[0]
    origin_text = data.index[0].strftime("%Y-%m-%d")
    data = data.copy()
    data *= scale

    vmin, vmax = data.quantile(0.01), data.quantile(0.95)
    data = pivot_for_hourmap(data)

    if ax is None:
        ax=plt.figure(figsize=(7,9)).add_subplot(111, projection='polar')

    ax.grid(color='w',
            # path_effects=[pe.withStroke(linewidth=2,foreground='w')]
            )
    ax.set_xticks(np.linspace(0, np.pi * 2 * 23 / 24, 24))
    ax.set_xticklabels(['%d h' % h for h in range(24)])
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    X = np.append(data.columns, 2 * np.pi)
    Y = np.append(data.index, data.index[-1] + 1)

    plt.pcolormesh(X, Y, data, vmax=vmax, vmin=vmin,
                   rasterized=True, antialiased=True)
    cb = plt.colorbar(orientation='horizontal', shrink=0.8)
    cb.ax.set_xlabel("Displacement (%s)" % unit)
    ax.set_rorigin(max(Y) / -4)
    ax.text(np.pi, max(Y) / -4,
            origin_text,
            ha='center', va='center')
    ax.set_xlabel(origin_text)
    ax.grid(color='w',)
    ax.set_rmax(max(Y))

    if bans is not None:
        rticks = [((UTCDateTime(ban).datetime - origin_time.to_pydatetime()).days) for iban, ban in enumerate(bans.keys())]
        xticks = [(UTCDateTime(ban).datetime.hour/24+UTCDateTime(ban).datetime.minute/60/24)*np.pi*2 for iban,ban in enumerate(bans.keys())]
        labels = [bans[iban] for iban in bans.keys()]
        xticks = [xticks[i] for i,d in enumerate(rticks) if d>0]
        labels = [labels[i] for i,d in enumerate(rticks) if d>0]
        rticks = [d for d in rticks if d>0]
        ax.set_rticks(rticks)
        for x,r,l,c in zip(xticks,
                           rticks,
                           labels,
                           range(len(labels))):
            ax.plot(x,r,'o',
                    label='\n'.join(wrapper.wrap(l)),
                    color='C%d'%c,
                    path_effects=[pe.withStroke(linewidth=5,
                                                foreground='w'),
                                  pe.withStroke(linewidth=3,
                                                foreground='k')])

    plt.legend(loc='lower left',
               bbox_to_anchor= (0.0, -0.2), 
               ncol=2,
               borderaxespad=0, 
               frameon=False)

    return ax

def gridmap(data,
            bans = {"2020-03-13":'Groups >100 banned',
                    "2020-03-20":'Groups >5 banned'},
            ax=None,
            scale = 1e9,
            unit = 'nm'):
    """
    Make a polar plot of rms

    :type data: dataframe.
    :param data: the rms.
    :type bans: dict.
    :param bans: some annotation, keys are date strings, fields are text desc strings.
    :type ax: axe.
    :param ax: use the provided exiting axe if provided.
    :type scale: float.
    :param scale: scale amplitudes (to nm by default).
    :type unit: string
    :param unit: units for amplitudes (to nm by default).
    :return: A axe with the plot.

    .. rubric:: Basic Usage

    You may omit bans, ax and scale parameters.

    >>> ax = gridmap(data[mseedid])
    """
    origin_time = data.index[0]
    origin_text = data.index[0].strftime("%Y-%m-%d")
    data = data.copy()
    data *= scale

    vmin, vmax = data.quantile(0.01), data.quantile(0.95)
    days = pd.DatetimeIndex(np.unique(data.index.strftime("%Y-%m-%d")))
    data = pivot_for_hourmap(data, columns='hours')

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(16, 5))

    X = pd.date_range(origin_text, periods=len(data) + 1).to_pydatetime()
    Y = np.append(data.columns, 24)

    plt.pcolormesh(X, Y, data.T, vmax=vmax, vmin=vmin,
                   rasterized=True, antialiased=True)
    plt.colorbar(shrink=0.7, pad=0.01).set_label("Displacement (%s)" % unit)
    ax.set_xticks(pd.date_range(X[0], X[-1], freq="W-MON").to_pydatetime())
    ax.set_yticks(np.arange(25))
    ax.set_yticklabels(['%d h' % h for h in range(25)])

    # fig.autofmt_xdate()
    plt.grid(True, which='both', c="k")
    plt.tight_layout()

    if bans is not None:
        yticks = [UTCDateTime(ban).hour + UTCDateTime(ban).minute/60. for ban in bans]
        xticks = [UTCDateTime(ban[:11]).datetime for ban in bans]
        labels = [bans[iban] for iban in bans.keys()]
        for x,y,l,c in zip(xticks,
                           yticks,
                           labels,
                           range(len(labels))):
            ax.plot(x,y,'o',
                    label='\n'.join(wrapper.wrap(l)),
                    color='C%d'%c,
                    path_effects=[pe.withStroke(linewidth=5,
                                                foreground='w'),
                                  pe.withStroke(linewidth=3,
                                                foreground='k')])

    plt.legend(loc='lower left',
               bbox_to_anchor= (0.0, -0.2),
               ncol=2,
               borderaxespad=0,
               frameon=False)
    plt.gcf().autofmt_xdate()
    return ax

days = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday','Saturday','Sunday']
# Just a bunch of helper functions
def stack_wday_time(df,scale):
    """Takes a DateTimeIndex'ed DataFrame and returns the unstaked table: hours vs day name"""
    return df.groupby(level=(0,1)).median().unstack(level=-1).T.droplevel(0)[days]*scale

def clock24_plot_commons(ax,unit='nm'):
    # Set the circumference labels
    ax.set_xticks(np.linspace(0, 2*np.pi, 24, endpoint=False))
    ax.set_xticklabels(["%i h"%i for i in range(24)], fontsize=8)
    ax.set_yticklabels(["%.2g %s" %(i,unit) for i in ax.get_yticks()], fontsize=7)
    ax.yaxis.set_tick_params(labelsize=8)
    ax.set_rlabel_position(0)

    # Make the labels go clockwise
    ax.set_theta_direction(-1)

    # Place 0 at the top
    ax.set_theta_offset(np.pi/2.0)
    plt.xlabel("Hour (local time)", fontsize=10)
    plt.grid(True)

def radial_hours(N):
    hours = np.deg2rad(np.linspace(0, 360, N-1, endpoint=False))
    hours = np.append(hours, hours[0])
    return hours

def localize_tz_and_reindex(df, freq="15Min", time_zone = "Europe/Brussels"):
    return df.copy().tz_localize("UTC").dropna().tz_convert(time_zone).tz_localize(None).resample(freq).mean().to_frame()
    
def plot(displacement_RMS,
         band = "4.0-14.0",
         logo = 'https://upload.wikimedia.org/wikipedia/commons/thumb/4/44/Logo_SED_2014.png/220px-Logo_SED_2014.png',
         bans = {"2020-03-13":'Groups >100 banned',
                 "2020-03-20":'Groups >5 banned'},
         type = '*',
         scale = 1e9,
         unit = 'nm',
         time_zone = "Europe/Brussels",
         sitedesc = "",# "in Uccle (Brussels, BE)", in original example
         show = True,
         save = None,
         format = 'pdf',
         self = None,
         data_provider='ETH',
         ):
    if save is not None and not os.path.isdir(save):
        os.makedirs(save)

    for channelcode in list(set([k[:-1] for k in displacement_RMS])):
        
        
        data={}
        for o in 'ZEN':
            if channelcode+o not in displacement_RMS :
                continue
            data[channelcode[-2:]+o] = displacement_RMS[channelcode+o][band]
            main=channelcode[-2:]+o
            
        if len(data.keys())>1:
            data[channelcode[-2:]+'*'] = data[main].copy().resample("30min").median().tshift(30, "min") # for the sum
            main=channelcode[-2:]+'*'
            for i,t in enumerate(data[main].index):
                data[main][i] = 0
            for o in data:
                if o == main:
                    continue
                data[o] = data[o].copy().resample("30min" ).median().tshift(30, "min")
                for i,t in enumerate(data[main].index):
                    if len(data[o].index)-1<i:
                        break
                    if True:#abs(data[o].index[i].timestamp()-data[main].index[i].timestamp())<60:
                        data[main][i] += data[o][i]**2
            for i,t in enumerate(data[main].index):
                data[main][i] = data[main][i]**.5

        data[main] = localize_tz_and_reindex(data[main], "30Min", time_zone = time_zone)
        basename = "%s%s-%s"%(save,
                              channelcode[:]+main[-1],
                              band)

        if type in ['*', 'all', 'sitemaps']:
            ax=sitemap(channelcode[:]+main[-1],
                       data_provider=data_provider,
                       self=self)
            if save is not None:
                ax.figure.savefig("%s-map.%s"%(basename,format),
                                  bbox_inches='tight')
            if show:
                plt.show()
                                
        if type in ['*', 'all', 'clockmaps']:
            ax = hourmap(data[main],
                         bans=bans,
                         scale=scale,
                         unit=unit)
            title = 'Seismic Noise for %s - Filter: [%s] Hz' % (channelcode[:]+main[-1],band)
            ax.set_title('Seismic Noise for %s - Filter: [%s] Hz' % (channelcode[:]+main[-1],band))
            if save is not None:
                ax.figure.savefig("%s-hourmap.%s"%(basename,format),
                                  bbox_inches='tight',
                                  facecolor='w')
            if show:
                plt.show()

        if type in ['*', 'all', 'gridmaps']:
            ax = gridmap(data[main],
                         bans=bans,
                         scale=scale,
                         unit=unit)
            title = 'Seismic Noise for %s - Filter: [%s] Hz' % (
            channelcode[:] + main[-1], band)
            ax.set_title('Seismic Noise for %s - Filter: [%s] Hz' % (
            channelcode[:] + main[-1], band))
            if save is not None:
                ax.figure.savefig("%s-gridmap.%s" % (basename, format),
                                  bbox_inches='tight',
                                  facecolor='w')
            if show:
                plt.show()

        if type in ['*', 'all', 'timeseries']:
            fig = plt.figure(figsize=(12,6))
            if logo is not None:
                fig.figimage(plt.imread(logo),
                             40, 40, alpha=.4, zorder=1)
            plt.plot(data[main].index, data[main], label = main)
            
            for o in data:
                rs = data[o].copy().between_time("6:00", "16:00")
                rs = rs.resample("1D" ).median().tshift(12, "H")
                plt.plot(rs.index, rs,
                         label="$\overline{%s}$ (6h-16h)"%o)#, c='purple')

            

            # Get normal business days and set their background color to green
            db = pd.bdate_range(min(data[main].index),
                                max(data[main].index))
            for dbi in db:
                plt.axvspan(dbi, dbi+datetime.timedelta(days=1),
                            facecolor='lightgreen', edgecolor="none",
                            alpha=0.2, zorder=-10)

            plt.ylim(0,np.nanpercentile(data[main],95)*1.5)
            plt.ylim(0,np.nanpercentile(data[main],95)*1.5)
            ticks = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x*scale))
            plt.gca().yaxis.set_major_formatter(ticks)
            plt.ylabel("Displacement (%s)"%unit)

            plt.title('Seismic Noise for %s - Filter: [%s] Hz' % (channelcode[:]+main[-1],
                                                                  band))
            plt.xlim(data[main].index.min(), data[main].index.max())
            fig.autofmt_xdate()
            plt.grid(True, zorder=-1)
            plt.gca().set_axisbelow(True)
            for iban,ban in enumerate(bans.keys()):
                plt.axvline(UTCDateTime(ban).datetime,
                            color='r',
                            linewidth=2,
                            linestyle=['-', '--', '-.', ':'][iban],
                            path_effects=[pe.withStroke(linewidth=4, foreground="k")],
                            zorder=-9,
                            label='\n'.join(wrapper.wrap(bans[ban])))
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            
            ## Idea: add map in an inset below the legend 
            #axins = inset_axes(ax, width="100%", height="100%",
            #                   bbox_to_anchor=(1.05, .6, .5, .4),
            #                   bbox_transform=ax.transAxes, loc=2, borderpad=0)
            #axins.tick_params(left=False, right=True, labelleft=False, labelright=True)
            if save is not None:
                fig.savefig("%s.%s"%(basename,format),
                            bbox_inches='tight',
                            facecolor='w')
            if show:
                plt.show()
        
        if type in ['*', 'all', 'clockplots', 'dailyplots']:
            preloc = data[main].loc[:max(list(bans.keys()))]
            preloc = preloc.set_index([preloc.index.day_name(), preloc.index.hour+preloc.index.minute/60.])
            postloc = data[main].loc[max(list(bans.keys())):]
            postloc = postloc.set_index([postloc.index.day_name(), postloc.index.hour+postloc.index.minute/60.])
            cmap = plt.get_cmap("tab20")

            if type in ['*', 'all', 'dailyplots']:
                ax = stack_wday_time(preloc,scale).plot(figsize=(14,8), cmap = cmap)
                if len(postloc):
                    stack_wday_time(postloc,scale).plot(ls="--", ax=ax, legend=False,cmap = cmap)
                
                plt.title("Daily Noise Levels in %s" % (channelcode[:]+main[-1]))
                plt.ylabel("Amplitude (%s)"%unit)
                plt.xlabel("Hour of day (local time)")
                plt.grid()
                plt.xlim(0,23)
                plt.ylim(0,np.nanpercentile(data[main],95)*1.5*scale)
                if save is not None:
                    ax.figure.savefig("%s-daily.%s"%(basename,format),
                                      bbox_inches='tight',
                                      facecolor='w')
                if show:
                    plt.show()

            if type in ['*', 'all', 'clockplots']:
                # Polar/clock Plot:
                _ = stack_wday_time(preloc,scale).copy()
                _.loc[len(_)+1] = _.iloc[0]
                _.index = radial_hours(len(_))
                
                #subplot_kw = {'polar':True}
                #opts={#'sharey':True,
                #      'figsize':(12,6),
                #      'subplot_kw':subplot_kw}
                #fig, axes  = plt.subplots(1,2,**opts)

                plt.figure(figsize=(12,6))
                ax = plt.subplot(121, polar=True)
                _.plot(ax=ax)#es[0])
    
                plt.title("Before Lockdown", fontsize=12)
                clock24_plot_commons(ax,unit=unit)#es[0])
                ax.set_rmax(np.nanpercentile(data[main],95)*1.5*scale)
                ax.set_rmin(0)
                ax = plt.subplot(122, polar=True, sharey=ax)
                if len(postloc):
                    _ = stack_wday_time(postloc,scale).copy()
                    _.loc[len(_)+1] = _.iloc[0]
                    _.index = radial_hours(len(_))
                    _.plot(ax=ax,#es[0], 
                           ls="--")
    
                plt.title("After Lockdown", fontsize=12)
                clock24_plot_commons(ax,unit=unit)#es[0])
                # ax.set_rmax(np.nanpercentile(data[main],95)*1.5*scale)
                
                suptitle = "Day/Hour Median Noise levels %s\n"
                suptitle += "Station %s - [%s] Hz"
                plt.suptitle(suptitle % (sitedesc,
                                         channelcode[:]+main[-1],
                                         band),
                             fontsize=16)
                plt.subplots_adjust(top=0.80)
                if save is not None:
                    fig = ax.figure
                    fig.savefig("%s-hourly.%s"%(basename,format),
                                bbox_inches='tight',
                                facecolor='w')
                if show:
                    plt.show()
   


if __name__ == "__main__":
    # Include standard modules
    import argparse
    # parse key pairs into a dictionary
    class StoreDictKeyPair(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            my_dict = {}
            for kv in values.split(","):
                k,v = kv.split("=")
                my_dict[k] = v
            setattr(namespace, self.dest, my_dict)
    def listoftup(s):
        try:
            lt = [ tuple([ float(f) for f in l.split('-')]) for l in s.split(',') ] 
            return lt
        except:
            raise 
    # Initiate the parser
    parser = argparse.ArgumentParser()
    # Add long and short argument
    parser.add_argument("--freqs", "-f", 
                        help="set freqs ('4.0-14.0')", 
                        metavar="fmin1-fmax2,fmin2-fmax2,...",
                        type=listoftup,
                        default=[(4.0,14.0)])
    parser.add_argument("--network", "-n", 
                        help="set network ('AA')", 
                        default='CH')
    parser.add_argument("--station", "-s", 
                        help="set station ('CCC,DDD')", 
                        default='SGEV')
    parser.add_argument("--location", "-l", 
                        help="set location ('EE')", 
                        default='')
    parser.add_argument("--channel", "-c", 
                        help="set channel ('FFF,GGG')", 
                        default='HGZ,HGE,HGN')
    parser.add_argument("--begin", "-b", 
                        help="set start time (days from now or date string '2020-03-04')",
                        #type=int, 
                        default=3)
    parser.add_argument("--end", "-e", 
                        help="set end time (days from now or date string '2020-03-07')", 
                        #type=int, 
                        default=0)
    # Arguments for the plots
    parser.add_argument("--type", "-t", 
                        help="set plot type ('*', 'timeseries', 'clockplots', 'clockmaps', 'gridmaps')",
                        default='timeseries')
    parser.add_argument("--output", "-o", 
                        help="save plot (can provide a path)", 
                        default='./')
    parser.add_argument("--extension", "-E", 
                        help="format of the file to save plot (e.g. 'png','pdf')", 
                        default='pdf')
    parser.add_argument("--band", "-F", 
                        help="frequency band for the plot", 
                        default='4.0-14.0')
    parser.add_argument("--logo", "-L", 
                        help="add logo on the plot (a url or path)", 
                        default='https://upload.wikimedia.org/wikipedia/commons/thumb/4/44/Logo_SED_2014.png/220px-Logo_SED_2014.png')
    parser.add_argument("--bans", "-B", 
                        dest="bans", 
                        help="provide dates and label of lockdowns",
                        default={"2020-03-20":'Groups >5 banned',
                                 "2020-03-13":'Groups >100 banned'},
                        action=StoreDictKeyPair, 
                        metavar="DATE1=LABEL1,DATE2=LABEL2...")
    parser.add_argument("--time_zone", "-z", 
                        help="time zone for station (e.g. Europe/Brussels)", 
                        default="Europe/Brussels")
    parser.add_argument("--sitedesc", "-D", 
                       help="site description e.g. 'in Uccle (Brussels, BE)'",                       default="")
    parser.add_argument("--show", "-y", 
                        help="show the plot (True)",
                        default=True, # In any case the default is changed
                        action="store_true")
    parser.add_argument("--noshow", "-Y", 
                        help="do not show the plot (False)", 
                        default=False,
                        action="store_true")
    # Arguments of the PQLX interface
    parser.add_argument("--pqlx", "-p", 
                        help="set PQLX mode", 
                        action="store_true")
    parser.add_argument("--tocsv", "-C", 
                        help="save to csv (False)",
                        default=False, # In any case the default is changed
                        action="store_true")
    parser.add_argument("--slow", "-w",
                        help="Slower RMS computation (False)",
                        default=True, # In any case the default is changed
                        action="store_true")
    parser.add_argument("--sshuserhost", "-S", 
                        help="set ssh parameter (login@hostname)", 
                        default='SQLX')
    parser.add_argument("--dbname", "-d", 
                        help="set dbname, pqlx mode", 
                        default='AllNetworks')
    parser.add_argument("--blocksize", "-x", 
                        help="set blocksize (number PSDs fetched at once)", 
                        type=int,
                        default=31*24*2)
    # Read arguments from the command line
    args = parser.parse_args()
    # Pre-process args
    show=True
    if args.noshow:
        args.show=False
        plt.switch_backend('Agg')
    if not isinstance(args.begin,int):
        args.begin=UTCDateTime(args.begin)
    else:
        args.begin=UTCDateTime()-60*60*24*int(args.begin)
    if not isinstance(args.end,int):
        args.end=UTCDateTime(args.end)
    else:
        args.end=UTCDateTime()-60*60*24*int(args.end)
    args.begin._set_minute(0)
    args.begin._set_second(0)
    args.begin._set_microsecond(0)
    args.end._set_minute(0)
    args.end._set_second(0)
    args.end._set_microsecond(0)
    # Check for --pqlx
    clientpqlx=False
    clientobspy=False
    if args.pqlx:
        clientpqlx=True
    myPSDs = PSDs()
    print(args)
    myPSDs.load(clientpqlx = clientpqlx,
                clientobspy = clientobspy,
                freqs = args.freqs,
                save = args.output,
                network = args.network,
                station = args.station,
                location = args.location,
                channel = args.channel,
                start = args.begin,
                end = args.end,
                sshuserhost=args.sshuserhost,
                dbname = args.dbname,
                blocksize = args.blocksize,
                tocsv = args.tocsv,
                slow = args.slow,
                )
    myPSDs.plot(type=args.type,
                save=args.output,
                band=args.band,
                logo=args.logo,
                bans=args.bans,
                scale=1e9, # Still hardcoded 
                unit='nm', # Idem
                time_zone=args.time_zone,
                sitedesc=args.sitedesc,
                show=args.show,
                format=args.extension,
                )


