import pandas as pd
import numpy as np
from obspy import UTCDateTime

# For sqlx
import subprocess,sys

# For hour map
from matplotlib import colors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patheffects as pe


class PSDs(object):
    def __init__(self,count={},psd={},per={},times=[],mseedids=[]):
        self.count=count
        self.psd=psd
        self.per=per
        self.times=times
        self.mseedids=mseedids
    def add(self,time,mseedid):
        if (mseedid,time) not in self.psd:
            self.count[(mseedid,time)]=[]
            self.psd[(mseedid,time)]=[]
            self.per[(mseedid,time)]=[]
            self.times+=[(mseedid,time)]
            self.mseedids+=[mseedid]
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
        self.displacement_RMS={}
        for mseedid in self.mseedids:
            index = pd.DatetimeIndex(times[mseedid])
            self.displacement_RMS[mseedid] = pd.DataFrame(displacement_RMS[mseedid], 
                                                          index=index)


def sqlx2drms(sshuserhost,
              network = 'CH',
              station = 'SGEV',
              location = '',
              channel = 'HGZ,HGE,HGN',
              dbname = 'AllNetworks',
              start = UTCDateTime()-3*24*60*60,#"2020-03-07")
              end = UTCDateTime(),# means "now"
              freqs = [(0.1,1.0),(1.0,20.0),(4.0,14.0),(4.0,20.0)]):
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
    psds=PSDs()
    datesplit = pd.date_range(start.datetime, 
                              end.datetime, 
                              freq="9D")
    #print(datesplit)
    for date in datesplit:
        if date>= end.datetime:
            break
        ssh = subprocess.Popen(["ssh", 
                                "-i .ssh/id_rsa", 
                                sshuserhost],#sys.argv[1]],
                               stdin =subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True,
                               bufsize=0)
         
        datelist = pd.date_range(date, 
                                 date+pd.Timedelta(9, unit='D'), 
                                 freq="30min")
        #print(datelist)
        # Send ssh commands to stdin
        for date1 in datelist:
            if date1>= end.datetime:
                break
            date2 = date1+pd.Timedelta(minutes=30)
            for n in network.split(','):
                for s in station.split(','):
                    for l in location.split(','):
                        for c in channel.split(','):
                            mseedid = '.'.join([n,s,l,c])
                            command = 'exPSDhour'
                            command += ' AllNetworks'
                            command += ' %s'%mseedid.replace('..','.--.').replace('.',' ')
                            command += ' %s'%date1.strftime("%Y-%m-%d")
                            command += ' %s'%date2.strftime("%Y-%m-%d")
                            command += ' %s'%date1.strftime('%X')
                            command += ' %s'%date2.strftime('%X')
                            command += ' P | sed "s/$/\t%s\tmyprecious/"\n'%(mseedid)
                            #print(date1, date2)
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
                mseedid = data[-1]
                time = UTCDateTime('%s %s'%(data[0],data[1])).datetime
                psds.add(time,mseedid)
                psds.count[(mseedid,time)] += [1]
                psds.psd[(mseedid,time)] += [float(data[3])]
                psds.per[(mseedid,time)] += [float(data[2])]
    return psds
    

def hourmap(data,
            bans=None,
            ax=None,
            scale = 1e9):
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
    :return: A axe with the plot.

    .. rubric:: Basic Usage

    You may omit bans, ax and scale parameters.

    >>> ax = hourmap(data[mseedid])
    """
    width = data.index[1]-data.index[0]
    width = np.pi * 2 / 24 / 60 /60 * width.seconds 
    theta = np.asarray([(d.hour/24+d.minute/60/24)*np.pi*2-width/2 for d in data.index])
    radii = np.asarray([int(d.to_julian_date()+0.5) for d in data.index])
    radii = radii-min(radii)
    norm = colors.Normalize(vmin=np.nanpercentile(data,1),
                            vmax=np.nanpercentile(data,95))
    c_m = plt.cm.viridis
    s_m = plt.cm.ScalarMappable(cmap=c_m, 
                                norm=norm)
    s_m.set_array([])
    valid = np.where(np.isfinite(data))[0][::-1]
    
    if ax is None:
        ax=plt.figure(figsize=(7,9)).add_subplot(111, projection='polar')
    ax.grid(color='w',
            #path_effects=[pe.withStroke(linewidth=2,foreground='w')]
            )
    ax.set_xticks(np.linspace(0,np.pi*2*23/24,24))
    ax.set_xticklabels(['%dh'%h for h in range(24)])
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rmax(max(radii))
    if bans is not None:
        rticks = [(UTCDateTime(ban).datetime - data.index.min().to_pydatetime()).days for iban,ban in enumerate(bans.keys())]
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
                    label=l,
                    color='C%d'%c,
                    path_effects=[pe.withStroke(linewidth=5,
                                                foreground='w'),
                                  pe.withStroke(linewidth=3,
                                                foreground='k')])
    ax.set_yticklabels([])
    ax.set_rorigin(max(radii[valid])/-2)
    ax.text(np.pi,max(radii[valid])/-2,
            data.index[0].strftime("%Y-%m-%d"),
            ha='center',va='center')    
    ax.set_xlabel(data.index[-1].strftime("%Y-%m-%d"))
    plt.legend(loc='lower left',
               bbox_to_anchor= (0.0, -0.2), 
               ncol=2,
               borderaxespad=0, 
               frameon=False)
    cb=plt.colorbar(s_m,orientation='horizontal')#,pad=0.07)
    ticks = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x*scale))
    cb.ax.xaxis.set_major_formatter(ticks)
    cb.ax.set_xlabel("Displacement (nm)")    
    
    ax.bar(theta[valid], radii[valid], 
           color=s_m.to_rgba(data[valid]),
           bottom=radii[valid]-1,
           width=width)
    
    return ax
