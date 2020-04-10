import subprocess,sys
import pandas as pd
import numpy as np
from obspy import UTCDateTime

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
            date2 = date1+pd.Timedelta(30, unit='min')
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
