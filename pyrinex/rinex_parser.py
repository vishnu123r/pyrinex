#!/usr/bin/env python
from __future__ import division
from . import Path
from time import time
from datetime import datetime,timedelta
from numpy import ceil, arange,genfromtxt,empty
from pandas import DataFrame,Panel
from io import BytesIO

USEPANDAS=False
"""
basic RINEX parsing
"""


class Rinex():

    def __init__(self,fn):

        tic = time() # DEBUG

        self.fn = Path(fn).expanduser()

        self.whichRinex()
        self.parse_obs_header()
        self.obstimes()

        self.readrinex()

    def readrinex(self):
        blocks = Panel(data = empty((1,1,len(self.obstypes))),
 #                       items=[], #time
#                        major_axis=self.svnames, #satellites
                       minor_axis=self.obstypes) # obs types

        with self.fn.open('r') as f:
            for line in f:
                if 'END OF HEADER' in line[60:]:
                    break

            for line in f:
                Nsv = int(line[29:32]) # how many SV's received at this interval

                obslinespersat = int(ceil(self.Nobstypes/self.maxobsperline))
                Nrows = Nsv*obslinespersat # how many rows to read for this SV at this interval
                satnames = line[32:68]


                for _ in range(int(ceil(Nsv/self.maxsvperline))-1): #if more than one sat line
                    l =   f.readline()
                    line += l
                    satnames+=line[32:68] #FIXME is this right end?

                blocksvnames = self.satnumfixer(self.grouper(satnames,3,Nsv))
                #%% read this INTERVAL's text block
                block = ''.join(f.readline() for _ in range(Nrows))  # Nrows of text as a big string
                btime = self._obstime(line[:26].split()) #always whitespace
                bdf = self._block2df(block,blocksvnames,Nsv)
                try:
                    blocks.loc[btime,blocksvnames,:] = bdf
                except KeyError: #a new SV has come into view
                    for i,n in enumerate(blocksvnames):
                        blocks.loc[btime,n] = bdf[i,:]

    def whichRinex(self):
        """
        input first (top) line of RINEX file to determine what type it is (O,N,G,M)
        """
#%% get first line of file
        with self.fn.open('r') as f:
            line = f.readline()
#%% what type of RINEX file
        self.filetype = line[20] # by RINEx 2.11 definition
        assert self.filetype in ('O','N','G','M'), 'unknown RINEX file type {}'.format(self.filetype)
#%% what RINEX version
        self.version = float(line[:9])

        print('{} is a RINEX {} file.'.format(self.fn,self.version))


    def parse_obs_header(self):
        if int(self.version)==2:
            self.obstypes=[] #might be more than one line
        elif int(self.version)==3:
            self.obstypes={}
            self.Nobstypes={}
            self.sattypes=[] # all sat types in file

        self.svnames = []

        with self.fn.open('r') as f:
            for line in f:
                self.grabfromhead(line)

                if 'END OF HEADER' in line[60:]:
                    break

        assert self.Nobstypes == len(self.obstypes),'number of observation types does not match header'

    def obstimes(self):
        try:
            interval_delta = timedelta(seconds=     int(self.interval),
                               microseconds=int(self.interval % 1)*100000)

            self.Nt = int(ceil((self.tend-self.tstart).total_seconds()/interval_delta.total_seconds()) + 1)

            self.t = self.tstart + interval_delta * arange(self.Nt)
        except AttributeError:
            pass

    def grabfromhead(self,line):
        """
        returns (list of) strings from header based on label, split on whitespace
        header: raw text of header (one big string)
        start: if unused set to None
        end: if unused set to None
        label: header text to match
        """
#%% code common to RINEX 2,3
        assert len(line)<=81,'line in {} was too long'.format(self.fn) # NOTE: s

        ltype = line[60:]

        if 'RINEX VERSION / TYPE' in ltype:
            pass #already handled the first line
        elif 'APPROX POSITION XYZ' in ltype:
            self.site_ecef = (float(line[:14]), float(line[14:28]), float(line[28:42]))
        elif 'INTERVAL' in ltype:
            self.interval = float(line[:10])
        elif '# OF SATELLITES' in ltype:
            self.Nsatellite = int(line[:6])
        elif "TIME OF FIRST OBS" in ltype:
            self.tstart = self._obstime(line[:60].split())
        elif "TIME OF LAST OBS" in ltype:
            self.tend   = self._obstime(line[:60].split())
        elif "PRN / # OF OBS" in ltype: #TODO validate
            self.linespersat = int(ceil(self.Nobstypes / self.maxobstypes))
            assert self.linespersat > 0, 'problem reading number of satellite lines from header'
            if line[3:6]:
                self.svnames += line[3:6]
#%% distinct for RINEX 2,3
        if '{:.2f}'.format(self.version)=='3.01':
            self.maxobstypes = 13
            if "SYS / # / OBS TYPES" in ltype:
                self.sattypes += line[0]
                self.Nobstypes[self.sattypes[-1]] = int(line[3:6])
                self.obstypes[self.sattypes[-1]] = line[6:60].split() # spec guarantees whitespace
        elif '{:.1f}'.format(self.version) == '2.1':
            self.maxobstypes = 9
            self.maxobsperline=5
            self.maxsvperline=12
            if '# / TYPES OF OBSERV' in ltype:
                self.Nobstypes = int(line[:6]) # no more than 9 per line
                self.obstypes += line[6:60].split() # NOTE: Rinex 2.11 spec guarantees whitespace
        else:
            raise NotImplementedError("RINEX version {} is not yet handled".format(self.version))

    def _obstime(self,fol):
        year = int(fol[0])

        if 80<= year <=99:
            year+=1900
        elif year<80: #because we might pass in four-digit year
            year+=2000

        return datetime(year=year, month=int(fol[1]), day= int(fol[2]),
                        hour= int(fol[3]), minute=int(fol[4]),
                        second=int(float(fol[5])),
                        microsecond=int(float(fol[5]) % 1) *100000
                        )

    def satnumfixer(self,satnames):
        return [s[0] + '{:02d}'.format(int(s[1:3])) for s in satnames]

    def grouper(self,txt,n,maxn):
        return [txt[n*i:n+n*i] for i in range(min(len(txt)//n,maxn))]

    def _block2df(self,block,svnames,Nsv):
        """
        input: block of text corresponding to one time increment INTERVAL of RINEX file
        output: 2-D array of float64 data from block. Future: consider whether best to use Numpy, Pandas, or Xray.
        """
        stride=3

        strio = BytesIO(block.encode())
        barr = genfromtxt(strio, delimiter=(14,1,1)*5).reshape((Nsv,-1), order='C')

        data = barr[:,0:self.Nobstypes*stride:stride]
       # lli  = barr[:,1:nobs*stride:stride]
       # ssi  = barr[:,2:nobs*stride:stride]

        if USEPANDAS:
            return DataFrame(index=svnames,columns=self.obstypes, data = data)
        else:
            return data
