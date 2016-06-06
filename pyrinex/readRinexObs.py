#!/usr/bin/env python
"""
RINEX 2 OBS reader
under testing
Michael Hirsch, Greg Starr
MIT License

Program overviw:
1) read the OBS file header:   readHead()
2) parse the OBS file header, obtaining the times, satellites, and data measurement types:   makeSvset()
3) read the OBS files in blocks, where each block is one time interval of data (all sats, all measurements):  makeBlock()

makeBlock dumps the results into a preallocated pandas 3-D Panel with axes:
items / page: time
rows / major_axis: SV
column / minor_axis: data type P1,P2, etc.
"""
from __future__ import division #yes this is needed for py2 here.
from . import Path

USEPANDAS=False # True: DataFrame  False: Numpy NDarray

def rinexobs(obsfn,odir=None,maxtimes=None):
    obsfn = Path(obsfn).expanduser()
    if odir: odir = Path(odir).expanduser()

    if obsfn.suffix.lower().endswith('o'): #raw text file
    #%% save to disk (optional)
        if odir:
            h5fn = odir/obsfn.name.with_suffix('.h5')
            print('saving OBS data to {}'.format(h5fn))
            blocks.to_hdf(h5fn,key='OBS',mode='a',complevel=6,append=False)
    elif obsfn.suffix.lower().endswith('.h5'):
        blocks = read_hdf(obsfn,key='OBS')
        print('loaded OBS data from {} to {}'.format(blocks.items[0],blocks.items[-1]))

    return blocks

def TEC(data,startTime):
    # TODO: update to use datetime()
    TECs=[]
    for d in data:
        difference = []
        for i in range(6):
            difference.append(d[0][i]-startTime[i])
        dt = difference[5]+60*difference[4]+3600*difference[3]+86400*difference[2]
        tec = (9.517)*(10**16)*(d[7]-d[2])
        TECs.append([dt,tec])

    return TECs