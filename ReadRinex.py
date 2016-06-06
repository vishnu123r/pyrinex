#!/usr/bin/env python
from matplotlib.pyplot import figure,show,close
#
from pyrinex.rinex_parser import Rinex

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description='example of reading a RINEX 2 Navigation file')
    p.add_argument('rinexfn',help='path to RINEX file')
    p.add_argument('-o','--odir',help='directory in which to write data as HDF5')
    p.add_argument('--noplot',help='Choose to not plot',action='store_false')
    p.add_argument('--profile',help='profile code for debugging',action='store_true')
    p = p.parse_args()

    Data = Rinex(p.rinexfn)


    if p.rinexfn.lower().endswith('n'):
        nav = readRinexNav(rinexfn,p.odir)
        print(nav.head())
    else:
        if p.profile:
            import cProfile
            from pstats import Stats
            profFN = 'RinexObsReader.pstats'
            cProfile.run('rinexobs(rinexfn,p.odir,p.maxtimes)',profFN)
            Stats(profFN).sort_stats('time','cumulative').print_stats(20)
        else:
            blocks = Data.readrinex()
            if not p.noplot:
                pkey = ('P1','C1')
                for k in pkey:
                    try:
                        ax = figure().gca()
                        ax.plot(blocks.items,blocks.ix[:,0,k])
                        ax.set_xlabel('time [UTC]')
                        ax.set_ylabel(k)
                        ax.set_title(k)
                    except KeyError:
                        close()


    show()
