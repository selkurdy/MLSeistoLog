"""
python logsrange.py


"""

# import  os.path
import argparse
# import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def logsminmax(logsdf,indepth=False,hideplot=True):
    allx =logsdf.dropna()
    allwn = allx.WELL.unique().tolist()
    allxcols = allx.columns.tolist()
    dmin = [allx[allx.WELL ==wn ][allxcols[1]].min() for wn in allwn]
    dmax = [allx[allx.WELL ==wn ][allxcols[1]].max() for wn in allwn]
    if not indepth:
        cols1= ['WELL','TIMEMIN','TIMEMAX']
        allwminmax = pd.DataFrame({'WELL':allwn,'TIMEMIN':dmin,'TIMEMAX':dmax})
        allwminmax = allwminmax[cols1].copy()
        xi = [i for i in range(len(allwn)) ]
        fig,ax = plt.subplots(figsize=(8,6))
        ax.plot(allwminmax.TIMEMIN,label='Min T')
        ax.plot(allwminmax.TIMEMAX,label='Max T')
        plt.xticks(xi,allwminmax.WELL,rotation='75')
        fig.legend()
        pdfcl = "logsminmax.pdf"
        if not hideplot:
            plt.show()
        fig.savefig(pdfcl)
        minmaxdf = 'logsminmax.csv'
        allwminmax.to_csv(minmaxdf,index=False)
        print(f'Sucessfully generated {minmaxdf}')
        tstartavg = allwminmax.TIMEMIN.mean()
        tendavg = allwminmax.TIMEMAX.mean()
        print(f'Average Minimum Time: {tstartavg:.2f},  Average Maximum Time: {tendavg:.2f} ')
    else:
        cols1= ['WELL','DEPTHMIN','DEPTHMAX']
        allwminmax = pd.DataFrame({'WELL':allwn,'DEPTHMIN':dmin,'DEPTHMAX':dmax})
        allwminmax = allwminmax[cols1].copy()
        xi = [i for i in range(len(allwn)) ]
        fig,ax = plt.subplots(figsize=(8,6))
        ax.plot(allwminmax.DEPTHMIN,label='Min Z')
        ax.plot(allwminmax.DEPTHMAX,label='Max Z')
        plt.xticks(xi,allwminmax.WELL,rotation='75')
        fig.legend()
        pdfcl = "logsminmax.pdf"
        if not hideplot:
            plt.show()
        fig.savefig(pdfcl)
        minmaxdf = 'logsminmax.csv'
        allwminmax.to_csv(minmaxdf,index=False)
        print(f'Sucessfully generated {minmaxdf}')
        zstartavg = allwminmax.DEPTHMIN.mean()
        zendavg = allwminmax.DEPTHMAX.mean()
        print(f'Average Minimum Depth: {zstartavg:.2f},  Average Maximum Depth: {zendavg:.2f} ')





def getcommandline():
    parser = argparse.ArgumentParser(description='Find and plot depth ranges of all logs')
    parser.add_argument('logfilecsv',help='csv file with all wells generated from preparelogs.py ')
    parser.add_argument('--indepth',action='store_true',default=False,help='logfile is in depth. default=False,i.e. in time')
    parser.add_argument('--hideplot',action='store_true',default=False,help='Hide plots.default=display')

    result=parser.parse_args()
    return result



def main():
    cmdl = getcommandline()
    if cmdl.logfilecsv:
        logsdf = pd.read_csv(cmdl.logfilecsv)
        logsminmax(logsdf,indepth=cmdl.indepth,hideplot= cmdl.hideplot)


if __name__ ==  '__main__':
    main()
