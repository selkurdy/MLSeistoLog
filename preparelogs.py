"""
python preparelogs.py "BA DEEP-001.las" "BA DEEP-001.dev" "BA DEEP-001TDR.txt" BA_DEEP-001 --logheader 32 --logcols 0 4 --logcolname GR --devheader 17 --devflipzsign

**new update : input is las file
python preparelogs.py BA_DEEP-001.las  BA_DEEP-001.dev --logcolname GR --devheader 17

** input is not las, you have to supply logcols and wellname
python preparelogs.py BA_DEEP-001.las  BA_DEEP-001.dev --logcolname GR --devheader 17 --notlas --logheader 32 --logcols 0 4 --wellname BA_DEEP-001
"""

import  os.path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy import interpolate
import lasio


def smooth(a, wlen=11, mode='valid') :
    if wlen % 2 == 0:  #window has to be odd
        wlen +=1
    asmth = savgol_filter(a, wlen, 3) # window size 51, polynomial order 3
    return asmth



def logsminmax(logsdf,hideplot=True):
    allx =all.dropna()
    allwn = allx.WNAME.unique().tolist()
    dmin = [allx[allx.WNAME ==wn ]['DEPTH'].min() for wn in allwn]
    dmax = [allx[allx.WNAME ==wn ]['DEPTH'].max() for wn in allwn]
    cols1= ['WNAME','DEPTHMIN','DEPTHMAX']
    allwminmax = pd.DataFrame({'WNAME':allwn,'DEPTHMIN':dmin,'DEPTHMAX':dmax})
    allwminmax = allwminmax[cols1].copy()
    xi = [i for i in range(len(allwn)) ]
    fig,ax = plt.subplots(figsize=(8,6))
    ax.plot(allwminmax.DEPTHMIN,label='Min Z')
    ax.plot(allwminmax.DEPTHMAX,label='Max Z')
    plt.xticks(xi,allwminmax.WNAME,rotation='75')
    fig.legend()
    pdfcl = "logsminmax.pdf"
    if not hideplot:
        plt.show()
    fig.savefig(pdfcl)
    minmaxdf = 'logsminmax.csv'
    allwminmax.to_csv(minmaxdf,index=False)
    print(f'Sucessfully generated {minmaxdf}')



def getcommandline():
    parser = argparse.ArgumentParser(description='Process logs for mlseistolog')
    parser.add_argument('logfilename',help='Las file name')
    parser.add_argument('devfilename',help='deviation file name')
    parser.add_argument('logcolname',help='log names. single word only. no default')

    parser.add_argument('--dtout',type=int,default=2,help='time sampling interval in ms. default= 2')
    parser.add_argument('--logendtime',type=int,default=2500,help='Maximum T2W to output log.default=2500')

    parser.add_argument('--dzout',type=int,default=2,help='depth sampling interval in ms. default= 2')
    parser.add_argument('--logenddepth',type=int,default=2500,help='Maximum depth to output log.default=2500')

    parser.add_argument('--notlas',action = 'store_true',default=False,help='Log file is not LAS format.default= las' )

    parser.add_argument('--logheader',type=int,default=33,help='header lines of logs file.default=33')
    parser.add_argument('--logcols',type=int,nargs=2,default=[0,1],
        help='columns of log file to read. First col has to be depth.default= 0 1')
    parser.add_argument('--lognull',type=float,default=-999.25,help='logs null value to remove.default= -999.25')
    parser.add_argument('--wellname',help='Well name without spaces. Use if not las format')

    parser.add_argument('--devheader',type=int,default=15,help='Deviation survey header lines. default=15' )
    parser.add_argument('--devcols',type=int,nargs=4,default =[0,1,2,3],help='Deviation Survey MD X Y Z columns.default = 0 1 2 3')
    parser.add_argument('--devflipzsign',action='store_true',default=False,help='Deviation survey flip sign of z values.default: leave as is ')

    parser.add_argument('--tdfilename',help='time depth file name')
    parser.add_argument('--tdcols',type=int,nargs=2,default=[0,1],help='column #s of depth time 1W  pairs.default=0 1 ')
    parser.add_argument('--tdheader',type=int,default=1,help='fb file header lines to skip.default=1')
    parser.add_argument('--flipsign',action='store_true',default=False,help='reverse sign of both d and t1w picks. default=keep as is')
    # parser.add_argument('--tmultiplier',type=float,default=1.0,help='Multiplier applied only to time column.default=1.0')

    parser.add_argument('--logsmoothwlen',type=int,default=61,
        help='smooth logs window length. default=61. freq = Vmax/Lambdamin')
    parser.add_argument('--hideplot',action='store_true',default=False,help='Only save to pdf. default =show and save')


    result=parser.parse_args()
    return result






def main():
    cmdl = getcommandline()
    if cmdl.logfilename:
        logcols = ['DEPTH',cmdl.logcolname]
        if cmdl.notlas:
            logsdf = pd.read_csv(cmdl.logfilename,header=None,usecols = cmdl.logcols,skiprows=cmdl.logheader,delim_whitespace=True)
            # logsdf = pd.DataFrame(data=lgs,columns = logcols)
            # print(logsdf.describe())
            logsdf.columns = logcols
            logsdfx = logsdf.replace(cmdl.lognull,np.nan)
        else:
            logsfile = lasio.read(cmdl.logfilename)
            cmdl.wellname = logsfile.well.WELL.value
            # print(f'well name from las {cmdl.wellname} ')
            alllogsdf = logsfile.df()
            # print(alllogsdf.head())
            print(list(enumerate(alllogsdf.columns)))
            logsdfx = alllogsdf[logcols[1]].copy()
            # logsdf[logcols[0]] = alllogsdf.index
            # logsdf[logcols[1]] = alllogsdf[logcols[1]]
            logsdfx = logsdfx.reset_index()
            logsdfx.columns = logcols

        # print(logsdfx.head())
        logsdfx.dropna(inplace=True)
        print(logsdfx.head())
        # print(logsdfx.describe())
        dirsplitlogs,fextsplitlogs= os.path.split(cmdl.logfilename)
        fnamelogs,fextnlogs= os.path.splitext(fextsplitlogs)

    if cmdl.devfilename:
        devdf = pd.read_csv(cmdl.devfilename,usecols=cmdl.devcols,skiprows=cmdl.devheader,header=None,delim_whitespace=True)

        devkb = devdf.iloc[0,3]
        print('KB {:4.2f}'.format(devkb))
        devcolnames= ['MD','X','Y','Z']
        devdf.columns =devcolnames
        if cmdl.devflipzsign:
            devdf[devcolnames[3]] = devdf[devcolnames[3]].apply(lambda x: x * (-1.0))
        mdx = interpolate.interp1d(devdf.MD.values,devdf.X.values,bounds_error=False,fill_value=np.nan)
        logsdfx['DEVX']= mdx(logsdfx.DEPTH.values)
        mdy = interpolate.interp1d(devdf.MD.values,devdf.Y.values,bounds_error=False,fill_value=np.nan)
        logsdfx['DEVY']= mdy(logsdfx.DEPTH.values)
        logsdfx['WELL'] = cmdl.wellname
        logsdfxcols = ['WELL','DEPTH','DEVX','DEVY',cmdl.logcolname]
        logsdfx = logsdfx[logsdfxcols].copy()
        print(logsdfx.head())
        if cmdl.logsmoothwlen:
            logname = cmdl.logcolname + '_smooth'
            logsdfx[logname] = smooth(logsdfx[cmdl.logcolname],cmdl.logsmoothwlen,mode='same')
            fig,ax = plt.subplots(figsize=(5,7))
            ax.invert_yaxis()
            ax.plot(logsdfx[cmdl.logcolname],logsdfx['DEPTH'],c='r')
            ax.plot(logsdfx[logname],logsdfx['DEPTH'],c='g')
            ax.set_xlabel(cmdl.logcolname)
            ax.set_ylabel('DEPTH')
            ax.set_title(cmdl.wellname)
            pdfcl = os.path.join(dirsplitlogs,fnamelogs) +"_log.pdf"
            if not cmdl.hideplot:
                plt.show()
            fig.savefig(pdfcl)

    if cmdl.tdfilename:
        tddf = pd.read_csv(cmdl.tdfilename,usecols = cmdl.tdcols,header=None,skiprows=cmdl.tdheader,delim_whitespace=True)
        tdcols = ['DEPTH','TIME']
        tddf.columns= tdcols
        tddf['TIME'] = tddf.TIME * 2.0
        # tddf['TIME'] = tddf.TIME * cmdl.tmultiplier
        # print(tddf.head())
        if cmdl.flipsign:
            tddf['TIME'] =tddf.TIME.apply(lambda x : x * (-1.0))
            tddf['DEPTH']=tddf.DEPTH.apply(lambda x : x * (-1.0))
            print(tddf.describe())

        tout = interpolate.interp1d(tddf['TIME'],tddf['DEPTH'],bounds_error=False,fill_value=np.nan)
        ti = np.arange(0,cmdl.logendtime,cmdl.dtout)
        zatt = tout(ti)
        mdx = interpolate.interp1d(devdf.MD.values,devdf.X.values,bounds_error=False,fill_value=np.nan)
        devx_att= mdx(zatt)
        mdy = interpolate.interp1d(devdf.MD.values,devdf.Y.values,bounds_error=False,fill_value=np.nan)
        devy_att= mdy(zatt)

        logv = interpolate.interp1d(logsdfx['DEPTH'],logsdfx[cmdl.logcolname] ,bounds_error=False,fill_value=np.nan)
        logv_att= logv(zatt)

        logstcols = ['WELL','TIME','DEVX','DEVY',cmdl.logcolname]
        logstdfx = pd.DataFrame({'WELL':cmdl.wellname,'TIME':ti,'DEVX':devx_att,'DEVY':devy_att,cmdl.logcolname:logv_att})
        logstdfx = logstdfx[logstcols].copy()
        print(logstdfx.head())

        logsfnamecsv = 'All'+ cmdl.logcolname +'.csv'
        if os.path.isfile(logsfnamecsv):
            logstdfx.to_csv(logsfnamecsv,header=False,index=False,mode='a')
            print('Sucessfully appended {}'.format(logsfnamecsv))
        else:
            logstdfx.to_csv(logsfnamecsv,index=False)
            print('Sucessfully generated {}'.format(logsfnamecsv))


    else: # no depth time pairs, i.e. all wells in depth
        zi = np.arange(0,cmdl.logenddepth,cmdl.dzout)
        mdx = interpolate.interp1d(devdf.MD.values,devdf.X.values,bounds_error=False,fill_value=np.nan)
        devx_atz= mdx(zi)
        mdy = interpolate.interp1d(devdf.MD.values,devdf.Y.values,bounds_error=False,fill_value=np.nan)
        devy_atz= mdy(zi)

        logv = interpolate.interp1d(logsdfx['DEPTH'],logsdfx[cmdl.logcolname] ,bounds_error=False,fill_value=np.nan)
        logv_atz= logv(zi)

        logszcols = ['WELL','DEPTH','DEVX','DEVY',cmdl.logcolname]
        logszdfx = pd.DataFrame({'WELL':cmdl.wellname,'DEPTH':zi,'DEVX':devx_atz,'DEVY':devy_atz,cmdl.logcolname:logv_atz})
        logszdfx = logszdfx[logszcols].copy()
        print(logszdfx.head())



        logsfnamecsv = 'All'+ cmdl.logcolname +'.csv'
        if os.path.isfile(logsfnamecsv):
            logszdfx.to_csv(logsfnamecsv,header=False,index=False,mode='a')
            print('Sucessfully appended {}'.format(logsfnamecsv))
        else:
            logszdfx.to_csv(logsfnamecsv,index=False)
            print('Sucessfully generated {}'.format(logsfnamecsv))







        # logsfnametxt = fnamelogs +'.txt'
        # logsdfx.to_csv(logsfnamecsv,index=False)
        # print('Sucessfully generated {}'.format(logsfnamecsv))
        # logsdfx.to_csv(logsfnametxt,index=False,sep = ' ')
        # print('Sucessfully generated {}'.format(logsfnametxt))
        # logsdfxs = MinMaxScalercaler().fit_transform(logsdfx.iloc[:,1:].values)


if __name__ == '__main__':
    main()
