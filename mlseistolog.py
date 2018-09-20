"""
Read in logs in ASCII flat file
Read in TDR in ASCII flat file
Extract depth interval of interest

python mlseistolog.py allsegy.txt --startendslice 200 202
python mlseistolog.py allsegy.txt vclay.csv --startendslice 1000 1020
# after changing from slice # to depth
python mlseistolog.py allsegy.txt vclay.csv --startendslice 1000 1020 --hideplots
python mlseistolog.py allsegy.txt vclay.csv --startendslice 1000 2000 --hideplots

python mlseistolog.py allsgy.txt allgr.csv

"""

import os.path
import argparse
from datetime import datetime
import numpy as np
import scipy.signal as sg
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import scipy.stats as sts
import segyio as sg
from shutil import copyfile
from scipy.spatial import cKDTree as KDTree
from scipy.interpolate import griddata,Rbf,LinearNDInterpolator,CloughTocher2DInterpolator
from matplotlib.backends.backend_pdf import PdfPages
import itertools as it
from scipy.signal import savgol_filter
from scipy import interpolate
from sklearn.preprocessing import StandardScaler
from sklearn import linear_model
from sklearn.metrics import r2_score, mean_squared_error

try:
    from catboost import CatBoostRegressor
except ImportError:
    print('***Warning:CatBoost is not installed')



def pip(x,y,poly):
   # check if point is a vertex
   """
   if (x,y) in poly:
        return True
    """
   # check if point is on a boundary
   for i in range(len(poly)):
      p1 = None
      p2 = None
      if i==0:
         p1 = poly[0]
         p2 = poly[1]
      else:
         p1 = poly[i-1]
         p2 = poly[i]
      if p1[1] == p2[1] and p1[1] == y and x > min(p1[0], p2[0]) and x < max(p1[0], p2[0]):
         return True

   n = len(poly)
   inside = False

   p1x,p1y = poly[0]
   for i in range(n+1):
      p2x,p2y = poly[i % n]
      if y > min(p1y,p2y):
         if y <= max(p1y,p2y):
            if x <= max(p1x,p2x):
               if p1y != p2y:
                  xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
               if p1x == p2x or x <= xints:
                  inside = not inside
      p1x,p1y = p2x,p2y

   if inside: return True
   else: return False


#...............................................................................
class Invdisttree:
    """ inverse-distance-weighted interpolation using KDTree:
    invdisttree = Invdisttree( X, z )  -- data points, values
    interpol = invdisttree( q, nnear=3, eps=0, p=1, weights=None, stat=0 )
    interpolates z from the 3 points nearest each query point q;
    For example, interpol[ a query point q ]
    finds the 3 data points nearest q, at distances d1 d2 d3
    and returns the IDW average of the values z1 z2 z3
        (z1/d1 + z2/d2 + z3/d3)
        / (1/d1 + 1/d2 + 1/d3)
        = .55 z1 + .27 z2 + .18 z3  for distances 1 2 3

    q may be one point, or a batch of points.
    eps: approximate nearest, dist <= (1 + eps) * true nearest
    p: use 1 / distance**p
    weights: optional multipliers for 1 / distance**p, of the same shape as q
    stat: accumulate wsum, wn for average weights

    How many nearest neighbors should one take ?
    a) start with 8 11 14 .. 28 in 2d 3d 4d .. 10d; see Wendel's formula
    b) make 3 runs with nnear= e.g. 6 8 10, and look at the results --
        |interpol 6 - interpol 8| etc., or |f - interpol*| if you have f(q).
        I find that runtimes don't increase much at all with nnear -- ymmv.

    p=1, p=2 ?
        p=2 weights nearer points more, farther points less.
        In 2d, the circles around query points have areas ~ distance**2,
        so p=2 is inverse-area weighting. For example,
            (z1/area1 + z2/area2 + z3/area3)
            / (1/area1 + 1/area2 + 1/area3)
            = .74 z1 + .18 z2 + .08 z3  for distances 1 2 3
        Similarly, in 3d, p=3 is inverse-volume weighting.

    Scaling:
        if different X coordinates measure different things, Euclidean distance
        can be way off.  For example, if X0 is in the range 0 to 1
        but X1 0 to 1000, the X1 distances will swamp X0;
        rescale the data, i.e. make X0.std() ~= X1.std() .

    A nice property of IDW is that it's scale-free around query points:
    if I have values z1 z2 z3 from 3 points at distances d1 d2 d3,
    the IDW average
        (z1/d1 + z2/d2 + z3/d3)
        / (1/d1 + 1/d2 + 1/d3)
    is the same for distances 1 2 3, or 10 20 30 -- only the ratios matter.
    In contrast, the commonly-used Gaussian kernel exp( - (distance/h)**2 )
    is exceedingly sensitive to distance and to h.

        """
    # anykernel( dj / av dj ) is also scale-free
    # error analysis, |f(x) - idw(x)| ? todo: regular grid, nnear ndim+1, 2*ndim

    def __init__( self, X, z, leafsize=10, stat=0 ):
        assert len(X) == len(z), "len(X) %d != len(z) %d" % (len(X), len(z))
        self.tree = KDTree( X, leafsize=leafsize )  # build the tree
        self.z = z
        self.stat = stat
        self.wn = 0
        self.wsum = None;

    def __call__( self, q, nnear=6, eps=0, p=1, weights=None ):
            # nnear nearest neighbours of each query point --
        q = np.asarray(q)
        qdim = q.ndim
        if qdim == 1:
            q = np.array([q])
        if self.wsum is None:
            self.wsum = np.zeros(nnear)

        self.distances, self.ix = self.tree.query( q, k=nnear, eps=eps )
        interpol = np.zeros( (len(self.distances),) + np.shape(self.z[0]) )
        jinterpol = 0
        for dist, ix in zip( self.distances, self.ix ):
            if nnear == 1:
                wz = self.z[ix]
            elif dist[0] < 1e-10:
                wz = self.z[ix[0]]
            else:  # weight z s by 1/dist --
                w = 1 / dist**p
                if weights is not None:
                    w *= weights[ix]  # >= 0
                w /= np.sum(w)
                wz = np.dot( w, self.z[ix] )
                if self.stat:
                    self.wn += 1
                    self.wsum += w
            interpol[jinterpol] = wz
            jinterpol += 1
        return interpol if qdim > 1  else interpol[0]


def idw(xy,vr,xyi):

    # N = vr.size
    # Ndim = 2
    # Nask = N  # N Nask 1e5: 24 sec 2d, 27 sec 3d on mac g4 ppc
    Nnear = 8  # 8 2d, 11 3d => 5 % chance one-sided -- Wendel, mathoverflow.com
    leafsize = 10
    eps = .1  # approximate nearest, dist <= (1 + eps) * true nearest
    p = 2  # weights ~ 1 / distance**p
    invdisttree = Invdisttree( xy, vr, leafsize=leafsize, stat=1 )
    interpol = invdisttree( xyi, nnear=Nnear, eps=eps, p=p )
    return interpol

def qhull(sample):
    link = lambda a,b: np.concatenate((a,b[1:]))
    edge = lambda a,b: np.concatenate(([a],[b]))
    def dome(sample,base):
        h, t = base
        dists = np.dot(sample-h, np.dot(((0,-1),(1,0)),(t-h)))
        outer = np.repeat(sample, dists>0, 0)
        if len(outer):
            pivot = sample[np.argmax(dists)]
            return link(dome(outer, edge(h, pivot)),
                        dome(outer, edge(pivot, t)))
        else:
            return base

    if len(sample) > 2:
        axis = sample[:,0]
        base = np.take(sample, [np.argmin(axis), np.argmax(axis)], 0)
        return link(dome(sample, base),dome(sample, base[::-1]))
    else:
        return sample

def gridlistin(fname,xyvcols=[0,1,2],nheader=0): #used for single coef per file
    xyv=np.genfromtxt(fname,usecols=xyvcols,skip_header=nheader)
    #filter surfer null values by taking all less than 10000, arbitrary!!
    xyv = xyv[xyv[:,2]<10000.0]
    #xya = xya[~xya[:,2]==  missing]
    return xyv[:,0],xyv[:,1],xyv[:,2]


def map2ddata(xy,vr,xyi,radius=5000.0,maptype='idw'):
    # stats=sts.describe(vr)
    # statsstd=sts.tstd(vr)
    if maptype == 'idw':
        vri=idw(xy,vr,xyi)
    elif maptype =='nearest':
        vri=griddata(xy,vr,(xyi[:,0],xyi[:,1]),method='nearest')
    elif maptype == 'linear':
        #                vri=griddata(xy,vr,(xyifhull[:,0],xyifhull[:,1]),method='linear')
        vri=griddata(xy,vr,(xyi[:,0],xyi[:,1]),method='linear')
    elif maptype == 'cubic':
        vri=griddata(xy,vr,(xyi[:,0],xyi[:,1]),method='cubic')
    elif maptype =='rbf':
        rbf=Rbf(xy[:,0],xy[:,1],vr)
        vri= rbf(xyi[:,0],xyi[:,1])
    # elif maptype =='avgmap':
    #     vri=dataavgmap(xy,vr,xyi,radius)
    elif maptype =='triang':
        linearnd=LinearNDInterpolator(xy,vr,stats[2])
        vri= linearnd(xyi)
    elif maptype == 'ct':
        ct=CloughTocher2DInterpolator(xy,vr,stats[2])
        vri=ct(xyi)
    return vri



def filterhullpolygon(x,y,polygon):
    xf=[]
    yf=[]
    for i in range(x.size):
        if pip(x[i],y[i],polygon):
            xf.append(x[i])
            yf.append(y[i])
    return np.array(xf),np.array(yf)

def filterhullpolygon_mask(x,y,polygon):
    ma=[]
    for i in range(x.size):
        if pip(x[i],y[i],polygon):
            ma.append(True)
        else:
            ma.append(False)
    return np.array(ma)


def process_segylist(segyflist):
    segylist = []
    with open(segyflist,'r') as f:
        for line in f:
            segylist.append(line.rstrip())
    return segylist

def get_slice(fname,slicenum):
    slc = list()
    with sg.open(fname,'r',ignore_geometry= True) as srcp:
        for trnum,tr in enumerate(srcp.trace):
            slc.append(tr[slicenum])
    return slc

def zero_segy(fname):
    with sg.open(fname,'r+',ignore_geometry= True) as srcp:
        for trnum,tr in enumerate(srcp.trace):
            srcp.trace[trnum] = tr * 0

def get_samplerate(fname):
    with sg.open(fname,'r',ignore_geometry= True) as srcp:
        hdrdict = dict(enumerate(srcp.header[1].items()))
    return hdrdict[39][1]/1000



def get_xy(fname,xhdr,yhdr,xyscalerhdr):
    xclst = list()
    yclst = list()
    with sg.open(fname,'r',ignore_geometry= True) as srcp:
        for trnum,tr in enumerate(srcp.trace):
            xysc = np.fabs(srcp.header[trnum][xyscalerhdr])
            xclst.append(srcp.header[trnum][xhdr]/ xysc)
            yclst.append(srcp.header[trnum][yhdr]/ xysc)
    return xclst,yclst

def process_sscalecols(seisdf,includexy=False):
    if includexy:
        seisdf['Xscaled'] = seisdf[seisdf.columns[0]]
        seisdf['Yscaled'] = seisdf[seisdf.columns[1]]
        # seisdf['Zscaled'] = seisdf[seisdf.columns[2]]
        # z is not added because it is a slice, i.e. constant z
    if seisdf.isnull().values.any():
        print('Warning: Null Values in the file will be dropped')
        seisdf.dropna(inplace=True)
    # xyzcols = [0,1,2]
    xyzcols = [0,1]
    xyz = seisdf[seisdf.columns[xyzcols]]
    seisdf.drop(seisdf.columns[[xyzcols]],axis=1,inplace=True)
    cols = seisdf.columns.tolist()
    seisdfs = StandardScaler().fit_transform(seisdf.values)

    seisdfsdf = pd.DataFrame(seisdfs,columns=cols)
    seisdfsxyz = pd.concat([xyz,seisdfsdf],axis=1)
    return seisdfsxyz

def process_seiswellattrib(sa,wa):
    # print(sa.head())
    # print(wa.head())
    xs = sa.iloc[:,0]
    ys = sa.iloc[:,1]
    xys = np.transpose(np.vstack((xs,ys)))
    xyhull = qhull(xys)

    xw = wa.iloc[:,2].values
    yw = wa.iloc[:,3].values
    # xyw = np.transpose(np.vstack((xw,yw)))
    wz = wa.iloc[:,1]  #z value
    wid = wa.iloc[:,0]  #wellname
    wida = wa.iloc[:,4] #well porosity or any other attribute

    ma = filterhullpolygon_mask(xw,yw,xyhull)
    # print('Remaining Wells after convex hull: ',len(ma))
    xwf= xw[ma]
    ywf = yw[ma]
    wzf = wz[ma]
    widf = wid[ma]
    wattribf = wida[ma]

    xywf = np.transpose(np.vstack((xwf,ywf)))
    welldfcols =['WELL','DEPTH','DEVX','DEVY']
    #wellsin_df = pd.DataFrame([widf,xwf,ywf,wattrib],columns=welldfcols)
    wellsin_df = pd.DataFrame(widf,columns = [welldfcols[0]])
    wellsin_df[welldfcols[1]] = wzf
    wellsin_df[welldfcols[2]] = xwf
    wellsin_df[welldfcols[3]] = ywf
    # print('wellsin df shape:',wellsin_df.shape)
    # print(wellsin_df.head(10))

    nattrib = sa.shape[1]
    # print('nattrib:',nattrib)
    #welldflist =[]
    for i in range(2,nattrib):
        vs = sa.iloc[:,i].values
        zwsa = map2ddata(xys,vs,xywf)
        # print('i:',i,'zwsa:',zwsa.size)
        #welldflist.append(zwsa)
        colname = sa.columns[i]
        wellsin_df[colname] = zwsa
    wa_col = wa.columns[4]
    wellsin_df[wa_col] = wattribf
    # print('Inside seiswellattrib ....')
    # print(wellsin_df.tail())
    # print('wellsin shape:',  wellsin_df.shape)
    wellsin_df.dropna(axis=0,inplace=True)
    # print(wellsin_df.tail())
    # print('wellsin shape after dropna:',  wellsin_df.shape)


    return wellsin_df



def plotwells(wdf,hideplots=True):
    wlst = wdf[wdf.columns[0]].unique().tolist()
    print(wlst)
    pdfcl = "AllWellsPsuedoLogs.pdf"
    with PdfPages(pdfcl) as pdf:
        for wname in wlst:
            wndf = wdf[wdf[wdf.columns[0]] == wname]
            fig,ax = plt.subplots(figsize = (4,7))
            ax.invert_yaxis()
            ax.plot(wndf[wndf.columns[-2]],wndf[wndf.columns[1]],label='Actual')
            ax.plot(wndf[wndf.columns[-1]],wndf[wndf.columns[1]],label ='Predicted')
            ax.set_xlabel(wndf.columns[-1])
            ax.set_ylabel(wndf.columns[1])
            ax.set_title(wname)
            ax.legend()
            fig.tight_layout()
            pdf.savefig()
            if not hideplots:
                fig.show()
            plt.close()

def  getcommandline():
    parser = argparse.ArgumentParser(description='ML to convert seismic to logs')
    parser.add_argument('segyfileslist',help='File that lists all segy files to be used for attributes')
    parser.add_argument('wellscsv',help='csv file with well depth devx devy log ')
    parser.add_argument('--segyxhdr',type=int,default=73,help='xcoord header.default=73')
    parser.add_argument('--segyyhdr',type=int,default=77,help='ycoord header. default=77')
    parser.add_argument('--xyscalerhdr',type=int,default=71,help='hdr of xy scaler to divide by.default=71')
    parser.add_argument('--startendslice',type=int,nargs=2,default=[1000,2000],help='Start end slice. default= 1000 2000')
    parser.add_argument('--cbriterations',type=int,default=500,help='Learning Iterations, default =500')
    parser.add_argument('--cbrlearningrate',type=float,default=0.03,help='learning_rate. default=0.03')
    parser.add_argument('--cbrdepth',type=int,default=2,help='depth of trees. default=2')
    parser.add_argument('--includexy',action='store_true',default=False,help='include x y coords in model.default= not to')
    parser.add_argument('--slicesout',action='store_true',default=False,
        help='Save individual unscaled slices to csv. default=false, i.e do not save')
    parser.add_argument('--plotincrement',type=int,default=100,
        help='increment for xplotting actual vs predicted. default= every 100th slice ')
    parser.add_argument('--outdir',help='output directory,default= same dir as input')
    parser.add_argument('--hideplots',action='store_true',default=False,
        help='Only save to pdf. default =show and save')

    result=parser.parse_args()
    if not result.segyfileslist:
        parser.print_help()
        exit()
    else:
        return result



def main():
    cmdl = getcommandline()
    if cmdl.wellscsv:
        allwells = pd.read_csv(cmdl.wellscsv)
        # dz = np.diff(allwells.DEPTH)[2]
        dz = np.diff(allwells[allwells.columns[1]])[2]
        print('Well Vertical increment {}'.format(dz))
        wdirsplit,wfextsplit= os.path.split(cmdl.wellscsv)
        wfname,wfextn= os.path.splitext(wfextsplit)
        # logname = allwells.columns[-1]
        wcols = allwells.columns.tolist()
        print(wcols)
        logname = wcols[-1]
        print('logname:',logname)
        lognamepred = logname + 'pred'
        wcols.append(lognamepred)



        if cmdl.outdir:
            outfw = os.path.join(cmdl.outdir,wfname) +"_pred.csv"
        else:
            outfw = os.path.join(wdirsplit,wfname) +"_pred.csv"


    if cmdl.segyfileslist:
        sflist = list()
        sflist = process_segylist(cmdl.segyfileslist)

        dirsplit,fextsplit= os.path.split(sflist[0])
        fname,fextn= os.path.splitext(fextsplit)
        if cmdl.outdir:
            outfsegy = os.path.join(cmdl.outdir,wfname) +"_p%s.sgy" % (logname)
        else:
            outfsegy = os.path.join(dirsplit,wfname) +"_p%s.sgy" %(logname)

        print('Copying file, please wait ........')
        start_copy = datetime.now()
        copyfile(sflist[0], outfsegy)
        end_copy = datetime.now()
        print('Duration of copying: {}'.format(end_copy - start_copy))

        sr = get_samplerate(outfsegy)
        print('Seismic Sample Rate: {}'.format(sr))

        print('Zeroing segy file, please wait ........')
        start_zero = datetime.now()
        zero_segy(outfsegy)
        end_zero = datetime.now()
        print('Duration of zeroing: {}'.format(end_zero - start_zero))


        xclst,yclst = get_xy(fextsplit,cmdl.segyxhdr,cmdl.segyyhdr,cmdl.xyscalerhdr)
        xydf = pd.DataFrame({'XC':xclst,'YC':yclst})
        preddf = xydf.copy()
        scols = list()
        for f in sflist:
            dirsplit,fextsplit= os.path.split(f)
            fname,fextn= os.path.splitext(fextsplit)
            scols.append(fname)

        sfname = 'allattrib'
        # slicerange = cmdl.startendslice[1] - cmdl.startendslice[0]
        sstart = int(cmdl.startendslice[0] // dz)
        send = int(cmdl.startendslice[1] // dz)
        start_process = datetime.now()
        slicelst   = list()
        slicenumlst =list()
        wnlst      = list()
        slicewnlst = list()
        coef0lst   = list()
        coef1lst   = list()
        r2lst      = list()
        for slicenum in range (sstart, send):
            if cmdl.outdir:
                outfslice = os.path.join(cmdl.outdir,sfname) +"_slice%d.csv"%slicenum
            else:
                outfslice = os.path.join(dirsplit,sfname) +"_slice%d.csv"%slicenum
            zslice = slicenum *dz
            wdf = allwells[allwells.DEPTH == zslice]
            c = wdf.columns[4] #log name
            nw = wdf[~ wdf[c].isnull() ].count()[4]
            print('# of wells for depth slice {} is {}'.format(zslice,nw))


            slicefiles = list()
            for i in range(len(sflist)):
                slicefiles.append(get_slice(sflist[i],slicenum))
            slicear = np.array(slicefiles).T
            slicedf = pd.DataFrame(slicear,columns=scols)

            alldata = pd.concat((xydf,slicedf),axis=1)
            print('Slice#: {} @ Time : {} ms'.format(slicenum,zslice) )
            # print(alldata.head())
            if cmdl.slicesout:
                alldata.to_csv(outfslice,index=False)
            alldatas = process_sscalecols(alldata,includexy=cmdl.includexy)
            # print('After Scaling .....')
            # print(alldatas.head())
            wdfsa = process_seiswellattrib(alldatas,wdf)
            print(wdfsa.tail())
            # lastcol = wdfsa.shape[1]
            X = wdfsa.iloc[:,4 : -1]
            y = wdfsa.iloc[:,-1]
            Xpred = alldatas.iloc[:,2:]
            # print(f'Xpred: {Xpred.shape}' )
            # print('# of wells used: ', X.shape[0], y.shape)
            # print(f'X shape: {X.shape} ')
            # print(X )

            model = CatBoostRegressor(iterations=cmdl.cbriterations, learning_rate= cmdl.cbrlearningrate,
                    depth=cmdl.cbrdepth,loss_function= 'RMSE',random_seed=42,logging_level='Silent')
            # Fit model
            model.fit(X, y)
            # Get predictions
            ypred = model.predict(X)
            # Calculating Mean Squared Error
            mse = np.mean((ypred - y)**2)
            print('Metrics on input data: ')
            print('MSE: %.4f' %(mse))
            r2 = r2_score(y,ypred)
            print('R2 : %10.3f' % r2)

            ccmdl = sts.pearsonr(y,ypred)
            if slicenum == sstart:
                wellsdf = wdfsa[wdfsa.columns[:4]].copy()
                wellsdf[logname] = wdfsa[wdfsa.columns[-1]].copy()
                wellsdf[lognamepred] = ypred
                # print(wellsdf.tail())
                # print(wellsdf.shape)
            else:
                wellsdf0 = wdfsa[wdfsa.columns[:4]].copy()
                wellsdf0[logname] = wdfsa[wdfsa.columns[-1]].copy()
                wellsdf0[lognamepred] = ypred
                allwellspred = wellsdf.append(wellsdf0)
                wellsdf = allwellspred[wcols].copy()
                print(allwellspred.tail())
                print(allwellspred.shape)

            pred = model.predict(Xpred)
            alldatas[wdfsa.columns[4]] = pred
            # print('After Prediction........')
            # print(alldatas.head())
            slicestr = '{:.0f}'.format(zslice)
            preddf[slicestr] = pred

            qc0 = np.polyfit(y,ypred,1)
            xrngmin,xrngmax = y.min(),y.max()
            xvi = np.linspace(xrngmin,xrngmax)
            yvi0 = np.polyval(qc0,xvi)


            if slicenum % cmdl.plotincrement == 0:
                slicedepth = slicenum * dz
                fig,ax = plt.subplots()

                plt.scatter(y,ypred,alpha=0.5,c='b',label='Model Predicted')
                plt.plot(xvi,yvi0,c='k',lw=2)



                ax.annotate('Model = %-.*f * Actual + %-.*f' %
                        (2,qc0[0],2,qc0[1]),xy =(xvi[0],yvi0[0]),xytext=(0.14,0.85),
                        textcoords='figure fraction', fontsize=10)
                ax.annotate ('Model Pearson cc = %-.*f   Pearson p = %-.*f' %
                        (2,ccmdl[0],3,ccmdl[1]),xy =(xvi[0],yvi0[0]),xytext=(0.14,0.81),
                        textcoords='figure fraction', fontsize=10)
                ax.set_title(f'CBR Slice {slicedepth:.0f} Pseudo {logname}' )
                ax.set_xlabel('Actual')
                ax.set_ylabel('Predicted')
                if not cmdl.hideplots:
                    plt.show()
                swfname = 'SWAttrib'
                if cmdl.outdir:
                    # pdfcl = os.path.join(cmdl.outdir,swfname)+"%d" %(slicenum) +"_cbr%s.pdf" %(logname)
                    # wsdf = os.path.join(cmdl.outdir,swfname)+"%d" %(slicenum) +"_cbr%s.csv" %(logname)
                    pdfcl = os.path.join(cmdl.outdir,swfname)+ f"{slicedepth:.0f}"  +"_cbr%s.pdf" %(logname)
                    wsdf = os.path.join(cmdl.outdir,swfname)+ f"{slicedepth:.0f}" +"_cbr%s.csv" %(logname)
                else:
                    pdfcl = os.path.join(dirsplit,swfname) + f"{slicedepth:.0f}" +"_cbr%s.pdf" %(logname)
                    wsdf = os.path.join(dirsplit,swfname)+ f"{slicedepth:.0f}" +"_cbr%s.csv" %(logname)
                fig.savefig(pdfcl)
                wdfsa.to_csv(wsdf,index=False)
                print(f'Successfully generated {wsdf}')


            slicelst.append(zslice)
            wnlst.append(nw)
            slicewnlst.append(wdfsa.shape[0])
            slicenumlst.append(slicenum)
            r2lst.append(r2)
            coef0lst.append(qc0[0])
            coef1lst.append(qc0[1])

        end_process = datetime.now()
        print('Duration of ML model building and prediction : {}'.format(end_process - start_process))

        qccols = ['SLICENUM','SLICEZ','WELLSFOUND','WELLSUSED','COEF0','COEF1','R2']
        qcdf = pd.DataFrame({'SLICENUM':slicenumlst,'SLICEZ':slicelst,'WELLSFOUND':wnlst,'WELLSUSED':slicewnlst,
            'COEF0':coef0lst,'COEF1' :coef1lst,'R2':r2lst})
        qcdf = qcdf[qccols].copy()

        if cmdl.outdir:
            outseispred = os.path.join(cmdl.outdir,wfname) +"_slices.csv"
            outqc = os.path.join(cmdl.outdir,wfname) +"_qc.csv"
        else:
            outseispred = os.path.join(dirsplit,wfname) +"_slices.csv"
            outqc = os.path.join(dirsplit,wfname) +"_qc.csv"

        preddf.to_csv(outseispred,index=False)
        print('Successfully generated {}'.format(outseispred))
        print('DataFrame size: ', preddf.shape)
        endsmpl = preddf.shape[1] - 2
        # print(preddf.head())

        qcdf.to_csv(outqc,index=False)
        print('Successfully generated {}'.format(outqc))



        with sg.open( outfsegy, "r+" ) as srcp:
            for trnum,tr in enumerate(srcp.trace):
                trplog = preddf.iloc[trnum,2:].values
                # lentrplog = trplog.size
                # print(trplog)
                tr[sstart: (sstart + endsmpl)] = trplog
                srcp.trace[trnum] = tr
        print('Successfully generated {}'.format(outfsegy))


        allwellspred.to_csv(outfw,index=False)
        print('Successfully generated {}'.format(outfw))
        plotwells(allwellspred,hideplots=cmdl.hideplots)





if __name__ == '__main__':
    main()
