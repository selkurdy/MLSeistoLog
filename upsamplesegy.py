"""
Program to upsample segy:
If segy is in time, i.e. sample interval is 2 ms it will be represented in the segy as
2000 microseconds. This will be upsampled to 500 microseconds, i.e. half a ms
If segy is in depth, i.e. sample interval is 3 m, it will be repreented in the segy as
3000 ?. This will be upsampled to 500, i.e. half a meter.
Samples has to be entered as integers.

The segy can be shortened by entering endtime in thousands, either ms or m as integers


"""
import os.path
import numpy as np
import argparse
import pandas as pd
import segyio
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



def upsample(segyfilein,segyfileout,
    newsi=500,
    endtout =None,
    outdir=None,
    hideplots=True):
    """
    upsample3 in notebook
    new sample interval is half a micro units,
    e.g. 2 ms is represented as 2000,
    2 m is represented as 2000
    0.5 m or ms is represented as 500 micro units,
    i.e. it is multiplied by 1000

    endtout is in ms
    """
    dirsplit,fextsplit= os.path.split(segyfilein)
    fname,fextn= os.path.splitext(fextsplit)



    with segyio.open(segyfilein,'r') as fi:
        spec = segyio.tools.metadata(fi)
        # allitems =list(fo.header[1].items())
        tsi = segyio.tools.dt(fi)
        print(f'tsi:{tsi}'  )
        nt = segyio.tools.sample_indexes(fi)
        # ns = len(spec.samples)
        # print(len(spec.samples))
        # tisource = spec.samples
        nta = np.array(nt)
        nta /= 1000.0
        endt = nta[-1]
        print(f'endt of original {endt} ')
        nsii = newsi
        if endtout:
            endt = endtout
        nsi = endt // (newsi / 1000 )
        ti = np.linspace(0,endt,int(nsi))
        # print(f'endt: {endt}, newsamples nsi: {nsi}')
        # print(f'len of ti: {ti.size}, last value of ti {ti[-1]}')
        # print(f'ti[50]: {ti[50]},ti[100]: {ti[100]}')
        # print(f'nta[50]: {nta[50]},nta[100]: {nta[100]}')
        # print(f'nsii: {nsii}')
        spec.samples = ti
        # print(spec)

        if outdir:
            pdfcl = os.path.join(outdir,fname)+ ".pdf"
        else:
            pdfcl = os.path.join(dirsplit,fname) +".pdf"

        with PdfPages(pdfcl) as pdf:

            with segyio.create(segyfileout, spec) as fo:
                fo.text[0] = fi.text[0]
                fo.bin = fi.bin
                #segyio.BinField.Samples = ti.size
                fo.header = fi.header
                fo.bin[3221] = ti.size
                for trnum,tr in enumerate(fi.trace):
                    ampi = interpolate.interp1d(nta,tr,bounds_error= False, fill_value =np.nan)
                    ai = ampi(ti)
                    if trnum % 50000 == 0 :
                        # print(ai.size)
                        fig,ax = plt.subplots(figsize=(5,10))
                        ax.invert_yaxis()
                        ax.plot(ai,ti,lw=3,label='interpolated')
                        ax.plot(tr,nta,label='actual')
                        ax.set_title(f'Tr#{trnum}')
                        ax.legend()

                        pdf.savefig()
                        if not hideplots:
                            fig.show()
                        plt.close()
                        # print(f'tr size: {tr.size}  ai size: {ai.size}')
                    fo.trace[trnum] = ai.astype('float32')
                    #fo.TraceField[115] = ti.size
                segyio.tools.resample(fo,rate=nsii,micro=True)


def getcommandline():
    parser = argparse.ArgumentParser(description='Upsample segy')
    parser.add_argument('segyfile', help='segy file to upsample')
    parser.add_argument('--outsampleinterval',type=int,default=500,
        help='output sample interval. default= 500, i.e. half of segy file unit')
    parser.add_argument('--endtime',type=int,default=None,
        help='Shorten the output traces to endtime in ms. default=None, i.e. keep same length as input')
    parser.add_argument('--outdir',help='output directory,default= same dir as input')
    parser.add_argument('--hideplots',action='store_true',default=False,
        help='Only save to pdf. default =show and save')

    result=parser.parse_args()
    if not result.segyfile:
        parser.print_help()
        exit()
    else:
        return result



def main():
    cmdl = getcommandline()
    dirsplit,fextsplit= os.path.split(cmdl.segyfile)
    fname,fextn= os.path.splitext(fextsplit)
    if cmdl.outdir:
        upsampledsegy = os.path.join(cmdl.outdir,fname) +"us%d.sgy" %cmdl.outsampleinterval
    else:
        upsampledsegy = os.path.join(dirsplit,fname) +"us%d.sgy" %cmdl.outsampleinterval



    upsample(cmdl.segyfile,upsampledsegy,
        newsi = cmdl.outsampleinterval,endtout=cmdl.endtime,
        outdir = cmdl.outdir,
        hideplots = cmdl.hideplots)


if __name__ == '__main__':
    main()
