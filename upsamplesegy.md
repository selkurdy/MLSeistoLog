
# UPSAMPLING OF SEGY

## CONCEPT  
>   *  A segy file in time most likely is sampled at 2 ms, while in depth it is probably sampled at 2 m. The idea is to upsample the seismic to half ms or half a m. The upsampling is nothing but interpolation and adds no new information to the segy. The resulting segy file will be considerably bigger in size.  
>  *  The value of the upsampled segy will be useful when using ```mlseistolog.py``` program which processes the segy by slicing and fitting an ML model per slice to well log data. In using `` preparelogs.py `` a comparable interval to that of the upsampled segy has to chosen. 
>  *  After running `` mlseistolog.py`` all the samples will be predicted according to well log data and the resulting segy will be to a  higher frequency than the original seismic.  
>  *  Please note that this program only upscales, i.e. no downscaling  



## COMMAND LINE INTERFACE

```
>python upsamplesegy.py -h
usage: upsamplesegy.py [-h] [--outsampleinterval OUTSAMPLEINTERVAL]
                       [--endtime ENDTIME] [--outdir OUTDIR] [--hideplots]
                       segyfile

Upsample segy

positional arguments:
  segyfile              segy file to upsample

optional arguments:
  -h, --help            show this help message and exit
  --outsampleinterval OUTSAMPLEINTERVAL
                        output sample interval. default= 500, i.e. half of
                        segy file unit
  --endtime ENDTIME     Shorten the output traces to endtime in ms.
                        default=None, i.e. keep same length as input
  --outdir OUTDIR       output directory,default= same dir as input
  --hideplots           Only save to pdf. default =show and save
````  
 
>`` python upsamplesegy.py segyfilename.sgy`` will simply read the segy file and outputs the upsampled file   
---  
>``--outsampleinterval`` option to change the sample interval of the output  
>``--endtime``  is used when you want to shorten the traces in the vertical dimension. The default value is to keep the same length as the input. Please note that the start of the file is assumed to remain at zero   
>``--outdir`` option allows the user to output the segy to a different directory than the working one.  
>``--hideplots``  option hides plots generated at regular increments to QC traces. Currently, the default of traceplots is hardcoded at 50K. If ``--hideplots``  is used then user will have to click on each generated plot to exit for the program to continue.  
