# PREPARELOGS

## CONCEPT

*__preparelogs.py__* is designed to read per well the logs (las format or just plain ASCII flat files), its associated deviation file and format and compile all wells in one file with only one log curve. That file will be the input to *__mlseistolog__*  

The program takes the desired log for each well, e.g. GR, and generates a regular increment that you specify in the domain you desire. The regular increment should match that of the seismic you desire to run *__mlseistolog.py__* on. The domain could either be time or depth. 

The logs are smoothed by a filter to upscale them to the seismic traces. The filter length is an optional argument. The longer the filter the smoother the log. A QC plot of the actual log superimposed by the smoothed version is generated for every well.

We then create a batch file of all the logs available and process them to extract the desired log curve and generate a one big csv file for all the wells in your project.

## COMMAND LINE ARGUMENTS




```
python preparelogs.py -h
usage: preparelogs.py [-h] [--dtout DTOUT] [--logendtime LOGENDTIME]
                      [--dzout DZOUT] [--logenddepth LOGENDDEPTH] [--notlas]
                      [--logcolname LOGCOLNAME] [--logheader LOGHEADER]
                      [--logcols LOGCOLS LOGCOLS] [--lognull LOGNULL]
                      [--wellname WELLNAME] [--devheader DEVHEADER]
                      [--devcols DEVCOLS DEVCOLS DEVCOLS DEVCOLS]
                      [--devflipzsign] [--tdfilename TDFILENAME]
                      [--tdcols TDCOLS TDCOLS] [--tdheader TDHEADER]
                      [--flipsign] [--logsmoothwlen LOGSMOOTHWLEN]
                      [--hideplot]
                      logfilename devfilename

Process logs for mlseistolog

positional arguments:
  logfilename           Las file name
  devfilename           deviation file name

optional arguments:
  -h, --help            show this help message and exit
  --dtout DTOUT         time sampling interval in ms. default= 2
  --logendtime LOGENDTIME
                        Maximum T2W to output log.default=2500
  --dzout DZOUT         depth sampling interval in ms. default= 2
  --logenddepth LOGENDDEPTH
                        Maximum depth to output log.default=2500
  --notlas              Log file is not LAS format.default= las
  --logcolname LOGCOLNAME
                        log names. single word only. no default
  --logheader LOGHEADER
                        header lines of logs file.default=33
  --logcols LOGCOLS LOGCOLS
                        columns of log file to read. First col has to be
                        depth.default= 0 1
  --lognull LOGNULL     logs null value to remove.default= -999.25
  --wellname WELLNAME   Well name without spaces. Use if not las format
  --devheader DEVHEADER
                        Deviation survey header lines. default=15
  --devcols DEVCOLS DEVCOLS DEVCOLS DEVCOLS
                        Deviation Survey MD X Y Z columns.default = 0 1 2 3
  --devflipzsign        Deviation survey flip sign of z values.default: leave
                        as is
  --tdfilename TDFILENAME
                        time depth file name
  --tdcols TDCOLS TDCOLS
                        column #s of depth time 1W pairs.default=0 1
  --tdheader TDHEADER   fb file header lines to skip.default=1
  --flipsign            reverse sign of t d picks. default=keep as is
  --logsmoothwlen LOGSMOOTHWLEN
                        smooth logs window length. default=61. freq =
                        Vmax/Lambdamin
  --hideplot            Only save to pdf. default =show and save
  
  
  
  ```
  The program is expecting 2 positional arguments (i.e. must have):
  *  the log file name: which could be either a flat ASCII file or an LAS format
  *  the deviation file of the same well
  
  The optional arguments cover several situations:  
  
  *  You have to supply the argument *--logcolname*, e.g. GR, exactly as is listed in the log file. That necessitates checking what the code is for that curve in the file.
  *  If the log file is *las* format then that is all what is needed to define the log data
  *  If the log file is just an ASCII flat file then you have to define:
  >  * *--notlas* flag to signify that log file  is not las. default is las
  >  *  *--logheader* which is the number of header lines in the file  
  >  *  *--logcols*  2 values are expected, the first is the depth column # and the second is the desired log column #. Default values are 0 and 1
  >  *  *--lognull* null values. Default is -999.25 
  >  *  *--wellname* well name of this file. It should not have any spaces. If it has spaces you have to enclose it in double quotes  
  ---
  >  *  *--devheader* option to specify the number of header lines in the deviation survey file
  >  *  *--devcols* deviation survey file column numbers for MD, X, Y, and Z columns
  >  *  *--devflipzsign* flip sign of deviation z column (Petrel outputs Z in negative)
  ---
  >  *  *--logsmoothwlen* smoothing window length. the default is 61. The equation to figure it out is maximum Velocity / minimum wave length. This should give the frequencies expected.
  ---
  >  *--hideplot* is a good option to use especially if you are running a batch file because each well will generate a plot of both the original well and the smoothed version superimposed. The user needs to close that plot to continue. Using the hideplot option the plot is only saved as a pdf file.
  ---
  >  *--tdfilename* time depth pairs file name. It has to be supplied if you are outputting logs in the time domain.  
  >  *--tdcols* depth and time columns in the supplied file. default is 0 and 1. Time has to be one way in ms  
  ---
  >  *--dtout* if outputting logs in the time domain, then specify the vertical increment of the logs in time in ms  
  >  *--logendtime* deepest level in time in ms for sampling the log. Values beyond the actual log range will be filled vith NAN  
  ---  
  >  *--dzout* if outputting logs in the time domain, then specify the vertical increment of the logs in time in ms  
  >  *--logenddepth* deepest level in depth for sampling the log. Values beyond the actual log range will be filled vith NAN  
  >  Please note that you either request *--dtout* or *--dzout* , you cannot have both on the command line. The trigger to output to time is by supplying the option *--tdfilename*. If it is not given then the process will output logs in depth.  
  
  
  
  
  
  
  
  
  
  
  