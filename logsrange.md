# logsrange.py

## CONCEPT  
>  After running *__preparelogs.py__* on each and every log in the project, a file is created with all the logs, one after the other. To get ranges of start and end log data for all the wells, run *__logsrange.py__*.  
>A plot is generated giving all the starting and end depths with data for all the wells. Also after you close the chart, an average value of both the start and end depth is printed out.  These should be the ranges you would use for *__mlseistolog.py__*  


```
>python logsrange.py -h
usage: logsrange.py [-h] [--indepth] [--hideplot] logfilecsv

Find and plot depth ranges of all logs

positional arguments:
  logfilecsv  csv file with all wells generated from preparelogs.py

optional arguments:
  -h, --help  show this help message and exit
  --indepth   logfile is in depth. default=False,i.e. in time
  --hideplot  Hide plots.default=display
  
```   
  
## COMMAND LINE INTERFACE  
>  You need to supply the logs file csv  
---  
### Optional arguments  
``--indepth``  use it if data is processed in depth. default is log data csv is in time  
``--hideplot`` to not display the plot but only to save it as pdf. default is to both display and save the plot.  


  