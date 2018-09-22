

# MLSEISTOLOG
>  *__mlseistolog.py__* attempts to use various seismic segy attribute volumes to fit an ML model to predict logs at every trace in the volume.

##  CONCEPT
>  There are various sources of seismic attributes, e.g. post stack instantaneous attributes, post stack structural attributes, post stack frequency attributes, post stack spectral decomposition attributes, prestack AVO/QI attributes.  

>  Any or all of these attributes can be generated as seismic segy volume(s) which can be input to influence the computation of the machine learning (ML) model that is generated to predict a specific log.  

>Using *__preparelogs.py__* program, we generate a file with all the wells of a specific log. This file should have the same sampling interval as the segy, i.e. if your seismic volumes are in depth at 2 m interval, then the logs should have been generated with the same interval.  

>  The result will be a segy volume of the same dimensions as those of the input but data will be that of the log.  
Models are built one slice at a time, i.e. all sliced volumes are sliced at the same level and used to build an ML model using well log data. A prediction at all traces is performed and inserted into the segy.  

>  The current version uses CatBoostRegression to build the models.



## COMMAND LINE INTERFACE  

```
>python mlseistolog.py -h
usage: mlseistolog.py [-h] [--segyxhdr SEGYXHDR] [--segyyhdr SEGYYHDR]
                      [--xyscalerhdr XYSCALERHDR]
                      [--startendslice STARTENDSLICE STARTENDSLICE]
                      [--cbriterations CBRITERATIONS]
                      [--cbrlearningrate CBRLEARNINGRATE]
                      [--cbrdepth CBRDEPTH] [--includexy] [--slicesout]
                      [--plotincrement PLOTINCREMENT] [--outdir OUTDIR]
                      [--hideplots]
                      segyfileslist wellscsv

ML to convert seismic to logs

positional arguments:
  segyfileslist         File that lists all segy files to be used for
                        attributes
  wellscsv              csv file with well depth devx devy log

optional arguments:
  -h, --help            show this help message and exit
  --segyxhdr SEGYXHDR   xcoord header.default=73
  --segyyhdr SEGYYHDR   ycoord header. default=77
  --xyscalerhdr XYSCALERHDR
                        hdr of xy scaler to divide by.default=71
  --startendslice STARTENDSLICE STARTENDSLICE
                        Start end slice. default= 1000 2000
  --cbriterations CBRITERATIONS
                        Learning Iterations, default =500
  --cbrlearningrate CBRLEARNINGRATE
                        learning_rate. default=0.03
  --cbrdepth CBRDEPTH   depth of trees. default=2
  --includexy           include x y coords in model.default= not to
  --slicesout           Save individual unscaled slices to csv. default=false,
                        i.e do not save
  --plotincrement PLOTINCREMENT
                        increment for xplotting actual vs predicted. default=
                        every 100th slice
  --outdir OUTDIR       output directory,default= same dir as input
  --hideplots           Only save to pdf. default =show and save
  ```
  
>  There are 2 positional arguments that have to be supplied:  
>*  segy files list: this is a simple text file listing all the segy volumes that will be used as attributes  
>*  well log file: This is the output from *__preparelogs.py__*. This file has only one log for all the wells in your project.  
---

>  Optional arguments:  
>*  *--segyxhdr* to specify the header byte location for the x coordinate. default is 73  
>*  *--segyyhdr* to specify the header byte location for the y coordinate. default is 77  
---
>*  *--startendslice* 2 values are expected: the start slice in the vertical units (depth or time in ms) and the end slice. default values are 1000 2000. The resulting segy will have log values only within those limits. Other levels will be zeros.  
---  
>  CatBoostRegression Model parameters:  
>*  *--cbriterations*: ML model parameter. CatBoostRegression # of interations. default = 500  
>*  *--cbrdepth*: ML model parameter. CatBoostRegression depth of trees. default= 2  
>*  *--cbrlearningrate* : ML model parameter. CatBoostRegression learning rate. default = 0.03  
---  
>*  *--includexy* : include the x and y coordinates as attributes for the model building. Default is not to include the x y coordinates. By including them we assume that our log prediction is also based on spatial distribution. We can think of that as that of geostatistics without the restriction of variogram and only one co attribute to consider.  
---  
>*  *--sliceout* save every seismic slice out as a csv file. default is false, i.e. do not save those slices. The csv will be extremely large, because it will be the seismic data of all attributes but saved as csv. These slices are the ones used for prediction after the model is built.  
>*  *--plotincrement* : The increment at which a cross plot of actual vs predicted is generated for QC purposes. default is 100, i.e. every 100 slices generate a cross plot and the corresponding slice csv. These could be later used in *__swattriblist.py__* to further hyperparameter tune the model.  
---   
>*  *--outdir* to output all these files to a different directory other than the current working dir. Useful when you want to seperate results from input.  
>*  *--hideplots* use this option to only output to pdf and not have to click exit for every plot generated.  



## SUGGESTED WORKFLOW

*  Run *__mlseistolog.py__* on a small set of slices, e.g. 1500 to 1520 using the option *--slicestartend*. Also use option *--sliceout to create a csv of every slice.  
*  Using both csv's, the seismic and the predicted log as input for *__swattriblist.py__* to tune the model parameters  
*  Run *__mlseistolog.py__* again with option *--slicestartend* to cover the full range required. This time, do not use *--sliceout* option.  
*  The resulting segy should represent a log with comparable values at every trace. The segy should be standard and imported into the workstation for further evaluation.  



  
  
  