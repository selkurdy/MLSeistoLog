# MLSeistoLog
Transform seismic segy to logs using Machine Learning

This is a three step workflow. 
###  Prepare logs:
Using *_preparelogs.py_* each well needs:
*  the las, or any flat file listing  depth and log values
*  the deviation survey file with MD X Y Z 
*  the time depth listing to convert the log to time

###  Average start and end times
Using *_logsrange.py_* compute the average start and end times for your data set


###  Transform seismic segy to logs
Using  *_mlseistolog.py_* submit the attributes segy files and the logsfile csv  
generated from *_preparelogs.py_*.  The result is a segy file of the same values as
the logs.

The approach is to slice the attributes segys and the logs file, create a data set then  
submit those to CatBoostRegression to compute a model per horizontal slice. The final product  
results in predicted log values at every trace location. The seismic is created one slice at a time
