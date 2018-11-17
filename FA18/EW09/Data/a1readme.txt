In these data files each column corresponds to a given variable, moving down the column 
goes forward in time. There are three types
statevariables: model evolved variables
radiances: proxy for observed sattelite radiances
radiances_error: radiences with gaussian white noise added with variance 
equal to half of the mean of the variable values. 

The number in the file name represents the value for the Fc (Co2/icebreakup/whatever) 
forcing paramater. 


The statevariable file colums correspond to:
column1=energy (E)
column2=maximum attainable albedo (am)
column3=actual albedo (a)
column4= ice concentration (Ci)
column5= pond fraction (Cp)

The radiance and radience_error columns correspond to

column1=|E*am|
column2= am-a
column3=|E*a|
column4=(.5+.4*tanh((-(E-50)/10)))*(E)+273.15)
column5=Cp*Ci

Non of them are 1-1!

All of the data generated here began with an inital condition of E=-250 am=0.75 

This is also the non-autonomous system so there is seasonality.  

