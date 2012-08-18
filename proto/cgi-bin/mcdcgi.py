#!/home/aymeric/Software/epd-7.0-2-rh5-x86/bin/python
### here the version used to f2py the MCD Fortran routines

##################################################
### A Python CGI for the Mars Climate Database ###
### ------------------------------------------ ###
### Aymeric SPIGA 18-19/04/2012 ~ 11/08/2012   ###
### ------------------------------------------ ###
### (see mcdtest.py for examples of use)       ###
##################################################

import cgi, cgitb 
import numpy as np
from mcd import mcd
import cStringIO
import os as daos
import matplotlib.pyplot as mpl
import Image

# for debugging in web browser
cgitb.enable()

# Create instance of FieldStorage 
form = cgi.FieldStorage() 

# create a MCD object
query = mcd()

# Get data from user-defined fields and define free dimensions
getlat = form.getvalue("latitude")
if getlat == "all":  islatfree = 1 ; query.lat = -9999.
else:                islatfree = 0 ; query.lat = float(getlat)
getlon = form.getvalue("longitude")
if getlon == "all":  islonfree = 1 ; query.lon = -9999.
else:                islonfree = 0 ; query.lon = float(getlon)
getloct = form.getvalue("localtime")
if getloct == "all": isloctfree = 1 ; query.loct = -9999.
else:                isloctfree = 0 ; query.loct = float(getloct)
getalt = form.getvalue("altitude")
if getalt == "all":  isaltfree = 1 ; query.xz = -9999.
else:                isaltfree = 0 ; query.xz = float(getalt)
sumfree = islatfree + islonfree + isloctfree + isaltfree 
if sumfree > 2: exit() ## only 1D or 2D plots for the moment
query.xdate = float(form.getvalue("ls"))
query.hrkey = int(form.getvalue("hrkey"))
query.dust = int(form.getvalue("dust"))
#        self.zkey      = 3  # specify that xz is the altitude above surface (m)
#        self.perturkey = 0  #integer perturkey ! perturbation type (0: none)
#        self.seedin    = 1  #random number generator seed (unused if perturkey=0)
#        self.gwlength  = 0. #gravity Wave wavelength (unused if perturkey=0)

# Get variables to plot
var1 = form.getvalue("var1")
var2 = form.getvalue("var2")
var3 = form.getvalue("var3")
var4 = form.getvalue("var4")
vartoplot = [var1]
if var2 != "none": vartoplot = np.append(vartoplot,var2)
if var3 != "none": vartoplot = np.append(vartoplot,var3)
if var4 != "none": vartoplot = np.append(vartoplot,var4)
iswind = form.getvalue("iswind")
if iswind == "on": iswindlog = True
else:              iswindlog = False
isfixedlt = form.getvalue("isfixedlt")
if isfixedlt == "on": input_fixedlt=True
else:                 input_fixedlt=False  

# reference name (to test which figures are already in the database)
reference = str(islatfree)+str(islonfree)+str(isloctfree)+str(isaltfree)+query.getnameset()+str(var1)+str(var2)+str(var3)+str(var4)+str(iswind)+str(isfixedlt)
figname = 'img/'+reference+'.jpg'
testexist = daos.path.isfile(figname)

# extract data from MCD if needed
if not testexist:

  ### 1D plots
  if sumfree == 1:

    ### getting data
    if isloctfree == 1:  	query.diurnal(nd=24)
    elif islonfree == 1: 	query.zonal()
    elif islatfree == 1: 	query.meridional()
    elif isaltfree == 1: 	query.profile()   
    else:			exit()  

    ### generic building of figure
    #query.plot1d(["t","p","u","v"],vertplot=isaltfree)
    query.plot1d(vartoplot,vertplot=isaltfree)
    mpl.savefig("img/temp.png",dpi=85,bbox_inches='tight',pad_inches=0.25)
    Image.open("img/temp.png").save(figname,'JPEG')

  ### 2D plots
  elif sumfree == 2:

    ### getting data
    if islatfree == 1 and islonfree == 1:  	query.latlon()   
    else:					exit()  

    ### figure    
    query.map2d(vartoplot,incwind=iswindlog,fixedlt=input_fixedlt) 
    mpl.savefig("img/temp.png",dpi=110,bbox_inches='tight',pad_inches=0.4)
    Image.open("img/temp.png").save(figname,'JPEG') ##lighter images   
    ## http://www.pythonware.com/library/pil/handbook/introduction.htm

## This is quite common
print "Content-type:text/html\r\n\r\n"
print "<html>"
print "<head>"
print "<title>MCD. Simple Python interface</title>"
print "</head>"
print "<body>"

## Now the part which differs
if sumfree == 0: 	query.update() ; query.htmlprinttabextvar(vartoplot)  #query.printmeanvar()
elif sumfree >= 1: 	print "<img src='../"+figname+"'><br />"
else:			exit()

## This is quite common
#print "Based on the <a href='http://www-mars.lmd.jussieu.fr'>Mars Climate Database</a> (c) LMD/OU/IAA/ESA/CNES.<br />"
print "<hr>"
print "<a href='../index.html'>Click here to start a new query</a>."
#query.printset()
print "<hr>"
print "</body>"
print "</html>"

##write to file object
#f = cStringIO.StringIO()
#mpl.savefig(f)
#f.seek(0)

##output to browser
#print "Content-type: image/png\n"
#print f.read()
#exit()

#print "Content-type:text/html\r\n\r\n"
#print "<html>"
#print "<head>"
#print "<title>MCD. Simple Python interface</title>"
#print "</head>"
#print "<body>"






## HTTP Header
#print "Content-Type:application/octet-stream; name=\"FileName\"\r\n";
#print "Content-Disposition: attachment; filename=\"FileName\"\r\n\n";
## Actual File Content will go hear.
#fo = open("foo.txt", "rb")
#str = fo.read();
#print str
## Close opend file
#fo.close()
