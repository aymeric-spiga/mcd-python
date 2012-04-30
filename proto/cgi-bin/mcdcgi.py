#!/home/aymeric/Software/epd-7.0-2-rh5-x86/bin/python
### here the version used to f2py the MCD Fortran routines

##################################################
### A Python CGI for the Mars Climate Database ###
### -------------------------------------------###
### Aymeric SPIGA 18-19/04/2012                ###
### -------------------------------------------###
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

# Get data from user-defined fields
query.lat = float(form.getvalue("latitude"))
query.lon = float(form.getvalue("longitude"))
query.loct = float(form.getvalue("localtime"))
query.xdate = float(form.getvalue("ls"))
query.xz = float(form.getvalue("altitude"))
query.hrkey = int(form.getvalue("hrkey"))
query.dust = int(form.getvalue("dust"))
#        self.zkey      = 3  # specify that xz is the altitude above surface (m)
#        self.perturkey = 0  #integer perturkey ! perturbation type (0: none)
#        self.seedin    = 1  #random number generator seed (unused if perturkey=0)
#        self.gwlength  = 0. #gravity Wave wavelength (unused if perturkey=0)

# Get free dimensions
islatfree  = float(form.getvalue("islatfree"))
islonfree  = float(form.getvalue("islonfree"))
isloctfree = float(form.getvalue("isloctfree"))
isaltfree  = float(form.getvalue("isaltfree"))
sumfree = islatfree + islonfree + isloctfree + isaltfree 
if sumfree > 2: exit() ## only 1D or 2D plots for the moment

# reference name (to test which figures are already in the database)
reference = str(islatfree)+str(islonfree)+str(isloctfree)+str(isaltfree)+query.getnameset()
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
    query.plot1d(["t","p","u","v"],vertplot=isaltfree)    
    mpl.savefig("img/temp.png",dpi=85,bbox_inches='tight',pad_inches=0.25)
    Image.open("img/temp.png").save(figname,'JPEG')

  ### 2D plots
  elif sumfree == 2:

    ### getting data
    if islatfree == 1 and islonfree == 1:  	query.latlon()   
    else:					exit()  
    
    query.map2d(["t","u"]) 
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
if sumfree == 0: 	query.update() ; query.printmeanvar()  
elif sumfree >= 1: 	print "<img src='../"+figname+"'><br />"
else:			exit()

## This is quite common
print "Based on the <a href='http://www-mars.lmd.jussieu.fr'>Mars Climate Database</a> (c) LMD/OU/IAA/ESA/CNES.<br />"
print "<hr>"
query.printset()
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
