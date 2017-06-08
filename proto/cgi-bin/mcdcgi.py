#!/usr/bin/python
###!/usr/bin/env python
###!/home/aymeric/Software/epd-7.0-2-rh5-x86/bin/python
####!/home/marshttp/EPD/epd-7.0-2-rh5-x86_64/bin/python
### here the version used to f2py the MCD Fortran routines

##################################################
### A Python CGI for the Mars Climate Database ###
### ------------------------------------------ ###
### Aymeric SPIGA 18-19/04/2012 ~ 11/08/2012   ###
### ------------------------------------------ ###
### (see mcdtest.py for examples of use)       ###
##################################################
### ajouts et corrections par Franck Guyon 09/2012
### ajouts suite a brainstorm equipe AS 10/2012

import cgi, cgitb 
import numpy as np
from modules import *

import cStringIO
import os as daos
import matplotlib.pyplot as mpl

import hashlib

### a function to read HTML arguments for coordinates
def gethtmlcoord(userinput,defmin,defmax):
   import string
   # accepted separators. the symbol - should always be last.
   separators = [":",";",",","/","_"," "] 
   # initial values
   val = -9999. ; vals = None ; vale = None ; foundinterv = False
   if userinput == None:   userinput = "1"
   # remove leading and trailing space in order to use space as separator
   userinput = userinput.strip()
   # the main work. either all, or an interval, or a single value.
   # ... all
   if userinput == "all":  isfree = 1 ; vals = defmin ; vale = defmax ; foundinterv = True
   # ... an interval
   else:
       for sep in separators:
         if not foundinterv:
           isfree = 1 ; ind = userinput.find(sep)
           if ind != -1: vals = float(userinput[:ind]) ; vale = float(userinput[ind+1:]) ; foundinterv = True
   # ... a single value (or an error)
   if not foundinterv: 
       # treat the difficult case of possible - separator
       test = userinput[1:].find("-") # at this stage:
                                      # * if - is found in the first position, it could not be a separator 
                                      # * if - found at positions > 0, it must be considered as a separator
       if test != -1: 
         isfree = 1 ; ind=test+1 
         vals = float(userinput[:ind]) ; vale = float(userinput[ind+1:]) ; foundinterv = True
       else:
         # check if input is valid (each character is numeric or -)
         for char in userinput: 
            if char not in string.digits: 
               if char not in ["-","."]: userinput="yorgl"
         # either we are OK. if we are not we set isfree to -1.
         if userinput != "yorgl":  isfree = 0 ; val = float(userinput)
         else:                     isfree = -1
   # return values
   return isfree, val, vals, vale

# set an errormess variable which must stay to None for interface to proceed
errormess = ""

# for debugging in web browser
cgitb.enable()

# Create instance of FieldStorage 
form = cgi.FieldStorage() 

# create a MCD object
# -- use in-code import to choose 
# --  between official and dev versions 
try:    dev = form.getvalue("dev")
except: dev = "off"
if dev == "on":
  from modules import mcd_dev
  query=mcd_dev.mcd()
else:
  from modules import mcd
  query=mcd.mcd()

# MCD version
query.toversion5(version="5.2")

## set MCD version changes if needed
#try:     betatest = form.getvalue("betatest")
#except:  betatest = "off"
#if betatest == "on": query.toversion5(version="5.2")
#else: query.toversion5(version="5.1")

# Get the kind of vertical coordinates and choose default behavior for "all"
try: query.zkey = int(form.getvalue("zkey"))
except: query.zkey = int(3)
if query.zkey == 2:    minxz = -5000.   ; maxxz = 150000.
elif query.zkey == 3:  minxz = 0.       ; maxxz = 250000.
elif query.zkey == 5:  minxz = -5000.   ; maxxz = 150000.
elif query.zkey == 4:  minxz = 1.e3     ; maxxz = 1.e-6
elif query.zkey == 1:  minxz = 3396000. ; maxxz = 3596000.

# Get data from user-defined fields and define free dimensions
getlat = form.getvalue("latitude")
if getlat is None: 
  errormess = errormess+"<li>A value or interval for latitude is missing. Please correct."
  islatfree = False ; query.lat = 0. ; query.lats = 0. ; query.late = 0.
else:
  islatfree,  query.lat,  query.lats,  query.late  = gethtmlcoord( getlat,   -90.,  90. )
islonfree,  query.lon,  query.lons,  query.lone  = gethtmlcoord( form.getvalue("longitude"), -180., 180. )
isloctfree, query.loct, query.locts, query.locte = gethtmlcoord( form.getvalue("localtime"),    0.,  24. )
isaltfree,  query.xz,   query.xzs,   query.xze   = gethtmlcoord( form.getvalue("altitude"),  minxz, maxxz)
if minxz < 0.1: minxz=0.1 # otherwise bug with values smaller than 0.1m

try: query.datekey = int(form.getvalue("datekeyhtml"))
except: query.datekey = float(1)
badlschar = False
if query.datekey == 1:
    try: query.xdate = float(form.getvalue("ls"))
    except: query.xdate = float(1) ; badlschar = True # comment the second part if in debug command line mode
else:
    try: query.xdate = float(form.getvalue("julian"))
    except: query.xdate = float(1)
    query.loct = 0.
try: query.dust = int(form.getvalue("dust"))
except: query.dust  = int(1)

# Prevent the user from doing bad
badinterv = (islatfree == -1) or (islonfree == -1) or (isloctfree == -1) or (isaltfree == -1)
if badinterv: 
    errormess = errormess+"<li>Bad syntax. Write a value (or) a range val1 val2 (or) 'all'. Separator shall be either ; : , / _ space"
badls = (query.datekey == 1 and (query.xdate < 0. or query.xdate > 360.))
if badls: 
    errormess = errormess+"<li>Solar longitude must be between 0 and 360."
if badlschar:
    errormess = errormess+"<li>Solar longitude is in the wrong format. It should be a positive number between 0 and 360. Intervals of solar longitude are not allowed."
badloct = (isloctfree == 0 and query.loct > 24.) \
       or (isloctfree == 1 and (query.locts > 24. or query.locte > 24.)) \
       or (isloctfree == 0 and query.loct < 0.) \
       or (isloctfree == 1 and (query.locts < 0. or query.locte < 0.))
if badloct: 
    errormess = errormess+"<li>Local time must be less than 24 martian hours (and not a negative number)."
badlat = (islatfree == 0 and abs(query.lat) > 90.) \
      or (islatfree == 1 and (abs(query.lats) > 90. or abs(query.late) > 90.))
if badlat: 
    errormess = errormess+"<li>Latitude coordinates must be between -90 and 90."
badlon = (islonfree == 0 and abs(query.lon) > 360.) \
      or (islonfree == 1 and (abs(query.lons) > 360. or abs(query.lone) > 360.))
if badlon: 
    errormess = errormess+"<li>Longitude coordinates must be between -360 and 360."
badalt = (isaltfree == 0 and (query.zkey in [3]) and query.xz < 0.) \
      or (isaltfree == 1 and (query.zkey in [3]) and (query.xzs < 0. or query.xze < 0.))
if badalt: 
    errormess = errormess+"<li>Vertical coordinates must be positive when requesting altitude above surface."
badalt2 = (isaltfree == 0 and (query.zkey in [1,4]) and query.xz <= 0.) \
      or (isaltfree == 1 and (query.zkey in [1,4]) and (query.xzs <= 0. or query.xze <= 0.))
if badalt2: 
    errormess = errormess+"<li>Vertical coordinates must be <b>strictly</b> positive when requesting pressure levels or altitude above Mars center."
badalt3 = (isaltfree == 0 and query.zkey == 4 and query.xz > 1500.) \
       or (isaltfree == 1 and query.zkey == 4 and min(query.xzs,query.xze) > 1500.)
if badalt3: 
    errormess = errormess+"<li>Pressure values larger than 1500 Pa are unlikely to be encountered in the Martian atmosphere."
badrange = (isloctfree == 1 and query.locts == query.locte) \
        or (islatfree == 1 and query.lats == query.late) \
        or (islonfree == 1 and query.lons == query.lone) \
        or (isaltfree == 1 and query.xzs == query.xze)
if badrange: 
    errormess = errormess+"<li>One or several coordinate intervals are not... intervals. Set either a real range or an unique value."
stormls = ( (query.dust in [4,5,6]) and (query.datekey == 1 and query.xdate < 180.))
if stormls:
    errormess = errormess+"<li>When a dust storm scenario is selected, available dates must be within the dust storm season (180 < Ls < 360)."
if query.xdate == 666.:
    errormess = "<li>CONGRATULATIONS! <br><img src='../surprise.jpg'><br> You reached secret mode.<br> You can <a href='http://www.youtube.com/watch?v=fTpQOZcNASw'>watch a nice video</a>."

# Get how many free dimensions we have
sumfree = islatfree + islonfree + isloctfree + isaltfree 
if sumfree >= 3: errormess = errormess + "<li>3 or more free dimensions are set... but only 1D and 2D plots are supported!"

# Get additional parameters
try: query.hrkey = int(form.getvalue("hrkey"))
except: query.hrkey = int(0)
#        self.perturkey = 0  #integer perturkey ! perturbation type (0: none)
#        self.seedin    = 1  #random number generator seed (unused if perturkey=0)
#        self.gwlength  = 0. #gravity Wave wavelength (unused if perturkey=0)
try: query.colorm = form.getvalue("colorm")
except: query.colorm = "jet"

try: query.min2d = float(form.getvalue("minval"))
except: query.min2d = None
try: query.max2d = float(form.getvalue("maxval"))
except: query.max2d = None

try: 
  proj = form.getvalue("proj")
  if proj == "": query.proj = None
  else: query.proj = proj
except:
  proj = "" 
  query.proj = None

try: query.trans = float(form.getvalue("trans"))/100.
except: query.trans = 0.0

try: query.plat = float(form.getvalue("plat"))
except: query.plat = 0.0
try: query.plon = float(form.getvalue("plon"))
except: query.plon = 0.0

try:
  query.latpoint = float(form.getvalue("latpoint"))
  query.lonpoint = float(form.getvalue("lonpoint"))
  strpoint = str(query.latpoint)+str(query.lonpoint)
except:
  query.latpoint = None
  query.lonpoint = None
  strpoint = ""

try: 
  query.dpi = form.getvalue("dpi")
  if query.dpi == "eps":  yeaheps = True  ; query.dpi = 300.
  else:                   yeaheps = False ; query.dpi = float(query.dpi)
except: 
  query.dpi = 80
  yeaheps = False

# Get variables to plot
var1 = form.getvalue("var1")
var2 = form.getvalue("var2")
var3 = form.getvalue("var3")
var4 = form.getvalue("var4")

# fg: init var as with form values
if var1 == None: var1="t"

vartoplot = []
if var1 != "none": vartoplot = np.append(vartoplot,var1)
if var2 != "none" and var2 != None: vartoplot = np.append(vartoplot,var2)
if var3 != "none" and var3 != None: vartoplot = np.append(vartoplot,var3)
if var4 != "none" and var4 != None: vartoplot = np.append(vartoplot,var4)

iswind = form.getvalue("iswind")
if iswind == "on": iswindlog = True
else:              iswindlog = False
isfixedlt = form.getvalue("isfixedlt")
if isfixedlt == "on": query.fixedlt=True
else:                 query.fixedlt=False  
iszonmean = form.getvalue("zonmean")
if iszonmean  == "on": query.zonmean=True
else:                  query.zonmean=False

islog = form.getvalue("islog")
if islog  == "on": query.islog=True
else:              query.islog=False

### now, proceed...
if errormess == "":

 # reference name (to test which figures are already in the database)
 try: reference = query.getnameset()+str(var1)+str(var2)+str(var3)+str(var4)+str(iswind)+str(isfixedlt)+str(iszonmean)+query.colorm+str(query.min2d)+str(query.max2d)+str(query.dpi)+str(islog)+str(proj)+str(query.trans)+str(query.plat)+str(query.plon)+strpoint
 except: reference = "test"
 if dev == "on": reference = 'dev_'+reference
 ## -- use a MD5 hash for a unique reference which avoids long names
 reference = hashlib.md5(reference).hexdigest()
 ##
 if yeaheps:  figname = '../img/'+reference+'.eps'
 else:        figname = '../img/'+reference+'.png'
 txtname = '../txt/'+reference+'.txt'
 testexist = daos.path.isfile(figname)

 # extract data from MCD if needed
 if not testexist:

  ### 1D plots
  if sumfree == 1:

    ### getting data
    if isloctfree == 1:  	query.diurnal(nd=25) 
    elif islonfree == 1: 	query.zonal(nd=64)
    elif islatfree == 1: 	query.meridional(nd=48)
    elif isaltfree == 1: 	query.profile(nd=35)   
    else:			exit()  

    ### generic building of figure
    query.htmlplot1d(vartoplot,figname=figname)

  ### 2D plots
  elif sumfree == 2:

    ### getting data
    if islatfree == 1 and islonfree == 1:     query.htmlmap2d(vartoplot,incwind=iswindlog,figname=figname) 
    else:                                     query.htmlplot2d(vartoplot,figname=figname)

  ### ASCII file outputs
  query.getascii(vartoplot,filename=txtname)

#### NOW WRITE THE HTML PAGE TO USER

## This is quite common
print "Content-type:text/html\n"
print "  "  #Apache needs a space after content-type

#entete="""<?xml version="1.0" encoding="UTF-8"?> 
#<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1 plus MathML 2.0//EN" "http://www.w3.org/Math/DTD/mathml2/xhtml-math11-f.dtd">
#<html xmlns="http://www.w3.org/1999/xhtml"> """

header="""<html><head><title>Mars Climate Database: The Web Interface</title></head><body>"""
#if betatest == "on": 
#    print "<b>!!! THIS IS A BETA VERSION. RESULTS ARE NOT VALIDATED !!!</b>"
#    if sumfree == 2:     print "<br>"

print header
#print query.printset()
#print "<br />"

## Now the part which differs
if errormess != "":
                       print "<h1>Ooops!</h1>"
                       print "Please correct the following problems before submitting again."
                       print "<ul>"
                       print errormess
                       print "</ul>"
else:
  if sumfree == 0: 	query.update() ; query.htmlprinttabextvar(vartoplot)
  elif sumfree >= 1:      
                        print "<a href='"+txtname+"'>Click here to download an ASCII file containing data</a><br />"
                        print "<hr>"
                        if yeaheps:  print "<hr><a href='"+figname+"'>!!!! Click here to download the EPS figure file !!!!</a><br /><hr>"
                        else:        print "<img src='"+figname+"'><br />"

## This is quite common
bottom = "</body></html>"
print bottom
