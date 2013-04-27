#! /usr/bin/env python

from mcd import mcd
from myplot import maplatlon
import numpy as np
import matplotlib.pyplot as mpl

from math import isnan

dafile = "/home/aymeric/Work/submitted/coauthorship/2012reiss/reiss_data.txt"

yorgl = np.loadtxt(dafile)
yorgl = np.transpose(yorgl)
ddwind = yorgl[5]
ddwinddir = yorgl[6]
ddwinderror = yorgl[9]

ddheight = yorgl[7]
ddheighterror = yorgl[8]

#ddheight = yorgl[7]*0. + 10.
#ddheighterror = yorgl[8]*0. + 0.5

stddev = np.std(ddwind[np.where(ddwind != -999.)])
print stddev

## missing points are NaN
ind = np.where(ddwind == -999.) ; ddwind[ind] = np.NaN
print "1. wind not available", yorgl[0][ind]
ind = np.where(ddheight == -999.) ; ddwind[ind] = np.NaN
print "1. height not available", yorgl[0][ind]
## no data point with rel error more than x%
#ind = np.where(100. * ddwinderror / ddwinddir > 5.)
#there was an error here!
tol = 50.
ind = np.where(100. * ddwinderror / ddwind > tol) ; ddwind[ind] = np.NaN 
print "2. error on wind too large"
print yorgl[0][ind]
print ddwinderror[ind]
ind = np.where(100. * ddheighterror / ddheight > tol) ; ddwind[ind] = np.NaN 
print "3. error on height too large"
print yorgl[0][ind]
print ddheighterror[ind]
## one record is doubtful according to Reiss
ind = np.where( yorgl[0] == 3 ) ; ddwind[ind] = np.NaN #; ddwinddir[ind] = np.NaN ; ddwinderror[ind] = np.NaN
print "4. doubtful", yorgl[0][ind]

#exit()

query = mcd()
mcdwind = [] ; mcdwindle = [] ; mcdwindhe = [] ; mcdangle = [] ; mcdanglele = [] ; mcdanglehe = [] 

for i in yorgl[0]:

   if not isnan(ddwind[i-1]):
    query.lat = yorgl[1][i-1]
    query.lon = yorgl[2][i-1]
    query.xdate = yorgl[3][i-1]
    query.loct = yorgl[4][i-1]
    query.xz = ddheight[i-1]
    #query.xz = 1000.  #ddheight is really the one that works best
    query.update() ; uu = query.zonwind ; vv = query.merwind
    if uu == -999. or vv == -999.: 
       speed = np.NaN ; speedmin = np.NaN ; speedmax = np.NaN
       angle = np.NaN ; anglemin = np.NaN ; anglemax = np.NaN
    else:
       speed = np.sqrt(uu*uu+vv*vv) ; speedmin = speed ; speedmax = speed

       angle = 90. - np.arctan2(vv,uu) * 180. / np.pi
       if angle < 0.: angle = angle + 360.

       query.xzs = max(ddheight[i-1] - 2.*ddheighterror[i-1],10.)
       query.xze = ddheight[i-1] + 2.*ddheighterror[i-1]
       query.profile()
       tab = np.sqrt(query.zonwindtab*query.zonwindtab + query.merwindtab*query.merwindtab)
       tab2 = 90. - np.arctan2(query.merwindtab,query.zonwindtab) * 180. / np.pi
       tab2[np.where(tab2 < 0.)] = tab2[np.where(tab2 < 0.)] + 360.
       speedmin = min(tab)
       speedmax = max(tab)
       anglemin = min(tab2)-5.
       anglemax = max(tab2)+5.
       if anglemin < 0.: anglemin = angle - (anglemax - angle)
       print int(i), speed, speedmin, speedmax 
       print int(i), angle, anglemin, anglemax
   else:
    #print "removed ",i,ddwinderror[i-1],ddheighterror[i-1]
    speed = -999. ; speedmin = -1000. ; speedmax = -998. ; angle = -999. ; anglemin = -1000. ; anglemax = -998.

   mcdwind = np.append(mcdwind,speed)
   mcdwindle = np.append(mcdwindle,speed-speedmin)
   mcdwindhe = np.append(mcdwindhe,speedmax-speed)
   mcdangle = np.append(mcdangle,angle)
   mcdanglele = np.append(mcdanglele,angle-anglemin)
   mcdanglehe = np.append(mcdanglehe,anglemax-angle)

mpl.figure(figsize=(12,10))
mpl.plot([0,30], [stddev,stddev+30], 'g--')
mpl.plot([0,30], [-stddev,-stddev+30], 'g--')
mpl.plot([0,30], [0,30], 'g-')

mpl.errorbar(mcdwind, ddwind, yerr=[ddwinderror,ddwinderror], xerr=[mcdwindle,mcdwindhe], fmt='bo')
mpl.xlim(xmin=0.,xmax=17.)
mpl.ylim(ymin=0.,ymax=30.)

mpl.xlabel('MCD horizontal wind speed (m/s)')
mpl.ylabel('Dust devil observed drift speed (m/s)')


for i in yorgl[0]:
    ## on multiplie stddev par 1.2 de facon a enlever les points tres au bord
    ##  qui sont in si on considere l error bar
    if abs(ddwind[i-1] - mcdwind[i-1]) > stddev*1.2:
        mpl.annotate(str(int(i)), xy=(mcdwind[i-1], ddwind[i-1]), xytext=(10,2), 
            textcoords='offset points', ha='center', va='bottom' )#,
            #bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3),
            #arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.3', 
            #color='red'))


##mpl.show()
mpl.savefig("comp.eps")
mpl.savefig("comp.png")
#mpl.savefig("comp4.eps")
#mpl.savefig("comp4.png")

ddwinddirerror = mcdanglele*0. + 15.

mpl.figure(figsize=(12,10))
mpl.errorbar(mcdangle,ddwinddir,yerr=[ddwinddirerror,ddwinddirerror],xerr=[mcdanglele,mcdanglehe],fmt='bo')
mpl.plot(mcdangle,ddwinddir,'bo')
mpl.xlim(xmin=0.,xmax=360.)
mpl.ylim(ymin=0.,ymax=360.)

mpl.xlabel('MCD horizontal wind orientation (deg)')
mpl.ylabel('Dust devil observed drift orientation (deg)')


stddev = 50.  ## voir e.g. 8 9 10
mpl.plot([0,360], [stddev,stddev+360], 'g--')
mpl.plot([0,360], [-stddev,-stddev+360], 'g--')
mpl.plot([0,360], [0,360], 'g-')


for i in yorgl[0]:
    ## on multiplie stddev par 1.2 de facon a enlever les points tres au bord
    ##  qui sont in si on considere l error bar
    if abs(ddwinddir[i-1] - mcdangle[i-1]) > stddev*1.2:
        mpl.annotate(str(int(i)), xy=(mcdangle[i-1], ddwinddir[i-1]), xytext=(10,2),
            textcoords='offset points', ha='center', va='bottom')#,
            #bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3),
            #arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.3',
            #color='red'))

##xytext=(18,30)


mpl.savefig("comp2.eps")
mpl.savefig("comp2.png")
#mpl.savefig("comp3.eps")
#mpl.savefig("comp3.png")


exit()


dd = mcd()





#dd.zkey = 4
#dd.xzs = 600.
#dd.xze = 0.1
#dd.profile()
#dd.plot1d(["t","p"])
#mpl.show()
#exit()

## CASE 1
dd.lat   = -14.6584
dd.lon   = 175.54
dd.xdate = 265.278792
dd.loct  = 14.8487
dd.xz    = 130.4

dd.update()
dd.printmcd()

dd.lats  = -5.  
dd.late  = -25.
dd.lons  = 160. 
dd.lone  = 190.
dd.map2d(["wind"],incwind=True,fixedlt=True)
## il y a un piege pour les cartes locales il faut fixedlt=True
mpl.savefig("case1.png",dpi=110,bbox_inches='tight',pad_inches=0.4)

## CASE 2
dd.lat   = -68.6125
dd.lon   = 11.4303
dd.xdate = 286.507293
dd.loct  = 15.129030
dd.xz    = 872.8
#dd.xz    = 10000.0

dd.update()
dd.printmcd()

dd.lats  = -75.
dd.late  = -60.
dd.lons  = 0.
dd.lone  = 20.
dd.map2d(["wind"],incwind=True,fixedlt=True)
mpl.savefig("case2.png",dpi=110,bbox_inches='tight',pad_inches=0.4)

##3	25,1306	314,4842	60,8946	14,9971	30,89575722	232	250,217961736111

## CASE 3
dd.lat   = 25.1306
dd.lon   = 314.4842
dd.xdate = 60.8946
dd.loct  = 14.9971
dd.xz    = 250.218

dd.update()
dd.printmcd()

##4	35,564	199,486	60,5883	14,9949	12,2972688101353	125	895,247012611329

## CASE 4
dd.lat   = 35.564
dd.lon   = 199.486
dd.xdate = 60.5883
dd.loct  = 14.9949
dd.xz    = 895.247

dd.update()
dd.printmcd()

##5	68,346	234,396	142,828	15,1777	13,4079899322222	83	1128,06581216829
##6	68,316	234,462	142,828	15,1819	16,1939071396022	85	582,624074786015
##7	68,302	234,144	142,828	15,1607	17,3354121022885	112	262,299040101764

## CASE 5
dd.lat   = 68.32
dd.lon   = 234.3
dd.xdate = 142.828
dd.loct  = 15.2

dd.xz    = 1128.066
dd.update()
dd.printmcd()

dd.xz    = 582.624
dd.update()
dd.printmcd()

dd.xz    = 262.299
dd.update()
dd.printmcd()


##19	27,055	129,614	13,5818	14,444	42,6488205391137	302	1395,96076100437

## CASE 19
dd.lat   = 27.055
dd.lon   = 129.614
dd.xdate = 13.5818
dd.loct  = 14.444

dd.xz    = 1395.961
dd.update()
dd.printmcd()

dd.lats  = 20.
dd.late  = 30.
dd.lons  = 125.
dd.lone  = 135.
dd.map2d(["wind"],incwind=True,fixedlt=True)
mpl.savefig("case19.png",dpi=110,bbox_inches='tight',pad_inches=0.4)


