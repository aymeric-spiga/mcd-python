####################################################
### A Python Class for the Mars Climate Database ###
### ---------------------------------------------###
### Aymeric SPIGA 17-21/04/2012                  ###
### ---------------------------------------------###
### (see mcdtest.py for examples of use)         ###
####################################################

import numpy as np
import fmcd
import matplotlib.pyplot as mpl
import myplot


class mcd():
 
    def __repr__(self):
    # print out a help string when help is invoked on the object
        whatprint = 'MCD object. \"help(mcd)\" for more information\n'
        return whatprint

########################
### Default settings ###
########################

    def __init__(self):
    # default settings
        ## 0. general stuff
        self.name      = "MCD v4.3"
        self.ack       = "Mars Climate Database (c) LMD/OU/IAA/ESA/CNES"
        #self.dset      = '/home/aymeric/Science/MCD_v4.3/data/'
        self.dset      = '/home/marshttp/MCD_v4.3/data/'
        ## 1. spatio-temporal coordinates
        self.lat       = 0.
        self.lats      = None
        self.late      = None
        self.lon       = 0.
        self.lons      = None
        self.lone      = None
        self.loct      = 0.
        self.locts     = None
        self.locte     = None
        self.xdate     = 0.  # see datekey
        self.xdates    = None
        self.xdatee    = None
        self.xz        = 10. # see zkey
        self.xzs       = None
        self.xze       = None
        ## 1bis. related settings
        self.zkey      = 3  # specify that xz is the altitude above surface (m)
                            # zkey  : <integer>   type of vertical coordinate xz
                            # 1 = radius from centre of planet (m)
                            # 2 = height above areoid (m) (MOLA zero datum)
                            # 3 = height above surface (m)
                            # 4 = pressure level (Pa)
                            # 5 = altitude above mean Mars Radius(=3396000m) (m)
        self.datekey   = 1  # 0 = "Earth time": xdate is given in Julian days (localtime must be set to zero)
                            # 1 = "Mars date": xdate is the value of Ls
        ## 2. climatological options
        self.dust      = 2  #our best guess MY24 scenario, with solar average conditions
        self.hrkey     = 1  #set high resolution mode on (hrkey=0 to set high resolution off)
        ## 3. additional settings for advanced use
        self.extvarkey = 1  #extra output variables (1: yes, 0: no)
        self.perturkey = 0  #integer perturkey ! perturbation type (0: none)
        self.seedin    = 1  #random number generator seed (unused if perturkey=0)
        self.gwlength  = 0. #gravity Wave wavelength (unused if perturkey=0)
        ## outputs. just to define attributes.
        ## --> in update
        self.pres = None ; self.dens = None ; self.temp = None ; self.zonwind = None ; self.merwind = None ; self.meanvar = None ; self.extvar = None
        self.seedout = None ; self.ierr = None
        ## --> in prepare
        self.xcoord = None ; self.ycoord = None
        self.prestab = None ; self.denstab = None ; self.temptab = None 
        self.zonwindtab = None ; self.merwindtab = None ; self.meanvartab = None ; self.extvartab = None
        ## plot stuff
        self.xlabel = None ; self.ylabel = None
        self.vertplot = False

    def viking1(self): self.name = "Viking 1 site. MCD v4.3 output" ; self.lat = 22.48 ; self.lon = -49.97 ; self.xdate = 97.
    def viking2(self): self.name = "Viking 2 site. MCD v4.3 output" ; self.lat = 47.97 ; self.lon = -225.74 ; self.xdate = 117.6

    def getdustlabel(self):
        if self.dust == 1: self.dustlabel = "MY24 minimum solar scenario"
        elif self.dust == 2: self.dustlabel = "MY24 average solar scenario"
        elif self.dust == 3: self.dustlabel = "MY24 maximum solar scenario"
        elif self.dust == 4: self.dustlabel = "dust storm minimum solar scenario"
        elif self.dust == 5: self.dustlabel = "dust storm average solar scenario"
        elif self.dust == 6: self.dustlabel = "dust storm maximum solar scenario"
        elif self.dust == 7: self.dustlabel = "warm scenario (dusty, maximum solar)"
        elif self.dust == 8: self.dustlabel = "cold scenario (low dust, minimum solar)"

    def gettitle(self):
        self.getdustlabel()
        self.title = self.name + " with " + self.dustlabel + "."
        if self.lats is None:  self.title = self.title + " Latitude " + str(self.lat) + "E"
        if self.lons is None:  self.title = self.title + " Longitude " + str(self.lon) + "N"
        if self.xzs is None:   
            self.vertunits()
            self.title = self.title + " Altitude " + str(self.xz) + " " + self.vunits
        if self.locts is None: self.title = self.title + " Local time " + str(self.loct) + "h"

    def getextvarlab(self,num):
        whichfield = { \
        91: "Pressure (Pa)", \
        92: "Density (kg/m3)", \
        93: "Temperature (K)", \
        94: "W-E wind component (m/s)", \
        95: "S-N wind component (m/s)", \
	1: "Radial distance from planet center (m)",\
	2: "Altitude above areoid (Mars geoid) (m)",\
	3: "Altitude above local surface (m)",\
	4: "orographic height (m) (surface altitude above areoid)",\
	5: "Ls, solar longitude of Mars (deg)",\
	6: "LST local true solar time (hrs)",\
	7: "Universal solar time (LST at lon=0) (hrs)",\
	8: "Air heat capacity Cp (J kg-1 K-1)",\
	9: "gamma=Cp/Cv Ratio of specific heats",\
	10: "density RMS day to day variations (kg/m^3)",\
        11: "[not defined]",\
        12: "[not defined]",\
	13: "scale height H(p) (m)",\
	14: "GCM orography (m)",\
	15: "surface temperature (K)",\
	16: "daily maximum mean surface temperature (K)",\
	17: "daily minimum mean surface temperature (K)",\
	18: "surf. temperature RMS day to day variations (K)",\
	19: "surface pressure (high resolution if hireskey=1)",\
	20: "GCM surface pressure (Pa)",\
	21: "atmospheric pressure RMS day to day variations (Pa)",\
	22: "surface pressure RMS day to day variations (Pa)",\
	23: "temperature RMS day to day variations (K)",\
	24: "zonal wind RMS day to day variations (m/s)",\
	25: "meridional wind RMS day to day variations (m/s)",\
	26: "vertical wind component (m/s) >0 when downwards!",\
	27: "vertical wind RMS day to day variations (m/s)",\
	28: "small scale perturbation (gravity wave) (kg/m^3)",\
	29: "q2: turbulent kinetic energy (m2/s2)",\
        30: "[not defined]",\
	31: "thermal IR flux to surface (W/m2)",\
	32: "solar flux to surface (W/m2)",\
	33: "thermal IR flux to space (W/m2)",\
	34: "solar flux reflected to space (W/m2)",\
	35: "surface CO2 ice layer (kg/m2)",\
	36: "DOD: Dust column visible optical depth",\
	37: "Dust mass mixing ratio (kg/kg)",\
	38: "DOD RMS day to day variations",\
	39: "DOD total standard deviation over season",\
	40: "Water vapor column (kg/m2)",\
	41: "Water vapor vol. mixing ratio (mol/mol)",\
	42: "Water ice column (kg/m2)",\
	43: "Water ice mixing ratio (mol/mol)",\
	44: "O3 ozone vol. mixing ratio (mol/mol)",\
	45: "[CO2] vol. mixing ratio (mol/mol)",\
	46: "[O] vol. mixing ratio (mol/mol)",\
	47: "[N2] vol. mixing ratio (mol/mol)",\
	48: "[CO] vol. mixing ratio (mol/mol)",\
	49: "R: Molecular gas constant (J K-1 kg-1)",\
	50: "Air viscosity estimation (N s m-2)"
        }
        if num not in whichfield: myplot.errormess("Incorrect subscript in extvar.")
        return whichfield[num]

    def convertlab(self,num):        
        ## a conversion from text inquiries to extvar numbers. to be completed.
        if num == "p": num = 91
        elif num == "rho": num = 92
        elif num == "t": num = 93
        elif num == "u": num = 94
        elif num == "v": num = 95
        elif num == "tsurf": num = 15
        elif num == "topo": num = 4
        elif num == "h": num = 13
        elif num == "ps": num = 19
        elif num == "tau": num = 36
        elif num == "mtot": num = 40
        elif num == "icetot": num = 42
        elif num == "ps_ddv": num = 22
        elif num == "h2ovap": num = 41
        elif num == "h2oice": num = 43
        elif num == "cp": num = 8
        elif num == "rho_ddv": num = 10
        elif num == "tsurfmx": num = 16
        elif num == "tsurfmn": num = 17
        elif num == "lwdown": num = 31
        elif num == "swdown": num = 32
        elif num == "lwup": num = 33
        elif num == "swup": num = 34
        elif num == "o3": num = 44
        elif num == "o": num = 46
        elif num == "co": num = 48
        elif num == "visc": num = 50
        elif num == "co2ice": num = 35
        elif not isinstance(num, np.int): myplot.errormess("field reference not found.")
        return num

###################
### One request ###
###################

    def update(self):
    # retrieve fields from MCD (call_mcd). more info in fmcd.call_mcd.__doc__
        (self.pres, self.dens, self.temp, self.zonwind, self.merwind, \
         self.meanvar, self.extvar, self.seedout, self.ierr) \
         = \
         fmcd.call_mcd(self.zkey,self.xz,self.lon,self.lat,self.hrkey, \
             self.datekey,self.xdate,self.loct,self.dset,self.dust, \
             self.perturkey,self.seedin,self.gwlength,self.extvarkey )
        ## we use the end of extvar (unused) to store meanvar. this is convenient for getextvar(lab)
        self.extvar[90] = self.pres ; self.extvar[91] = self.dens
        self.extvar[92] = self.temp ; self.extvar[93] = self.zonwind ; self.extvar[94] = self.merwind
        ## treat missing values 
        if self.temp == -999: self.extvar[:] = np.NaN ; self.meanvar[:] = np.NaN

    def printset(self):
    # print main settings
        print "zkey",self.zkey,"xz",self.xz,"lon",self.lon,"lat",self.lat,"hrkey",self.hrkey, \
              "xdate",self.xdate,"loct",self.loct,"dust",self.dust

    def getnameset(self):
    # set a name referring to settings [convenient for databases]
        strlat = str(self.lat)+str(self.lats)+str(self.late)
        strlon = str(self.lon)+str(self.lons)+str(self.lone)
        strxz = str(self.xz)+str(self.xzs)+str(self.xze)
        strloct = str(self.loct)+str(self.locts)+str(self.locte)
        name = str(self.zkey)+strxz+strlon+strlat+str(self.hrkey)+str(self.datekey)+str(self.xdate)+strloct+str(self.dust)
        return name

    def printcoord(self):
    # print requested space-time coordinates
        print "LAT",self.lat,"LON",self.lon,"LOCT",self.loct,"XDATE",self.xdate

    def printmeanvar(self):
    # print mean MCD variables
        print "Pressure = %5.3f pascals. " % (self.pres)
        print "Density = %5.3f kilograms per cubic meter. " % (self.dens)
        print "Temperature = %3.0f kelvins (%4.0f degrees celsius)." % (self.temp,self.temp-273.15)
        print "Zonal wind = %5.3f meters per second." % (self.zonwind)
        print "Meridional wind = %5.3f meters per second." % (self.merwind)
        print "Total horizontal wind = %5.3f meters per second." % ( np.sqrt(self.zonwind**2 + self.merwind**2) )

    def printextvar(self,num):
    # print extra MCD variables
        num = self.convertlab(num)
        print self.getextvarlab(num) + " ..... " + str(self.extvar[num-1])

    def printallextvar(self):
    # print all extra MCD variables    
        for i in range(50): self.printextvar(i+1)

    def htmlprinttabextvar(self,tabtodo):
        print "Results from the Mars Climate Database"
        print "<ul>"
        for i in range(len(tabtodo)): print "<li>" ; self.printextvar(tabtodo[i]) ; print "</li>"
        print "</ul>"
        print "<hr>"
        print "SETTINGS<br />"
        self.printcoord()
        self.printset()

    def printmcd(self):
    # 1. call MCD 2. print settings 3. print mean vars
        self.update()
        self.printcoord()
        print "-------------------------------------------"
        self.printmeanvar()

########################
### Several requests ###
########################

    def prepare(self,ndx=None,ndy=None):
    ### prepare I/O arrays for 1d slices
      if ndx is None:  print "No dimension in prepare. Exit. Set at least ndx." ; exit()
      else:            self.xcoord = np.ones(ndx)
      if ndy is None:  dashape = (ndx)     ; dashapemean = (ndx,6)     ; dashapeext = (ndx,101)     ; self.ycoord = None
      else:            dashape = (ndx,ndy) ; dashapemean = (ndx,ndy,6) ; dashapeext = (ndx,ndy,101) ; self.ycoord = np.ones(ndy)
      self.prestab = np.ones(dashape) ; self.denstab = np.ones(dashape) ; self.temptab = np.ones(dashape)
      self.zonwindtab = np.ones(dashape) ; self.merwindtab = np.ones(dashape) 
      self.meanvartab = np.ones(dashapemean) ; self.extvartab = np.ones(dashapeext)

    def getextvar(self,num):
    ### get a given var in extvartab
      try: field=self.extvartab[:,:,num] 
      except: field=self.extvartab[:,num]
      return field

    def definefield(self,choice):
    ### for analysis or plot purposes, set field and field label from user-defined choice
      choice = self.convertlab(choice)
      field = self.getextvar(choice); fieldlab = self.getextvarlab(choice)
      return field,fieldlab

    def ininterv(self,dstart,dend,nd,start=None,end=None,yaxis=False,vertcoord=False):
    ### user-defined start and end are used to create xcoord (or ycoord) vector
      if start is not None and end is not None:  first, second = self.correctbounds(start,end,vertcoord)
      else:                                      first, second = self.correctbounds(dstart,dend,vertcoord)  
      if self.zkey != 4 or not vertcoord:   tabtab = np.linspace(first,second,nd)
      else:                                 tabtab = np.logspace(first,second,nd)
      if not yaxis:      self.xcoord = tabtab
      else:              self.ycoord = tabtab

    def correctbounds(self,start,end,vertcoord):
      if self.zkey != 4 or not vertcoord:
        # regular altitudes
        if start > end: first = end ; second = start
        else:           first = start ; second = end
      else:
        # pressure: reversed avis
        if start < end: first = np.log10(end) ; second = np.log10(start)
        else:           first = np.log10(start) ; second = np.log10(end)
      return first, second

    def vertlabel(self):
      if self.zkey == 1:   self.xlabel = "radius from centre of planet (m)"
      elif self.zkey == 2: self.xlabel = "height above areoid (m) (MOLA zero datum)"
      elif self.zkey == 3: self.xlabel = "height above surface (m)"
      elif self.zkey == 4: self.xlabel = "pressure level (Pa)"
      elif self.zkey == 5: self.xlabel = "altitude above mean Mars Radius(=3396000m) (m)"

    def vertunits(self):
      if self.zkey == 1:   self.vunits = "m CP"
      elif self.zkey == 2: self.vunits = "m AMR"
      elif self.zkey == 3: self.vunits = "m ALS"
      elif self.zkey == 4: self.vunits = "Pa"
      elif self.zkey == 5: self.vunits = "m AMMRad"

    def vertaxis(self,number,yaxis=False):
      if self.zkey == 2:   self.ininterv(-5000.,100000.,number,start=self.xzs,end=self.xze,yaxis=yaxis,vertcoord=True)
      elif self.zkey == 3: self.ininterv(0.,120000.,number,start=self.xzs,end=self.xze,yaxis=yaxis,vertcoord=True)
      elif self.zkey == 5: self.ininterv(-5000.,100000.,number,start=self.xzs,end=self.xze,yaxis=yaxis,vertcoord=True)
      elif self.zkey == 4: self.ininterv(1000.,0.001,number,start=self.xzs,end=self.xze,yaxis=yaxis,vertcoord=True)
      elif self.zkey == 1: self.ininterv(3396000,3596000,number,start=self.xzs,end=self.xze,yaxis=yaxis,vertcoord=True)

###################
### 1D analysis ###
###################

    def put1d(self,i):
    ## fill in subscript i in output arrays
    ## (arrays must have been correctly defined through prepare)
      if self.prestab is None:  myplot.errormess("arrays must be prepared first through self.prepare")
      self.prestab[i] = self.pres ; self.denstab[i] = self.dens ; self.temptab[i] = self.temp
      self.zonwindtab[i] = self.zonwind ; self.merwindtab[i] = self.merwind
      self.meanvartab[i,1:5] = self.meanvar[0:4]  ## note: var numbering according to MCD manual is kept
      self.extvartab[i,1:100] = self.extvar[0:99] ## note: var numbering according to MCD manual is kept

    def diurnal(self,nd=13):
    ### retrieve a local time slice
      save = self.loct
      self.xlabel = "Local time (Martian hour)"
      self.prepare(ndx=nd) ; self.ininterv(0.,24.,nd,start=self.locts,end=self.locte) 
      for i in range(nd): self.loct = self.xcoord[i] ; self.update() ; self.put1d(i)
      self.loct = save

    def zonal(self,nd=37):
    ### retrieve a longitude slice
      save = self.lon
      self.xlabel = "East longitude (degrees)"
      self.prepare(ndx=nd) ; self.ininterv(-180.,180.,nd,start=self.lons,end=self.lone)
      for i in range(nd): self.lon = self.xcoord[i] ; self.update() ; self.put1d(i)
      self.lon = save

    def meridional(self,nd=19):
    ### retrieve a latitude slice
      save = self.lat
      self.xlabel = "North latitude (degrees)"
      self.prepare(ndx=nd) ; self.ininterv(-90.,90.,nd,start=self.lats,end=self.late)
      for i in range(nd): self.lat = self.xcoord[i] ; self.update() ; self.put1d(i)
      self.lat = save

    def profile(self,nd=20,tabperso=None):
    ### retrieve an altitude slice (profile)
      save = self.xz
      self.vertlabel()
      self.vertplot = True
      if tabperso is not None: nd = len(tabperso)
      correct = False
      self.prepare(ndx=nd) ; self.vertaxis(nd)
      if tabperso is not None: self.xcoord = tabperso
      for i in range(nd): self.xz = self.xcoord[i] ; self.update() ; self.put1d(i)
      self.xz = save

    def seasonal(self,nd=12):
    ### retrieve a seasonal slice
      save = self.xdate
      self.xlabel = "Areocentric longitude (degrees)"
      self.prepare(ndx=nd) ; self.ininterv(0.,360.,nd,start=self.xdates,end=self.xdatee)
      for i in range(nd): self.xdate = self.xcoord[i] ; self.update() ; self.put1d(i)
      self.xdate = save

    def makeplot1d(self,choice):
    ### one 1D plot is created for the user-defined variable in choice. 
      (field, fieldlab) = self.definefield(choice)
      if not self.vertplot:  absc = self.xcoord ; ordo = field ; ordolab = fieldlab ; absclab = self.xlabel
      else:                  ordo = self.xcoord ; absc = field ; absclab = fieldlab ; ordolab = self.xlabel
      mpl.plot(absc,ordo,'-bo') ; mpl.ylabel(ordolab) ; mpl.xlabel(absclab) #; mpl.xticks(query.xcoord)
      if self.zkey == 4: mpl.semilogy() ; ax = mpl.gca() ; ax.set_ylim(ax.get_ylim()[::-1])
      mpl.figtext(0.5, 0.01, self.ack, ha='center')

    def plot1d(self,tabtodo):
    ### complete 1D figure with possible multiplots
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      fig = mpl.figure() ; subv,subh = myplot.definesubplot( len(tabtodo) , fig ) 
      for i in range(len(tabtodo)): mpl.subplot(subv,subh,i+1).grid(True, linestyle=':', color='grey') ; self.makeplot1d(tabtodo[i])

    def htmlplot1d(self,tabtodo,figname="temp.png",title=""):
    ### complete 1D figure with possible multiplots
    ### added in 09/2012 for online MCD
    ### see http://www.dalkescientific.com/writings/diary/archive/2005/04/23/matplotlib_without_gui.html
      from matplotlib.figure import Figure
      from matplotlib.backends.backend_agg import FigureCanvasAgg
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      fig = Figure(figsize=(8,8)) ; subv,subh = myplot.definesubplot( len(tabtodo) , fig )
      for i in range(len(tabtodo)):
        yeah = fig.add_subplot(subv,subh,i+1) #.grid(True, linestyle=':', color='grey') 
        choice = tabtodo[i]
        (field, fieldlab) = self.definefield(choice)
        if not self.vertplot:  absc = self.xcoord ; ordo = field ; ordolab = fieldlab ; absclab = self.xlabel
        else:                  ordo = self.xcoord ; absc = field ; absclab = fieldlab ; ordolab = self.xlabel
        yeah.plot(absc,ordo,'-bo') #; mpl.xticks(query.xcoord)
        ax = fig.gca() ; ax.set_ylabel(ordolab) ; ax.set_xlabel(absclab)
        if self.zkey == 4: ax.set_yscale('log') ; ax.set_ylim(ax.get_ylim()[::-1])
      self.gettitle()
      fig.text(0.5, 0.95, self.title, ha='center')
      fig.text(0.5, 0.01, self.ack, ha='center')
      canvas = FigureCanvasAgg(fig)
      # The size * the dpi gives the final image size
      #   a4"x4" image * 80 dpi ==> 320x320 pixel image
      canvas.print_figure(figname, dpi=80)

###################
### 2D analysis ###
###################

    def latlon(self,ndx=37,ndy=19,fixedlt=False):
    ### retrieve a latitude/longitude slice
    ### default is: local time is not fixed. user-defined local time is at longitude 0.
      save1 = self.lon ; save2 = self.lat ; save3 = self.loct
      self.xlabel = "East longitude (degrees)" ; self.ylabel = "North latitude (degrees)"
      self.prepare(ndx=ndx,ndy=ndy)
      self.ininterv(-180.,180.,ndx,start=self.lons,end=self.lone)
      self.ininterv(-90.,  90.,ndy,start=self.lats,end=self.late,yaxis=True)
      if not fixedlt: umst = self.loct
      for i in range(ndx):
       for j in range(ndy):
         self.lon = self.xcoord[i] ; self.lat = self.ycoord[j]
         if not fixedlt: 
           if self.lons is not None and self.lone is not None: self.loct = (umst + (self.lons+self.lone)/30.) % 24
           else:                                               self.loct = (umst + self.lon/15.) % 24
         self.update() ; self.put2d(i,j)
      if not fixedlt: self.loct = umst
      self.lon = save1 ; self.lat = save2 ; self.loct = save3

    def lonalt(self,ndx=37,ndy=20,fixedlt=False):
    ### retrieve a longitude/altitude slice
      save1 = self.lon ; save2 = self.xz ; save3 = self.loct
      self.vertlabel() ; self.ylabel = self.xlabel
      self.xlabel = "East longitude (degrees)"
      self.prepare(ndx=ndx,ndy=ndy)
      self.ininterv(-180.,180.,ndx,start=self.lons,end=self.lone)
      self.vertaxis(ndy,yaxis=True)
      if not fixedlt: umst = self.loct
      for i in range(ndx):
       for j in range(ndy):
         self.lon = self.xcoord[i] ; self.xz = self.ycoord[j] 
         if not fixedlt: 
           if self.lons is not None and self.lone is not None: self.loct = (umst + (self.lons+self.lone)/30.) % 24
           else:                                               self.loct = (umst + self.lon/15.) % 24
         self.update() ; self.put2d(i,j)
      if not fixedlt: self.loct = umst
      self.lon = save1 ; self.xz = save2 ; self.loct = save3

    def latalt(self,ndx=19,ndy=20,fixedlt=False):
    ### retrieve a latitude/altitude slice
      save1 = self.lat ; save2 = self.xz ; save3 = self.loct
      self.vertlabel() ; self.ylabel = self.xlabel
      self.xlabel = "North latitude (degrees)"
      self.prepare(ndx=ndx,ndy=ndy)
      self.ininterv(-90.,90.,ndx,start=self.lats,end=self.late)
      self.vertaxis(ndy,yaxis=True)
      if not fixedlt: umst = self.loct
      for i in range(ndx):
       for j in range(ndy):
         self.lat = self.xcoord[i] ; self.xz = self.ycoord[j] 
         if not fixedlt: self.loct = (umst + self.lon/15.) % 24
         self.update() ; self.put2d(i,j)
      if not fixedlt: self.loct = umst
      self.lat = save1 ; self.xz = save2 ; self.loct = save3

    def put2d(self,i,j):
    ## fill in subscript i,j in output arrays
    ## (arrays must have been correctly defined through prepare)
      if self.prestab is None:  myplot.errormess("arrays must be prepared first through self.prepare")
      self.prestab[i,j] = self.pres ; self.denstab[i,j] = self.dens ; self.temptab[i,j] = self.temp
      self.zonwindtab[i,j] = self.zonwind ; self.merwindtab[i,j] = self.merwind
      self.meanvartab[i,j,1:5] = self.meanvar[0:4]  ## note: var numbering according to MCD manual is kept
      self.extvartab[i,j,1:100] = self.extvar[0:99] ## note: var numbering according to MCD manual is kept

    def makemap2d(self,choice,incwind=False,fixedlt=False,proj="cyl"):
    ### one 2D map is created for the user-defined variable in choice.
      self.latlon(fixedlt=fixedlt) ## a map is implicitely a lat-lon plot. otherwise it is a plot (cf. makeplot2d)
      if choice == "wind" or incwind:
          (windx, fieldlabwx) = self.definefield("u")
          (windy, fieldlabwy) = self.definefield("v")
      if choice == "wind":
          field = np.sqrt(windx*windx + windy*windy)
          fieldlab = "Horizontal wind speed (m/s)"
      else:    
          (field, fieldlab) = self.definefield(choice)
      if incwind:   myplot.maplatlon(self.xcoord,self.ycoord,field,title=fieldlab,proj=proj,vecx=windx,vecy=windy) #,stride=1)
      else:         myplot.maplatlon(self.xcoord,self.ycoord,field,title=fieldlab,proj=proj)
      mpl.figtext(0.5, 0.0, self.ack, ha='center')

    def map2d(self,tabtodo,incwind=False,fixedlt=False,proj="cyl"):
    ### complete 2D figure with possible multiplots
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      fig = mpl.figure()
      subv,subh = myplot.definesubplot( len(tabtodo) , fig ) 
      for i in range(len(tabtodo)): mpl.subplot(subv,subh,i+1) ; self.makemap2d(tabtodo[i],incwind=incwind,fixedlt=fixedlt,proj=proj)

    def htmlmap2d(self,tabtodo,incwind=False,fixedlt=False,figname="temp.png",title=""):
    ### complete 2D figure with possible multiplots
    ### added in 09/2012 for online MCD
    ### see http://www.dalkescientific.com/writings/diary/archive/2005/04/23/matplotlib_without_gui.html
      from matplotlib.figure import Figure
      from matplotlib.backends.backend_agg import FigureCanvasAgg
      from matplotlib.cm import get_cmap
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      fig = Figure(figsize=(8,8)) ; subv,subh = myplot.definesubplot( len(tabtodo) , fig )

      ### topocontours
      fieldc = self.getextvar(self.convertlab("topo"))

      for i in range(len(tabtodo)):
        yeah = fig.add_subplot(subv,subh,i+1)
        choice = tabtodo[i]
        self.latlon(fixedlt=fixedlt) 
        ## a map is implicitely a lat-lon plot. otherwise it is a plot (cf. makeplot2d)
        (field, fieldlab) = self.definefield(choice)
        if incwind: (windx, fieldlabwx) = self.definefield("u") ; (windy, fieldlabwy) = self.definefield("v")

        proj="cyl" ; colorb="jet" ; ndiv=20 ; zeback="molabw" ; trans=1.0 #0.6
        title="" ; vecx=None ; vecy=None ; stride=2
        lon = self.xcoord
        lat = self.ycoord

        ### get lon and lat in 2D version. get lat/lon intervals
        #numdim = len(np.array(lon).shape)
        #if numdim == 2:     [lon2d,lat2d] = [lon,lat]
        #elif numdim == 1:   [lon2d,lat2d] = np.meshgrid(lon,lat)
        #else:               errormess("lon and lat arrays must be 1D or 2D")
        #[wlon,wlat] = myplot.latinterv()
        ### define projection and background. define x and y given the projection
        #m = basemap.Basemap(projection='moll') marche pas
        #m = myplot.define_proj(proj,wlon,wlat,back=zeback,blat=None,blon=None)
        #x, y = m(lon2d, lat2d)
        ### TEMP
        x = lon ; y = lat
        ## define field. bound field.
        what_I_plot = np.transpose(field)
        zevmin, zevmax = myplot.calculate_bounds(what_I_plot)  ## vmin=min(what_I_plot_frame), vmax=max(what_I_plot_frame))
        what_I_plot = myplot.bounds(what_I_plot,zevmin,zevmax)
        ## define contour field levels. define color palette
        ticks = ndiv + 1
        zelevels = np.linspace(zevmin,zevmax,ticks)
        palette = get_cmap(name=colorb)
        ## contours topo
        zelevc = np.linspace(-8000.,20000.,20)
        yeah.contour( x, y, np.transpose(fieldc), zelevc, colors='black',linewidths = 0.4)
        # contour field
        c = yeah.contourf( x, y, what_I_plot, zelevels, cmap = palette, alpha = trans )
        clb = Figure.colorbar(fig,c,orientation='vertical',format="%.1e")
        clb.set_label(fieldlab)
        ax = fig.gca() ; ax.set_title(fieldlab) ; ax.set_ylabel("Latitude") ; ax.set_xlabel("Longitude")
        ax.set_xticks(np.arange(-180,181,45)) ; ax.set_xbound(lower=self.lons, upper=self.lone)
        ax.set_yticks(np.arange(-90,91,30)) ; ax.set_ybound(lower=self.lats, upper=self.late)
        if incwind:
          [x2d,y2d] = np.meshgrid(x,y)
          yeah.quiver(x2d,y2d,np.transpose(windx),np.transpose(windy))
      self.gettitle()
      fig.text(0.5, 0.95, self.title, ha='center')
      fig.text(0.5, 0.01, self.ack, ha='center')
      canvas = FigureCanvasAgg(fig)
      # The size * the dpi gives the final image size
      #   a4"x4" image * 80 dpi ==> 320x320 pixel image
      canvas.print_figure(figname, dpi=80)

    def htmlplot2d(self,tabtodo,fixedlt=False,figname="temp.png",title=""):
    ### complete 2D figure with possible multiplots
    ### added in 10/2012 for online MCD
    ### see http://www.dalkescientific.com/writings/diary/archive/2005/04/23/matplotlib_without_gui.html
      from matplotlib.figure import Figure
      from matplotlib.backends.backend_agg import FigureCanvasAgg
      from matplotlib.cm import get_cmap
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      fig = Figure(figsize=(8,8)) ; subv,subh = myplot.definesubplot( len(tabtodo) , fig )

      for i in range(len(tabtodo)):
        yeah = fig.add_subplot(subv,subh,i+1)
        choice = tabtodo[i]

        if self.lons is not None:    self.lonalt(fixedlt=fixedlt)
        elif self.lats is not None:  self.latalt(fixedlt=fixedlt)

        (field, fieldlab) = self.definefield(choice)

        colorb="jet" ; ndiv=20 ; title=""

        ## define field. bound field.
        what_I_plot = np.transpose(field)
        zevmin, zevmax = myplot.calculate_bounds(what_I_plot)  ## vmin=min(what_I_plot_frame), vmax=max(what_I_plot_frame))
        what_I_plot = myplot.bounds(what_I_plot,zevmin,zevmax)
        ## define contour field levels. define color palette
        ticks = ndiv + 1
        zelevels = np.linspace(zevmin,zevmax,ticks)
        palette = get_cmap(name=colorb)
        # contour field
        c = yeah.contourf( self.xcoord, self.ycoord, what_I_plot, zelevels, cmap = palette )
        clb = Figure.colorbar(fig,c,orientation='vertical',format="%.1e")
        clb.set_label(fieldlab)
        ax = fig.gca() ; ax.set_ylabel(self.ylabel) ; ax.set_xlabel(self.xlabel)
        if self.zkey == 4: ax.set_yscale('log') ; ax.set_ylim(ax.get_ylim()[::-1])

        #ax.set_xticks(np.arange(-180,181,45)) ; ax.set_xbound(lower=self.lons, upper=self.lone)
        #ax.set_yticks(np.arange(-90,91,30)) ; ax.set_ybound(lower=self.lats, upper=self.late)

      self.gettitle()
      fig.text(0.5, 0.95, self.title, ha='center')
      fig.text(0.5, 0.01, self.ack, ha='center')
      canvas = FigureCanvasAgg(fig)
      # The size * the dpi gives the final image size
      #   a4"x4" image * 80 dpi ==> 320x320 pixel image
      canvas.print_figure(figname, dpi=80)


    ### TODO: makeplot2d, plot2d, passer plot settings

