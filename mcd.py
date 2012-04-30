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

class mcd:

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
        self.name      = "MCD v4.3 output"
        self.dset      = '/home/aymeric/Science/MCD_v4.3/data/'
        ## 1. spatio-temporal coordinates
        self.lat       = 0.
        self.lon       = 0.
        self.loct      = 0.
        self.xdate     = 0.  # see datekey
        self.xz        = 10. # see zkey
        ## 1bis. related settings
        self.zkey      = 3  # specify that xz is the altitude above surface (m)
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

    def viking1(self): self.name = "Viking 1 site. MCD v4.3 output" ; self.lat = 22.48 ; self.lon = -49.97 ; self.xdate = 97.
    def viking2(self): self.name = "Viking 2 site. MCD v4.3 output" ; self.lat = 47.97 ; self.lon = -225.74 ; self.xdate = 117.6

    def getextvarlab(self,num):
        whichfield = { \
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
        if num not in whichfield: errormess("Incorrect subscript in extvar.")
        return whichfield[num]

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

    def printset(self):
    # print main settings
        print "zkey",self.zkey,"xz",self.xz,"lon",self.lon,"lat",self.lat,"hrkey",self.hrkey, \
              "xdate",self.xdate,"loct",self.loct,"dust",self.dust

    def getnameset(self):
    # set a name referring to settings [convenient for databases]
        name = str(self.zkey)+str(self.xz)+str(self.lon)+str(self.lat)+str(self.hrkey)+str(self.datekey)+str(self.xdate)+str(self.loct)+str(self.dust)
        return name

    def printcoord(self):
    # print requested space-time coordinates
        print "----------------------------------------------------------------"
        print "LAT",self.lat,"LON",self.lon,"LOCT",self.loct,"XDATE",self.xdate
        print "----------------------------------------------------------------"

    def printmeanvar(self):
    # print mean MCD variables
        print "Pressure = %5.3f pascals. " % (self.pres)
        print "Density = %5.3f kilograms per cubic meter. " % (self.dens)
        print "Temperature = %3.0f kelvins (%4.0f degrees celsius)." % (self.temp,self.temp-273.15)
        print "Zonal wind = %5.3f meters per second." % (self.zonwind)
        print "Meridional wind = %5.3f meters per second." % (self.merwind)

    def printextvar(self,num):
    # print extra MCD variables
        print self.getextvarlab(num) + " ---> " + str(self.extvar[num-1])

    def printallextvar(self):
    # print all extra MCD variables    
        for i in range(50): self.printextvar(i+1)

    def printmcd(self):
    # 1. call MCD 2. print settings 3. print mean vars
        self.update()
        self.printcoord()
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
    ### --- choice can be a MCD number for extvar
      if isinstance(choice, np.int):    field = self.getextvar(choice); fieldlab = self.getextvarlab(choice)
      else:
       if choice == "t": 	field = self.temptab ; fieldlab="Temperature (K)"
       elif choice == "p":  	field = self.prestab ; fieldlab="Pressure (Pa)"
       elif choice == "rho":	field = self.denstab ; fieldlab="Density (kg/m3)"
       elif choice == "u":	field = self.zonwindtab ; fieldlab="W-E wind component (m/s)"
       elif choice == "v":      field = self.merwindtab ; fieldlab="S-N wind component (m/s)"
       elif choice == "tsurf": 	field = self.getextvar(15); fieldlab="Surface temperature (K)"
       elif choice == "topo":	field = self.getextvar(4) ; fieldlab="Topography (m)"
       elif choice == "h":      field = self.getextvar(13); fieldlab = "Scale height (m)"
       elif choice == "ps":     field = self.getextvar(19); fieldlab = "Surface pressure (Pa)"
       elif choice == "olr":    field = self.getextvar(33); fieldlab = "Outgoing longwave radiation (W/m2)"
       elif choice == "tau":	field = self.getextvar(36); fieldlab = "Dust optical depth"
       elif choice == "mtot":   field = self.getextvar(40); fieldlab = "Water vapor column (kg/m2)"
       elif choice == "icetot": field = self.getextvar(42); fieldlab = "Water ice column (kg/m2)"
       elif choice == "ps_ddv": field = self.getextvar(22); fieldlab = "Surface pressure RMS day to day variations (Pa)"
       else:                    errormess("field reference not found.")
      return field,fieldlab

###################
### 1D analysis ###
###################

    def put1d(self,i):
    ## fill in subscript i in output arrays
    ## (arrays must have been correctly defined through prepare)
      if self.prestab is None:  errormess("arrays must be prepared first through self.prepare")
      self.prestab[i] = self.pres ; self.denstab[i] = self.dens ; self.temptab[i] = self.temp
      self.zonwindtab[i] = self.zonwind ; self.merwindtab[i] = self.merwind
      self.meanvartab[i,1:5] = self.meanvar[0:4]  ## note: var numbering according to MCD manual is kept
      self.extvartab[i,1:100] = self.extvar[0:99] ## note: var numbering according to MCD manual is kept

    def diurnal(self,nd=13,start=0.,end=24.):
    ### retrieve a local time slice
      self.xlabel = "Local time (Martian hour)"
      self.prepare(ndx=nd) ; self.xcoord = np.linspace(start,end,nd)
      for i in range(nd): self.loct = self.xcoord[i] ; self.update() ; self.put1d(i)

    def zonal(self,nd=37,start=-180.,end=180.):
    ### retrieve a longitude slice
      self.xlabel = "East longitude (degrees)"
      self.prepare(ndx=nd) ; self.xcoord = np.linspace(start,end,nd)
      for i in range(nd): self.lon = self.xcoord[i] ; self.update() ; self.put1d(i)

    def meridional(self,nd=19,start=-90.,end=90.):
    ### retrieve a latitude slice
      self.xlabel = "North latitude (degrees)"
      self.prepare(ndx=nd) ; self.xcoord = np.linspace(start,end,nd)
      for i in range(nd): self.lat = self.xcoord[i] ; self.update() ; self.put1d(i)

    def profile(self,nd=20,start=0.,end=100000.):
    ### retrieve an altitude slice (profile)
      self.xlabel = "Altitude (m)"
      self.prepare(ndx=nd) ; self.xcoord = np.linspace(start,end,nd)
      for i in range(nd): self.xz = self.xcoord[i] ; self.update() ; self.put1d(i)

    def latlon(self,ndx=37,startx=-180.,endx=180.,ndy=19,starty=-90.,endy=90.):
    ### retrieve a latitude/longitude slice
      self.xlabel = "East longitude (degrees)" ; self.ylabel = "North latitude (degrees)"
      self.prepare(ndx=ndx,ndy=ndy)
      self.xcoord = np.linspace(startx,endx,ndx) ; self.ycoord = np.linspace(starty,endy,ndy)
      for i in range(ndx): 
       for j in range(ndy):
         self.lon = self.xcoord[i] ; self.lat = self.ycoord[j] ; self.update() ; self.put2d(i,j)

    def makeplot1d(self,choice,vertplot=0):
    ### one 1D plot is created for the user-defined variable in choice. 
      (field, fieldlab) = self.definefield(choice)
      if vertplot != 1:  absc = self.xcoord ; ordo = field ; ordolab = fieldlab ; absclab = self.xlabel
      else:              ordo = self.xcoord ; absc = field ; absclab = fieldlab ; ordolab = self.xlabel
      mpl.plot(absc,ordo,'-bo') ; mpl.ylabel(ordolab) ; mpl.xlabel(absclab) #; mpl.xticks(query.xcoord)

    def plot1d(self,tabtodo,vertplot=0):
    ### complete 1D figure with possible multiplots
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      fig = mpl.figure() ; subv,subh = myplot.definesubplot( len(tabtodo) , fig ) 
      for i in range(len(tabtodo)): mpl.subplot(subv,subh,i+1).grid(True, linestyle=':', color='grey') ; self.makeplot1d(tabtodo[i],vertplot)

###################
### 2D analysis ###
###################

    def put2d(self,i,j):
    ## fill in subscript i,j in output arrays
    ## (arrays must have been correctly defined through prepare)
      if self.prestab is None:  errormess("arrays must be prepared first through self.prepare")
      self.prestab[i,j] = self.pres ; self.denstab[i,j] = self.dens ; self.temptab[i,j] = self.temp
      self.zonwindtab[i,j] = self.zonwind ; self.merwindtab[i,j] = self.merwind
      self.meanvartab[i,j,1:5] = self.meanvar[0:4]  ## note: var numbering according to MCD manual is kept
      self.extvartab[i,j,1:100] = self.extvar[0:99] ## note: var numbering according to MCD manual is kept

    def makemap2d(self,choice):
    ### one 2D map is created for the user-defined variable in choice.
      (field, fieldlab) = self.definefield(choice)
      myplot.maplatlon(self.xcoord,self.ycoord,field,title=fieldlab,proj="moll")

    def map2d(self,tabtodo):
    ### complete 2D figure with possible multiplots
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      fig = mpl.figure() ; subv,subh = myplot.definesubplot( len(tabtodo) , fig ) 
      for i in range(len(tabtodo)): mpl.subplot(subv,subh,i+1) ; self.makemap2d(tabtodo[i])

    ### TODO: makeplot2d, plot2d, passer plot settings, vecteurs, plot loct pas fixe
