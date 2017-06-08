###################################################
### A Python Class for the Mars Climate Database ###
### ---------------------------------------------###
### Aymeric SPIGA 17-21/04/2012                  ###
### ---------------------------------------------###
### (see mcdtest.py for examples of use)         ###
####################################################

import numpy as np
import matplotlib.pyplot as mpl

def errormess(text,printvar=None):
    print text
    if printvar is not None: print printvar
    exit()
    return

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
        self.xlabel = None ; self.ylabel = None ; self.title = ""
        self.vertplot = False
        self.fmt = "%.1e" 
        self.colorm = "jet"
        self.fixedlt = False
        self.zonmean = False
        self.min2d = None
        self.max2d = None
        self.dpi = 80.
        self.islog = False
        self.proj = False
        self.trans = 0.0
        self.iscontour = False
        self.plat = 0.0
        self.plon = 0.0
        self.latpoint = None
        self.lonpoint = None

    def toversion5(self,version="5.1"):
        self.name      = "MCD v"+version
        self.dset      = '/home/marshttp/MCD_v'+version+'/data/'
        self.extvarkey = np.ones(100)

    def viking1(self): self.name = "Viking 1 site. MCD v4.3 output" ; self.lat = 22.48 ; self.lon = -49.97 ; self.xdate = 97.
    def viking2(self): self.name = "Viking 2 site. MCD v4.3 output" ; self.lat = 47.97 ; self.lon = -225.74 ; self.xdate = 117.6

    def getdustlabel(self):
        if self.dust == 1: 
            self.dustlabel = "MY24 minimum solar scenario"
            if "v5" in self.name: self.dustlabel = "climatology average solar scenario"
        elif self.dust == 2: 
            self.dustlabel = "MY24 average solar scenario"
            if "v5" in self.name: self.dustlabel = "climatology minimum solar scenario"
        elif self.dust == 3: 
            self.dustlabel = "MY24 maximum solar scenario"
            if "v5" in self.name: self.dustlabel = "climatology maximum solar scenario"
        elif self.dust == 4: self.dustlabel = "dust storm minimum solar scenario"
        elif self.dust == 5: self.dustlabel = "dust storm average solar scenario"
        elif self.dust == 6: self.dustlabel = "dust storm maximum solar scenario"
        elif self.dust == 7: self.dustlabel = "warm scenario (dusty, maximum solar)"
        elif self.dust == 8: self.dustlabel = "cold scenario (low dust, minimum solar)"
        elif self.dust > 20: self.dustlabel = "Martian Year "+str(self.dust)+" scenario"

    def gettitle(self,oneline=False):
        self.getdustlabel()
        self.title = self.name + " with " + self.dustlabel + "."
        if self.datekey == 1:    self.title = self.title + " Ls " + str(self.xdate) + "deg."
        elif self.datekey == 0:  self.title = self.title + " JD " + str(self.xdate) + "."
        if not oneline: self.title = self.title + "\n"
        if self.lats is None:  self.title = self.title + " Latitude " + str(self.lat) + "N"
        if self.zonmean and self.lats is not None and self.xzs is not None: 
            self.title = self.title + "Zonal mean over all longitudes."
        elif self.lons is None: 
            self.title = self.title + " Longitude " + str(self.lon) + "E"
        if self.xzs is None:   
            self.vertunits()
            self.title = self.title + " Altitude " + str(self.xz) + " " + self.vunits
        if self.datekey == 1:
          if self.locts is None:
            self.title = self.title + " Local time " + str(self.loct) + "h"
            if not self.fixedlt:  self.title = self.title + " (at longitude 0) "
            else: self.title = self.title + " (at all longitudes) "

    def getextvarlab(self,num):
        whichfield = { \
        91: "Pressure (Pa)", \
        92: "Density (kg/m3)", \
        93: "Temperature (K)", \
        94: "W-E wind component (m/s)", \
        95: "S-N wind component (m/s)", \
        96: "Horizontal wind speed (m/s)", \
	1: "Radial distance from planet center (m)",\
	2: "Altitude above areoid (Mars geoid) (m)",\
	3: "Altitude above local surface (m)",\
	4: "orographic height (m) (surf alt above areoid)",\
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
	16: "daily max mean surface temperature (K)",\
	17: "daily min mean surface temperature (K)",\
	18: "surf. temperature RMS day to day variations (K)",\
	19: "surface pressure (Pa)",\
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
        ### MCD version 5 new variables. AS 12/2012.
        if "v5" in self.name:
          whichfield[30] = whichfield[34]
          whichfield[34] = "surface H2O ice layer (kg/m2, 0.5: perennial)"
          whichfield[29] = "Surface roughness length z0 (m)"
          whichfield[37] = "DOD RMS day to day variations"
          whichfield[38] = "Dust mass mixing ratio (kg/kg)"
          whichfield[39] = "Dust effective radius (m)"
          whichfield[44] =  whichfield[43]
          whichfield[43] =  whichfield[42]
          whichfield[42] =  whichfield[41]
          whichfield[41] =  whichfield[40]
          whichfield[40] = "Dust deposition on flat surface (kg m-2 s-1)"
          whichfield[45] = "Water ice effective radius (m)"
          whichfield[46] = "Convective PBL height (m)"
          whichfield[47] = "Max. upward convective wind within the PBL (m/s)"
          whichfield[48] = "Max. downward convective wind within the PBL (m/s)"
          whichfield[49] = "Convective vertical wind variance at level z (m2/s2)"
          whichfield[50] = "Convective eddy vertical heat flux at level z (m/s/K)"
          whichfield[51] = "Surface wind stress (Kg/m/s2)"
          whichfield[52] = "Surface sensible heat flux (W/m2) (<0 when flux from surf to atm.)"
          whichfield[53] = "R: Molecular gas constant (J K-1 kg-1)"
          whichfield[54] = "Air viscosity estimation (N s m-2)"
          whichfield[55] = "not used (set to zero)"
          whichfield[56] = "not used (set to zero)"
          whichfield[57] = "[CO2] vol. mixing ratio (mol/mol)"
          whichfield[58] = "[N2] vol. mixing ratio (mol/mol)"
          whichfield[59] = "[Ar] vol. mixing ratio (mol/mol)"
          whichfield[60] = "[CO] vol. mixing ratio (mol/mol)"
          whichfield[61] = "[O] vol. mixing ratio (mol/mol)"
          whichfield[62] = "[O2] vol. mixing ratio (mol/mol)"
          whichfield[63] = "[O3] vol. mixing ratio (mol/mol)"
          whichfield[64] = "[H] vol. mixing ratio (mol/mol)"
          whichfield[65] = "[H2] vol. mixing ratio (mol/mol)"
          whichfield[66] = "electron number density (cm-3)"
          whichfield[67] = "CO2 column (kg/m2)"
          whichfield[68] = "N2 column (kg/m2)"
          whichfield[69] = "Ar column (kg/m2)"
          whichfield[70] = "CO column (kg/m2)"
          whichfield[71] = "O column (kg/m2)"
          whichfield[72] = "O2 column (kg/m2)"
          whichfield[73] = "O3 column (kg/m2)"
          whichfield[74] = "H column (kg/m2)"
          whichfield[75] = "H2 column (kg/m2)"
          whichfield[76] = "Total Electronic Content (TEC) (m-2)"
        if num not in whichfield: errormess("Incorrect subscript in extvar.")
        dastuff = whichfield[num]
        expf = "%.1e"
        if "(K)" in dastuff:      self.fmt="%.0f"
        elif "effective radius" in dastuff: self.fmt=expf
        elif "(Pa)" in dastuff:   self.fmt=expf
        elif "(W/m2)" in dastuff: self.fmt="%.0f"
        elif "(m/s)" in dastuff:  self.fmt="%.1f"
        elif "(mol/mol)" in dastuff: self.fmt=expf
        elif "(kg/m2)" in dastuff: self.fmt=expf
        elif "(m)" in dastuff:    self.fmt="%.0f"
        else:                     self.fmt=expf
        return dastuff

    def convertlab(self,num):        
        ## a conversion from text inquiries to extvar numbers. to be completed.
        if num == "p": num = 91
        elif num == "rho": num = 92
        elif num == "t": num = 93
        elif num == "u": num = 94
        elif num == "v": num = 95
        elif num == "wind": num = 96
        elif num == "tsurf": num = 15
        elif num == "topo": num = 4
        elif num == "h": num = 13
        elif num == "ps": num = 19
        elif num == "tau": num = 36
        elif num == "mtot": 
            if "v5" in self.name:  num = 41 
            else:                  num = 40
        elif num == "icetot": 
            if "v5" in self.name:  num = 43
            else:                  num = 42
        elif num == "h2ovap": 
            if "v5" in self.name:  num = 42
            else:                  num = 41
        elif num == "h2oice": 
            if "v5" in self.name:  num = 44
            else:                  num = 43
        elif num == "cp": num = 8
        elif num == "rho_ddv": num = 10
        elif num == "ps_ddv": num = 22
        elif num == "p_ddv": num = 21
        elif num == "t_ddv": num = 23
        elif num == "u_ddv": num = 24
        elif num == "v_ddv": num = 25
        elif num == "w": num = 26
        elif num == "tsurfmx": num = 16
        elif num == "tsurfmn": num = 17
        elif num == "lwdown": num = 31
        elif num == "swdown": num = 32
        elif num == "lwup": num = 33
        elif num == "swup":
            if "v5" in self.name:  num = 30
            else:                  num = 34
        elif num == "tau": num = 36
        elif num == "tau_ddv":
            if "v5" in self.name:  num = 37
            else:                  num = 38
        elif num == "qdust":
            if "v5" in self.name:  num = 38
            else:                  num = 37
        elif num == "co2":
            if "v5" in self.name:  num = 57
            else:                  num = 45
        elif num == "o3": 
            if "v5" in self.name:  num = 63
            else:                  num = 44
        elif num == "o": 
            if "v5" in self.name:  num = 61
            else:                  num = 46
        elif num == "co": 
            if "v5" in self.name:  num = 60
            else:                  num = 48
        elif num == "visc": 
            if "v5" in self.name:  num = 54
            else:                  num = 50
        elif num == "co2ice": num = 35
        elif num == "n2":
            if "v5" in self.name:  num = 58
            else:                  num = 47
        elif num == "n2col":
            if "v5" in self.name:  num = 68
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "rdust":
            if "v5" in self.name:  num = 39
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "sdust":
            if "v5" in self.name:  num = 40
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "pbl":
            if "v5" in self.name:  num = 46
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "updraft":
            if "v5" in self.name:  num = 47
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "downdraft":
            if "v5" in self.name:  num = 48
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "pblwvar":
            if "v5" in self.name:  num = 49
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "pblhvar":
            if "v5" in self.name:  num = 50
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "stress":
            if "v5" in self.name:  num = 51
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "ar":
            if "v5" in self.name:  num = 59
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "o2":
            if "v5" in self.name:  num = 62
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "co2col":
            if "v5" in self.name:  num = 67
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "arcol":
            if "v5" in self.name:  num = 69
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "cocol":
            if "v5" in self.name:  num = 70
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "o3col":
            if "v5" in self.name:  num = 73
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "hydro":
            if "v5" in self.name:  num = 64
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "hydro2":
            if "v5" in self.name:  num = 65
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "e":
            if "v5" in self.name:  num = 66
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "ecol":
            if "v5" in self.name:  num = 76
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif num == "groundice":
            if "v5" in self.name:  num = 34
            else:                  num = 11 # an undefined variable to avoid misleading output
        elif not isinstance(num, np.int): errormess("field reference not found.")
        return num

###################
### One request ###
###################

    def update(self):
    # retrieve fields from MCD (call_mcd). more info in fmcd.call_mcd.__doc__
        ## sanity first
        self.loct = abs(self.loct)%24
        if self.locts is not None and self.locte is not None: 
            self.locts = abs(self.locts)%24
            self.locte = abs(self.locte)%24
            if self.locts == self.locte: self.locte = self.locts + 24
        ## ensure that local time = 0 when using Earth dates
        if self.datekey == 0: self.loct = 0.
        ## now MCD request
        if "v5.1" in self.name: from fmcd51 import call_mcd
        elif "v5.2" in self.name: from fmcd52 import call_mcd
        else: from fmcd import call_mcd
        (self.pres, self.dens, self.temp, self.zonwind, self.merwind, \
         self.meanvar, self.extvar, self.seedout, self.ierr) \
         = \
         call_mcd(self.zkey,self.xz,self.lon,self.lat,self.hrkey, \
             self.datekey,self.xdate,self.loct,self.dset,self.dust, \
             self.perturkey,self.seedin,self.gwlength,self.extvarkey )
        ## we use the end of extvar (unused) to store meanvar. this is convenient for getextvar(lab)
        self.extvar[90] = self.pres ; self.extvar[91] = self.dens
        self.extvar[92] = self.temp ; self.extvar[93] = self.zonwind ; self.extvar[94] = self.merwind
        self.extvar[95] = np.sqrt(self.extvar[93]**2 + self.extvar[94]**2) # calculate wind modulus
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
        if "v5.1" in self.name: name = "v51_" + name
        elif "v5.2" in self.name: name = "v52_" + name
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
        dastr = str(self.extvar[num-1])
        if dastr == "nan":   print "!!!! There is a problem, probably a value is requested below the surface !!!!"
        else:                print self.getextvarlab(num) + " ..... " + dastr

    def printallextvar(self):
    # print all extra MCD variables    
        if "v5" in self.name:  limit=76
        else:                  limit=50
        for i in range(limit): self.printextvar(i+1)

    def htmlprinttabextvar(self,tabtodo):
        self.fixedlt = True ## local time is real local time
        self.gettitle()
        print "<hr>"
        print self.title
        print "<hr>"
        print "<ul>"
        for i in range(len(tabtodo)): print "<li>" ; self.printextvar(tabtodo[i]) ; print "</li>"
        print "</ul>"
        print "<hr>"
        print self.ack
        print "<hr>"
        #print "SETTINGS<br />"
        #self.printcoord()
        #self.printset()

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
      ## fix for possibly slightly negative tracers
      if "(mol/mol)" in fieldlab or "(kg/kg)" in fieldlab or "(kg/m2)" in fieldlab or "(W/m2)" in fieldlab:
         ind = np.where(field < 1.e-30)
         if ind != -1: field[ind] = 0.e0 #1.e-30  ## 0 does not work everywhere.
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
#      if self.zkey == 1:   self.xlabel = "radius from centre of planet (m)"
#      elif self.zkey == 2: self.xlabel = "height above areoid (m) (MOLA zero datum)"
#      elif self.zkey == 3: self.xlabel = "height above surface (m)"
#      elif self.zkey == 4: self.xlabel = "pressure level (Pa)"
#      elif self.zkey == 5: self.xlabel = "altitude above mean Mars Radius(=3396000m) (m)"
      if self.zkey == 1:   self.xlabel = "radius from centre of planet (m)"
      elif self.zkey == 2: self.xlabel = "altitude above MOLA$_0$ (m)"
      elif self.zkey == 3: self.xlabel = "height above surface (m)"
      elif self.zkey == 4: self.xlabel = "pressure (Pa)"
      elif self.zkey == 5: self.xlabel = "altitude above mean Mars radius (m)"

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
      if self.prestab is None:  errormess("arrays must be prepared first through self.prepare")
      self.prestab[i] = self.pres ; self.denstab[i] = self.dens ; self.temptab[i] = self.temp
      self.zonwindtab[i] = self.zonwind ; self.merwindtab[i] = self.merwind
      self.meanvartab[i,1:5] = self.meanvar[0:4]  ## note: var numbering according to MCD manual is kept
      self.extvartab[i,1:100] = self.extvar[0:99] ## note: var numbering according to MCD manual is kept

    def diurnal(self,nd=13):
    ### retrieve a local time slice
      self.fixedlt = True  ## local time is real local time
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
      if not self.fixedlt: umst = self.loct
      for i in range(nd): 
          self.lon = self.xcoord[i]
          if not self.fixedlt: self.loct = (umst + self.lon/15.) % 24
          self.update() ; self.put1d(i)
      self.lon = save

    def meridional(self,nd=19):
    ### retrieve a latitude slice
      self.fixedlt = True  ## local time is real local time
      save = self.lat
      self.xlabel = "North latitude (degrees)"
      self.prepare(ndx=nd) ; self.ininterv(-90.,90.,nd,start=self.lats,end=self.late)
      for i in range(nd): self.lat = self.xcoord[i] ; self.update() ; self.put1d(i)
      self.lat = save

    def profile(self,nd=20,tabperso=None):
    ### retrieve an altitude slice (profile)
      self.fixedlt = True  ## local time is real local time
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

    def getascii(self,tabtodo,filename="output.txt"):
    ### print out values in an ascii file
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      asciifile = open(filename, "w")
      for i in range(len(tabtodo)):  
          (field, fieldlab) = self.definefield(tabtodo[i])
          self.gettitle(oneline=True)
          asciifile.write("### " + self.title + "\n")
          asciifile.write("### " + self.ack + "\n")
          asciifile.write("### Column 1 is " + self.xlabel + "\n")
          dim = field.ndim
          if (dim == 1):
            asciifile.write("### Column 2 is " + fieldlab + "\n")
            for ix in range(len(self.xcoord)):
              asciifile.write("%15.5e%15.5e\n" % ( self.xcoord[ix], field[ix] ) )
          elif (dim == 2):
            asciifile.write("### Columns 2+ are " + fieldlab + "\n")
            asciifile.write("### Line 1 is " + self.ylabel + "\n")
            asciifile.write("---- ||")
            for iy in range(len(self.ycoord)):
              asciifile.write("%15.5e" % ( self.ycoord[iy] ) )
            asciifile.write("\n-----------------------------------\n") 
            for ix in range(len(self.xcoord)):
             zestr = "%+.03d ||" % (self.xcoord[ix])
             for iy in range(len(self.ycoord)):
               zestr = zestr + "%15.5e" % (field[ix,iy])
             asciifile.write(zestr+"\n")
      asciifile.close()
      return 

    def makeplot1d(self,choice):
    ### one 1D plot is created for the user-defined variable in choice. 
      (field, fieldlab) = self.definefield(choice)
      if not self.vertplot:  absc = self.xcoord ; ordo = field ; ordolab = fieldlab ; absclab = self.xlabel
      else:                  ordo = self.xcoord ; absc = field ; absclab = fieldlab ; ordolab = self.xlabel
      mpl.plot(absc,ordo,'-bo') ; mpl.ylabel(ordolab) ; mpl.xlabel(absclab) #; mpl.xticks(query.xcoord)
      # cases with log axis
      if self.zkey == 4: mpl.semilogy() ; ax = mpl.gca() ; ax.set_ylim(ax.get_ylim()[::-1])
      if not self.vertplot and self.islog: mpl.semilogy()
      if self.vertplot and self.islog: mpl.semilogx()
      mpl.figtext(0.5, 0.01, self.ack, ha='center')

    def plot1d(self,tabtodo):
    ### complete 1D figure with possible multiplots
      import mcdcomp as mcdcomp
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      fig = mpl.figure() ; subv,subh = mcdcomp.definesubplot( len(tabtodo) , fig ) 
      for i in range(len(tabtodo)): mpl.subplot(subv,subh,i+1).grid(True, linestyle=':', color='grey') ; self.makeplot1d(tabtodo[i])
      mpl.show()

    def htmlplot1d(self,tabtodo,figname="temp.png",title=""):
    ### complete 1D figure with possible multiplots
    ### added in 09/2012 for online MCD
    ### see http://www.dalkescientific.com/writings/diary/archive/2005/04/23/matplotlib_without_gui.html
      import mcdcomp as mcdcomp
      from matplotlib.figure import Figure
      from matplotlib.backends.backend_agg import FigureCanvasAgg
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.

      howmanyplots = len(tabtodo)
      if howmanyplots == 1: fig = Figure(figsize=(16,8))
      elif howmanyplots == 2: fig = Figure(figsize=(8,8))
      elif howmanyplots == 3: fig = Figure(figsize=(8,16))
      elif howmanyplots == 4: fig = Figure(figsize=(16,8))

      subv,subh = mcdcomp.definesubplot( len(tabtodo) , fig )
      for i in range(len(tabtodo)):
        yeah = fig.add_subplot(subv,subh,i+1) #.grid(True, linestyle=':', color='grey') 
        choice = tabtodo[i]
        (field, fieldlab) = self.definefield(choice)

        # in log plots we do not show the negative values
        if self.islog: field[np.where(field <= 0.e0)] = np.nan

        if not self.vertplot:  absc = self.xcoord ; ordo = field ; ordolab = fieldlab ; absclab = self.xlabel
        else:                  ordo = self.xcoord ; absc = field ; absclab = fieldlab ; ordolab = self.xlabel

        yeah.plot(absc,ordo,'-bo') #; mpl.xticks(query.xcoord)
        ax = fig.gca() ; ax.set_ylabel(ordolab) ; ax.set_xlabel(absclab)

        if self.xzs is not None and self.zkey == 4: ax.set_yscale('log') ; ax.set_ylim(ax.get_ylim()[::-1])
        if not self.vertplot and self.islog: ax.set_yscale('log')
        if self.vertplot and self.islog: ax.set_xscale('log')

        if self.lats is not None:      ax.set_xticks(np.arange(-90,91,15)) ; ax.set_xbound(lower=self.lats, upper=self.late)
        elif self.lons is not None:    ax.set_xticks(np.arange(-360,361,30)) ; ax.set_xbound(lower=self.lons, upper=self.lone)
        elif self.locts is not None:   ax.set_xticks(np.arange(0,26,2)) ; ax.set_xbound(lower=self.locts, upper=self.locte)

        ## does not work
        #ax.ticklabel_format(useOffset=False,axis='x')
        #ax.ticklabel_format(useOffset=False,axis='y')

        ax.grid(True, linestyle=':', color='grey')

      self.gettitle()
      fig.text(0.5, 0.95, self.title, ha='center')
      fig.text(0.5, 0.01, self.ack, ha='center')
      canvas = FigureCanvasAgg(fig)
      # The size * the dpi gives the final image size
      #   a4"x4" image * 80 dpi ==> 320x320 pixel image
      canvas.print_figure(figname, dpi=self.dpi)

###################
### 2D analysis ###
###################

    def latlon(self,ndx=37,ndy=19):
    ### retrieve a latitude/longitude slice
    ### default is: local time is not fixed. user-defined local time is at longitude 0.
      save1 = self.lon ; save2 = self.lat ; save3 = self.loct
      self.xlabel = "East longitude (degrees)" ; self.ylabel = "North latitude (degrees)"
      self.prepare(ndx=ndx,ndy=ndy)
      self.ininterv(-180.,180.,ndx,start=self.lons,end=self.lone)
      self.ininterv(-90.,  90.,ndy,start=self.lats,end=self.late,yaxis=True)
      if not self.fixedlt: umst = self.loct
      for i in range(ndx):
       for j in range(ndy):
         self.lon = self.xcoord[i] ; self.lat = self.ycoord[j]
         if not self.fixedlt: self.loct = (umst + self.lon/15.) % 24
         self.update() ; self.put2d(i,j)
      if not self.fixedlt: self.loct = umst
      self.lon = save1 ; self.lat = save2 ; self.loct = save3

    def secalt(self,ndx=37,ndy=20,typex="lat"):
    ### retrieve a coordinate/altitude slice
      save1 = self.lon ; save2 = self.xz ; save3 = self.loct ; save4 = self.lat
      self.prepare(ndx=ndx,ndy=ndy)
      self.vertlabel() ; self.ylabel = self.xlabel
      self.vertaxis(ndy,yaxis=True)
      if "lat" in typex:
          self.xlabel = "North latitude (degrees)"
          self.ininterv(-90.,90.,ndx,start=self.lats,end=self.late)
      elif typex == "lon":
          self.xlabel = "East longitude (degrees)"
          self.ininterv(-180.,180.,ndx,start=self.lons,end=self.lone)
      if not self.fixedlt: umst = self.loct
      for i in range(ndx):
       for j in range(ndy):
         if typex == "lat":   self.lat = self.xcoord[i]
         elif typex == "lon": self.lon = self.xcoord[i]
         self.xz = self.ycoord[j]
         if not self.fixedlt: self.loct = (umst + self.lon/15.) % 24
         self.update() ; self.put2d(i,j)
      if not self.fixedlt: self.loct = umst
      self.lon = save1 ; self.xz = save2 ; self.loct = save3 ; self.lat = save4

    def zonalmean(self,ndx=37,ndy=20,ndmean=32):
    ### retrieve a zonalmean lat/altitude slice
      self.fixedlt = False
      save1 = self.lon ; save2 = self.xz ; save3 = self.loct ; save4 = self.lat
      self.prepare(ndx=ndx,ndy=ndy)
      self.vertlabel() ; self.ylabel = self.xlabel
      self.vertaxis(ndy,yaxis=True)
      self.xlabel = "North latitude (degrees)"
      self.ininterv(-180.,180.,ndmean)
      coordmean = self.xcoord
      self.ininterv(-90.,90.,ndx,start=self.lats,end=self.late)
      umst = self.loct #fixedlt false for this case
      for i in range(ndx):
       self.lat = self.xcoord[i]
       for j in range(ndy):
        self.xz = self.ycoord[j]
        meanpres = 0. ; meandens = 0. ; meantemp = 0. ; meanzonwind = 0. ; meanmerwind = 0. ; meanmeanvar = np.zeros(5) ; meanextvar = np.zeros(100)        
        for m in range(ndmean):
           self.lon = coordmean[m]
           self.loct = (umst + self.lon/15.) % 24 #fixedlt false for this case
           self.update() 
           meanpres = meanpres + self.pres/float(ndmean) ; meandens = meandens + self.dens/float(ndmean) ; meantemp = meantemp + self.temp/float(ndmean)
           meanzonwind = meanzonwind + self.zonwind/float(ndmean) ; meanmerwind = meanmerwind + self.merwind/float(ndmean)
           meanmeanvar = meanmeanvar + self.meanvar/float(ndmean) ; meanextvar = meanextvar + self.extvar/float(ndmean)
        self.pres=meanpres ; self.dens=meandens ; self.temp=meantemp ; self.zonwind=meanzonwind ; self.merwind=meanmerwind
        self.meanvar=meanmeanvar ; self.extvar=meanextvar
        self.put2d(i,j)
      self.loct = umst #fixedlt false for this case
      self.lon = save1 ; self.xz = save2 ; self.loct = save3 ; self.lat = save4

    def hovmoller(self,ndtime=25,ndcoord=20,typex="lat"):
    ### retrieve a time/other coordinate slice
      save1 = self.lat ; save2 = self.xz ; save3 = self.loct ; save4 = self.lon
      if typex == "lat": 
          ndx = ndcoord ; self.xlabel = "North latitude (degrees)" 
          ndy = ndtime ; self.ylabel = "Local time (Martian hour)"
          self.prepare(ndx=ndx,ndy=ndy)
          self.ininterv(-90.,90.,ndx,start=self.lats,end=self.late)
          self.ininterv(0.,24.,ndy,start=self.locts,end=self.locte,yaxis=True)
      elif typex == "lon":
          ndx = ndcoord ; self.xlabel = "East longitude (degrees)"
          ndy = ndtime ; self.ylabel = "Local time (Martian hour)"
          self.prepare(ndx=ndx,ndy=ndy)
          self.ininterv(-180.,180.,ndx,start=self.lons,end=self.lone)
          self.ininterv(0.,24.,ndy,start=self.locts,end=self.locte,yaxis=True)
      elif typex == "alt":
          ndy = ndcoord ; self.vertlabel() ; self.ylabel = self.xlabel
          ndx = ndtime ; self.xlabel = "Local time (Martian hour)"
          self.prepare(ndx=ndx,ndy=ndy)
          self.vertaxis(ndy,yaxis=True)
          self.ininterv(0.,24.,ndx,start=self.locts,end=self.locte)
      for i in range(ndx):
       for j in range(ndy):
         if typex == "lat":   self.lat = self.xcoord[i] ; self.loct = self.ycoord[j]
         elif typex == "lon": self.lon = self.xcoord[i] ; self.loct = self.ycoord[j]
         elif typex == "alt": self.xz = self.ycoord[j] ; self.loct = self.xcoord[i]
         self.update() ; self.put2d(i,j)
      self.lat = save1 ; self.xz = save2 ; self.loct = save3 ; self.lon = save4

    def put2d(self,i,j):
    ## fill in subscript i,j in output arrays
    ## (arrays must have been correctly defined through prepare)
      if self.prestab is None:  errormess("arrays must be prepared first through self.prepare")
      self.prestab[i,j] = self.pres ; self.denstab[i,j] = self.dens ; self.temptab[i,j] = self.temp
      self.zonwindtab[i,j] = self.zonwind ; self.merwindtab[i,j] = self.merwind
      self.meanvartab[i,j,1:5] = self.meanvar[0:4]  ## note: var numbering according to MCD manual is kept
      self.extvartab[i,j,1:100] = self.extvar[0:99] ## note: var numbering according to MCD manual is kept

    def makemap2d(self,choice,incwind=False,proj="cyl"):
    ### one 2D map is created for the user-defined variable in choice.
      import mcdcomp as mcdcomp
      self.latlon() ## a map is implicitely a lat-lon plot. otherwise it is a plot (cf. makeplot2d)
      if choice == "wind" or incwind:
          (windx, fieldlabwx) = self.definefield("u")
          (windy, fieldlabwy) = self.definefield("v")
      if choice == "wind":
          field = np.sqrt(windx*windx + windy*windy)
          fieldlab = "Horizontal wind speed (m/s)"
      else:    
          (field, fieldlab) = self.definefield(choice)
      if incwind:   mcdcomp.maplatlon(self.xcoord,self.ycoord,field,title=fieldlab,proj=proj,vecx=windx,vecy=windy) #,stride=1)
      else:         mcdcomp.maplatlon(self.xcoord,self.ycoord,field,title=fieldlab,proj=proj)
      mpl.figtext(0.5, 0.0, self.ack, ha='center')

    def map2d(self,tabtodo,incwind=False,proj="cyl"):
    ### complete 2D figure with possible multiplots
      import mcdcomp as mcdcomp
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      fig = mpl.figure()
      subv,subh = mcdcomp.definesubplot( len(tabtodo) , fig ) 
      for i in range(len(tabtodo)): mpl.subplot(subv,subh,i+1) ; self.makemap2d(tabtodo[i],incwind=incwind,proj=proj)
      mpl.show() 

    def makeinterv(self):
      self.latinterv = 30.
      self.loninterv = 45.
      if self.lats is not None:
        if (abs(self.late-self.lats) < 90.): self.latinterv = 10.
        if (abs(self.late-self.lats) < 10.): self.latinterv = 1.
      if self.lons is not None:
        if (abs(self.lone-self.lons) < 135.): self.loninterv = 15.
        if (abs(self.lone-self.lons) < 15.): self.loninterv = 1.

    def htmlmap2d(self,tabtodo,incwind=False,figname="temp.png",back="zMOL"):
    ### complete 2D figure with possible multiplots
    ### added in 09/2012 for online MCD
    ### see http://www.dalkescientific.com/writings/diary/archive/2005/04/23/matplotlib_without_gui.html
      import mcdcomp as mcdcomp
      from matplotlib.figure import Figure
      from matplotlib.backends.backend_agg import FigureCanvasAgg
      from matplotlib.cm import get_cmap
      from matplotlib import rcParams 
      import matplotlib.pyplot as plt
      ####
      isproj = (self.proj is not None)
      if isproj:
        from mpl_toolkits.basemap import Basemap
      ####
      if (self.trans > 0.): isback = True
      else: isback = False
      if (self.trans > 0.99): self.iscontour = True
      ####

      try:
        from Scientific.IO import NetCDF
        filename = "/home/marshttp/surface.nc"
        zefile = NetCDF.NetCDFFile(filename, 'r') 
        fieldc = zefile.variables[back]
        yc = zefile.variables['latitude']
        xc = zefile.variables['longitude']
      except:
        print "Trouble with netCDF or surface.nc file. Continue without topo lines."

      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.

      howmanyplots = len(tabtodo)
      if howmanyplots == 1: fig = Figure(figsize=(16,8)) 
      elif howmanyplots == 2: fig = Figure(figsize=(8,8)) 
      elif howmanyplots == 3: fig = Figure(figsize=(8,16)) 
      elif howmanyplots == 4: fig = Figure(figsize=(16,8)) 

      subv,subh = mcdcomp.definesubplot( len(tabtodo) , fig )

      for i in range(len(tabtodo)):
        yeah = fig.add_subplot(subv,subh,i+1)
        choice = tabtodo[i]
        self.latlon(ndx=64,ndy=48) 
        ## a map is implicitely a lat-lon plot. otherwise it is a plot (cf. makeplot2d)
        (field, fieldlab) = self.definefield(choice)
        if incwind: (windx, fieldlabwx) = self.definefield("u") ; (windy, fieldlabwy) = self.definefield("v")

        ndiv = 40
        vecx=None ; vecy=None ; stride=2
        lon = self.xcoord ; lat = self.ycoord
     
        if isproj: 
          ##
          if self.plat is None: self.plat = 0.5*(self.lats+self.late)
          if self.plon is None: self.plon = 0.5*(self.lons+self.lone)
          ##
          if self.proj == "cyl": yeah = Basemap(projection=self.proj,\
                                                llcrnrlat=self.lats,urcrnrlat=self.late,\
                                                llcrnrlon=self.lons,urcrnrlon=self.lone,\
                                                ax=yeah,resolution=None)
          elif self.proj == "laea": yeah = Basemap(projection=self.proj,ax=yeah,resolution=None,\
                                                   lon_0=self.plon,lat_0=self.plat,\
                                                   llcrnrlat=self.lats,urcrnrlat=self.late,\
                                                   llcrnrlon=self.lons,urcrnrlon=self.lone)
          elif self.proj == "npstere": yeah = Basemap(projection=self.proj,boundinglat=self.plat,lon_0=self.plon,resolution=None,ax=yeah)
          elif self.proj == "spstere": yeah = Basemap(projection=self.proj,boundinglat=self.plat,lon_0=self.plon,resolution=None,ax=yeah)
          elif self.proj == "ortho": yeah = Basemap(projection=self.proj,lat_0=self.plat,lon_0=self.plon,resolution=None,ax=yeah)
          elif self.proj == "robin": yeah = Basemap(projection=self.proj,lon_0=0.,resolution=None,ax=yeah)
          ## NB: resolution=None is here to avoid loading coastlines which caused problems with some (plat,plon) couples
          ### background
          if isback:
            img="/home/marshttp/www-mars/mcd_python/MarsMap_2500x1250.jpg"
            yeah.warpimage(img,scale=0.75)
          mertab,partab = np.r_[-180.:180.:30.],np.r_[-90.:90.:15.]
          merlab,parlab = [0,0,0,1],[1,0,0,0]
          #format = '%.1f'

          yeah.drawmeridians(mertab,color="grey",labels=merlab)
          yeah.drawparallels(partab,color="grey",labels=parlab)
          [lon2d,lat2d] = np.meshgrid(lon,lat)
          x, y = yeah(lon2d, lat2d)
        else:
          x = lon ; y = lat

        ## define field. bound field.
        what_I_plot = np.transpose(field)
        zevmin, zevmax = mcdcomp.calculate_bounds(what_I_plot,vmin=self.min2d,vmax=self.max2d)  
        what_I_plot = mcdcomp.bounds(what_I_plot,zevmin,zevmax)
        ## define contour field levels. define color palette
        ticks = ndiv + 1
        zelevels = np.linspace(zevmin,zevmax,ticks)
        palette = get_cmap(name=self.colorm)

        ## topography contours
        try:
            rcParams['contour.negative_linestyle'] = 'solid' # negative contours solid instead of dashed
            zelevc = np.linspace(-9.,20.,11,0.)
            if isproj:
              [xc2,yc2] = np.meshgrid(xc,yc)
              xc,yc = yeah(xc2,yc2)
            else:
              yeah.contour( np.array(xc) + 360., yc, fieldc, zelevc, colors='black',linewidths = 0.4)
              yeah.contour( np.array(xc) - 360., yc, fieldc, zelevc, colors='black',linewidths = 0.4)
            ##
            yeah.contour( xc, yc, fieldc, zelevc, colors='black',linewidths = 0.4 )
        except:
            pass

        # contour field
        if self.iscontour: c = yeah.contour( x, y, what_I_plot, zelevels, cmap = palette )
        else: c = yeah.contourf( x, y, what_I_plot, zelevels, cmap = palette, alpha = 1.-self.trans )

        # colorbar
        if not isproj:               orientation='vertical'   ; frac = 0.15  ; pad = 0.04 ; lu = 0.5
        elif self.proj in ['moll']:  orientation="horizontal" ; frac = 0.08  ; pad = 0.03 ; lu = 1.0
        elif self.proj in ['robin']: orientation="horizontal" ; frac = 0.07  ; pad = 0.1  ; lu = 1.0
        elif self.proj in ['cyl']:   orientation="vertical"   ; frac = 0.023 ; pad = 0.03 ; lu = 0.5
        else:                        orientation='vertical'   ; frac = 0.05  ; pad = 0.03 ; lu = 0.5
        zelevpal = np.linspace(zevmin,zevmax,num=min([ticks/2+1,21]))
        clb = Figure.colorbar(fig,c,orientation=orientation,format=self.fmt,ticks=zelevpal,\
             fraction=frac,pad=pad,extend='both',spacing='proportional')
        clb.set_label(fieldlab)

        # wind vectors
        if incwind:
          if isproj: x2d,y2d = x,y
          else: [x2d,y2d] = np.meshgrid(lon,lat)
          wcolor = str(self.trans) # trans=0 black, trans=100 white
          yeah.quiver(x2d,y2d,np.transpose(windx),np.transpose(windy),color=wcolor)

        # add a point (TBD: text, use ax.annotate)
        if (self.lonpoint is not None) and (self.latpoint is not None):
          xpt,ypt = yeah(self.lonpoint,self.latpoint) # compute the native map projection coordinates
          yeah.plot(xpt,ypt,'go',markersize=8) # plot filled circle at the location

        # operation on axis
        ax = fig.gca()
        if not isproj:
          ax.set_ylabel("Latitude") ; ax.set_xlabel("Longitude")
          # make intervals 
          self.makeinterv()
          ax.set_xticks(np.arange(-360,361,self.loninterv)) ; ax.set_xbound(lower=self.lons, upper=self.lone)
          ax.set_yticks(np.arange(-90,91,self.latinterv)) ; ax.set_ybound(lower=self.lats, upper=self.late)

      ## titles and final production
      self.gettitle()
      #ax.set_title(self.title,x=0.5,y=1.05)
      #ax.set_xlabel('\n'+self.ack,x=0.5,y=0.05)

      fig.text(0.5, 0.95, self.title, ha='center')
      fig.text(0.5, 0.01, self.ack, ha='center')
      canvas = FigureCanvasAgg(fig)
      # The size * the dpi gives the final image size
      #   a4"x4" image * 80 dpi ==> 320x320 pixel image
      canvas.print_figure(figname, dpi=self.dpi) 
             #, bbox_inches='tight') removes title. and ax.set_title cannot set a global title for multiplots.



    def htmlplot2d(self,tabtodo,figname="temp.png"):
    ### complete 2D figure with possible multiplots
    ### added in 10/2012 for online MCD
    ### see http://www.dalkescientific.com/writings/diary/archive/2005/04/23/matplotlib_without_gui.html
      import mcdcomp as mcdcomp
      from matplotlib.figure import Figure
      from matplotlib.backends.backend_agg import FigureCanvasAgg
      from matplotlib.cm import get_cmap
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.

      howmanyplots = len(tabtodo)
      if howmanyplots == 1: fig = Figure(figsize=(16,8))
      elif howmanyplots == 2: fig = Figure(figsize=(8,8))
      elif howmanyplots == 3: fig = Figure(figsize=(8,16))
      elif howmanyplots == 4: fig = Figure(figsize=(16,8))

      subv,subh = mcdcomp.definesubplot( len(tabtodo) , fig )

      for i in range(len(tabtodo)):
        yeah = fig.add_subplot(subv,subh,i+1)
        choice = tabtodo[i]

        if self.lons is not None:    
           if self.locts is None:  self.secalt(ndx=64,ndy=35,typex="lon")
           else:                   self.hovmoller(ndcoord=64,typex="lon")
        elif self.lats is not None:  
           if self.locts is None:  
               if self.zonmean:   self.zonalmean()
               else:         self.secalt(ndx=48,ndy=35,typex="lat")
           else:                   self.hovmoller(ndcoord=48,typex="lat")
        else:
           self.hovmoller(ndcoord=35,typex="alt")

        (field, fieldlab) = self.definefield(choice)

        colorb=self.colorm ; ndiv=20 

        ## define field. bound field.
        what_I_plot = np.transpose(field)
        zevmin, zevmax = mcdcomp.calculate_bounds(what_I_plot,vmin=self.min2d,vmax=self.max2d)  
        what_I_plot = mcdcomp.bounds(what_I_plot,zevmin,zevmax)
        ## define contour field levels. define color palette
        ticks = ndiv + 1
        zelevels = np.linspace(zevmin,zevmax,ticks)
        palette = get_cmap(name=colorb)
        # contour field
        c = yeah.contourf( self.xcoord, self.ycoord, what_I_plot, zelevels, cmap = palette )
        clb = Figure.colorbar(fig,c,orientation='vertical',format=self.fmt,ticks=np.linspace(zevmin,zevmax,num=min([ticks/2+1,21])))
        clb.set_label(fieldlab)
        ax = fig.gca() ; ax.set_ylabel(self.ylabel) ; ax.set_xlabel(self.xlabel)

        self.makeinterv()
        if self.lons is not None:   ax.set_xticks(np.arange(-360,361,self.loninterv)) ; ax.set_xbound(lower=self.lons, upper=self.lone)
        elif self.lats is not None: ax.set_xticks(np.arange(-90,91,self.latinterv)) ; ax.set_xbound(lower=self.lats, upper=self.late)

        if self.locts is not None: 
            if self.xzs is not None: ax.set_xticks(np.arange(0,26,2)) ; ax.set_xbound(lower=self.locts, upper=self.locte)
            else:                    ax.set_yticks(np.arange(0,26,2)) ; ax.set_ybound(lower=self.locts, upper=self.locte)

        if self.zkey == 4 and self.xzs is not None: 
            ax.set_yscale('log') ; ax.set_ylim(ax.get_ylim()[::-1])
        else:
            #ax.set_yticks(np.arange(self.xzs,self.xze,10000.)) ; 
            ax.set_ybound(lower=self.xzs, upper=self.xze)

      self.gettitle()
      fig.text(0.5, 0.95, self.title, ha='center')
      fig.text(0.5, 0.01, self.ack, ha='center')
      canvas = FigureCanvasAgg(fig)
      # The size * the dpi gives the final image size
      #   a4"x4" image * 80 dpi ==> 320x320 pixel image
      canvas.print_figure(figname, dpi=self.dpi)

    ### TODO: makeplot2d, plot2d, passer plot settings

