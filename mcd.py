####################################################
### A Python Class for the Mars Climate Database ###
### ---------------------------------------------###
### Aymeric SPIGA 17-21/04/2012                  ###
### ---------------------------------------------###
### (see mcdtest.py for examples of use)         ###
####################################################

###
from fmcd import mcd      # MCD compiled with f2py
from fmcd import dataloc  # location of MCD data file
from fmcd import dataver  # compiled version of MCD 
###
import numpy as np
import matplotlib.pyplot as mpl
###

## default number of points for each dimension
dfzon = 64 #37 # zonal dimension
dfmer = 48 #19 # meridional dimension
dfver = 35 #20 # vertical dimension
dflct = 25 #13 # local time dimension
dfsea = 25 # solar long dimension
## default names
lslab = "Areocentric longitude (degrees)"
latlab = "North latitude (degrees)"
lonlab = "East longitude (degrees)"
ltlab = "Local time (Martian hour)"
#NB: vertical labels are treated by .vertlabel()

def errormess(text,printvar=None):
    print text
    if printvar is not None: print printvar
    exit()
    return

class mcd_class():
 
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
        self.name      = "MCD_v"+dataver # this is coming from fmcd
        self.dset      = dataloc # this is coming from fmcd
        self.ack       = "Mars Climate Database (c) LMD/OU/IAA/ESA/CNES"
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
          
        self.dust      = 1 # climatological average scenario       
        self.hrkey     = 1  #set high resolution mode on (hrkey=0 to set high resolution off)
     
        ## 3. additional settings for advanced use
        
        self.extvarkey = np.ones(100) # now a table since MCD version 5
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
        self.typex = None ; self.typey = None
        self.xlabel = None ; self.ylabel = None ; self.title = ""
        self.vertplot = False
        self.fmt = "%.1e" 
        self.colorm = "jet"
        self.fixedlt = False
        self.averaging = None
        self.min2d = None
        self.max2d = None
        self.dpi = 80.
        self.islog = False
        self.proj = None
        self.trans = 0.0
        self.iscontour = False
        self.plat = 0.0
        self.plon = 0.0
        self.palt = None             
        self.latpoint = None
        self.lonpoint = None

    def viking1(self): self.name = "Viking 1 site. MCD v4.3 output" ; self.lat = 22.48 ; self.lon = -49.97 ; self.xdate = 97.
    def viking2(self): self.name = "Viking 2 site. MCD v4.3 output" ; self.lat = 47.97 ; self.lon = -225.74 ; self.xdate = 117.6

    def getdustlabel(self):
        if self.dust == 1: 
            self.dustlabel = "climatology average solar scenario"                        
        elif self.dust == 2: 
            self.dustlabel = "climatology minimum solar scenario"
        elif self.dust == 3: 
            self.dustlabel = "climatology maximum solar scenario"
        elif self.dust == 4: self.dustlabel = "dust storm minimum solar scenario"
        elif self.dust == 5: self.dustlabel = "dust storm average solar scenario"
        elif self.dust == 6: self.dustlabel = "dust storm maximum solar scenario"
        elif self.dust == 7: self.dustlabel = "warm scenario (dusty, maximum solar)"
        elif self.dust == 8: self.dustlabel = "cold scenario (low dust, minimum solar)"
        elif self.dust > 20: self.dustlabel = "Martian Year "+str(self.dust)+" scenario"

    def gettitle(self,oneline=False):
        self.getdustlabel()
        self.title = self.name + " with " + self.dustlabel + "."
        if self.datekey == 1:    
          if self.xdates is None:
           self.title = self.title + " Ls " + str(self.xdate) + "deg."
        elif self.datekey == 0:  
          self.title = self.title + " JD " + str(self.xdate) + "."
        if not oneline: self.title = self.title + "\n"
        if self.lats is None:  
            self.title = self.title + " Latitude " + str(self.lat) + "N."
        if self.averaging == "lon": 
            self.title = self.title + " Zonal mean over all longitudes."
        elif self.lons is None: 
            self.title = self.title + " Longitude " + str(self.lon) + "E."
        if self.xzs is None:   
            self.vertunits()
            self.title = self.title + " Altitude " + str(self.xz) + " " + self.vunits + "."
        if self.datekey == 1:
          if self.averaging == "loct":
            self.title = self.title + " Diurnal mean over all local times."
          else:
            if self.locts is None and self.averaging != "lon":
              self.title = self.title + " Local time " + str(self.loct) + "h"
              if self.lons is not None: # if longitude is a free dimension
                if not self.fixedlt:  self.title = self.title + " (at longitude 0) "
                else: self.title = self.title + " (fixed at all longitudes) "
        if self.proj == "nsper" : 
            self.title = self.title + "[view from lon "+str(self.plon)+"$\degree$E. lat "+str(self.plat)+"$\degree$N. alt "+str(self.palt)+ " km]"

    def getextvarlab(self,num):
        # MCD version 6.1 variables
        whichfield = { \
        91: "Pressure (Pa)", \
        92: "Density (kg/m3)", \
        93: "Temperature (K)", \
        94: "W-E wind component (m/s)", \
        95: "S-N wind component (m/s)", \
        96: "Horizontal wind speed (m/s)", \
	1: "Radial distance from planet center (m)", \
	2: "Altitude above areoid (Mars geoid) (m)", \
	3: "Altitude above local surface (m)", \
	4: "orographic height (m) (surface altitude above areoid)", \
	5: "GCM orography (m)", \
	6: "Local slope inclination (deg) (HR mode only)", \
	7: "Local slope orientation (deg) (0 deg Northward) (HR mode only)", \
	8: "Sun-Mars distance (in Astronomical Unit AU)", \
	9: "Ls, solar longitude of Mars (deg)", \
	10: "LST:Local true solar time (hrs)", \
	11: "LMT:Local mean time (hrs) at sought longitude", \
	12: "Universal solar time (LST at lon=0) (hrs)", \
	13: "Solar zenith angle (deg)", \
	14: "Surface temperature (K)", \
	15: "Surface pressure (Pa)", \
	16: "GCM surface pressure (Pa)", \
	17: "Potential temperature (K) (reference pressure=610Pa)", \
	18: "Vertical wind component (m/s) (Up-Down)", \
	19: "Zonal slope wind component (m/s) (HR mode only)", \
	20: "Meridional slope wind component (m/s) (HR mode only)", \
	21: "Surface pressure RMS day to day variations (Pa)", \
	22: "Surface temperature RMS day to day variations (K)", \
	23: "Atmospheric pressure RMS day to day variations (Pa)", \
	24: "Density RMS day to day variations (kg/m^3)", \
	25: "Temperature RMS day to day variations (K)", \
	26: "Zonal wind RMS day to day variations (m/s)", \
	27: "Meridional wind RMS day to day variations (m/s)", \
	28: "Vertical wind RMS day to day variations (m/s)", \
	29: "Incident solar flux at top of the atmosphere (W/m2)", \
	30: "solar flux reflected to space (W/m2)", \
	31: "Incident solar flux on horizontal surface (W/m2)", \
	32: "Incident solar flux on local slope (W/m2) (HR mode only)", \
	33: "Reflected solar flux on horizontal surface (W/m2)", \
	34: "thermal IR flux to space (W/m2)", \
	35: "thermal IR flux on surface (W/m2)", \
	36: "GCM surface roughness length z0 (m)", \
	37: "GCM surface thermal inertia", \
	38: "GCM surface bare ground albedo", \
	39: "Monthly mean dust column visible optical depth above surface", \
	40: "Daily mean dust column visible optical depth above surface", \
	41: "Dust mass mixing ratio (kg/kg)", \
	42: "Dust effective radius (m)", \
	43: "Daily mean dust deposition rate on horizontal surface (kg m-2 s-1)", \
	44: "Monthly mean surface CO2 ice layer (kg/m2)", \
	45: "Monthly mean surface H2O layer (kg/m2) (non perennial frost)", \
	46: "GCM perennial surface water ice (0 or 1)", \
	47: "Water vapor column (kg/m2)", \
	48: "Water vapor vol. mixing ratio (mol/mol)", \
	49: "Water ice column (kg/m2)", \
	50: "Water ice mixing ratio (mol/mol)", \
	51: "Water ice effective radius (m)", \
	52: "Convective Planetary Boundary Layer (PBL) height (m)", \
	53: "Max. upward convective wind within the PBL (m/s)", \
	54: "Max. downward convective wind within the PBL (m/s)", \
	55: "Convective vertical wind variance at level z (m2/s2)", \
	56: "Convective eddy vertical heat flux at level z (m/s/K)", \
	57: "Surface wind stress (kg/m/s2)", \
	58: "Surface sensible heat flux (W/m2) (<0 when flux from surf to atm.)", \
	59: "Air heat capacity Cp (J kg-1 K-1)", \
	60: "gamma=Cp/Cv Ratio of specific heats", \
	61: "R:Molecular gas constant (J K-1 kg-1)", \
	62: "Air viscosity estimation (N s m-2)", \
	63: "Scale height H(p) (m)", \
	64: "[CO2] volume mixing ratio (mol/mol)", \
	65: "[N2] volume mixing ratio  (mol/mol)", \
	66: "[Ar] volume mixing ratio  (mol/mol)", \
	67: "[CO] volume mixing ratio  (mol/mol)", \
	68: "[O] volume mixing ratio   (mol/mol)", \
	69: "[O2] volume mixing ratio  (mol/mol)", \
	70: "[O3] volume mixing ratio  (mol/mol)", \
	71: "[H] volume mixing ratio   (mol/mol)", \
	72: "[H2] volume mixing ratio  (mol/mol)", \
	73: "[He] volume mixing ratio  (mol/mol)", \
	74: "CO2 column (kg/m2)", \
	75: "N2 column  (kg/m2)", \
	76: "Ar column  (kg/m2)", \
	77: "CO column  (kg/m2)", \
	78: "O column   (kg/m2)", \
	79: "O2 column  (kg/m2)", \
	80: "O3 column  (kg/m2)", \
	81: "H column   (kg/m2)", \
	82: "H2 column  (kg/m2)", \
	83: "He column  (kg/m2)", \
	84: "Electron number density (particules/cm3)", \
	85: "Total electonic content (TEC) (particules/m2)"
        }

        if num not in whichfield: errormess("Incorrect subscript in extvar.")
        dastuff = whichfield[num]
        expf = "%.1e"
        if "variations (K)" in dastuff: self.fmt="%.1f"
        elif "(K)" in dastuff:      self.fmt="%.0f"
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


        if   num == "p":      num = 91
        elif num == "rho":  num = 92
        elif num == "t":    num = 93
        elif num == "u":    num = 94
        elif num == "v":    num = 95
        elif num == "wind": num = 96   
        elif num == "zradius": num=1
        elif num == "zareoid" : num =2         
        elif num == "zsurface" : num = 3         
        elif num == "oroheight" : num = 4         
        elif num == "oro_gcm" : num = 5          
        elif num == "theta_s" : num = 6         
        elif num == "psi_s" : num = 7          
        elif num == "marsau" : num = 8         
        elif num == "ls" : num = 9         
        elif num == "loctime" : num = 10         
        elif num == "lmeantime" : num = 11          
        elif num == "utime" : num = 12          
        elif num == "solzenang" : num = 13         
        elif num == "tsurf" : num = 14         
        elif num == "ps" : num = 15         
        elif num == "ps_gcm" : num = 16        
        elif num == "potential_temp" : num = 17         
        elif num == "w_l" : num = 18         
        elif num == "zonal_slope_wind" : num = 19          
        elif num == "merid_slope_wind" : num = 20         
        elif num == "rmsps" : num = 21         
        elif num == "rmstsurf" : num = 22          
        elif num == "altrmsp" : num = 23          
        elif num == "rmsrho" : num = 24         
        elif num == "rmst" : num = 25          
        elif num == "rmsu" : num = 26         
        elif num == "rmsv" : num = 27          
        elif num == "rmsw" : num = 28         
        elif num == "fluxtop_dn_sw" : num = 29          
        elif num == "fluxtop_up_sw" : num = 30          
        elif num == "fluxsurf_dn_sw" : num = 31        
        elif num == "fluxsurf_dn_sw_hr" : num = 32          
        elif num == "fluxsurf_up_sw" : num = 33        
        elif num == "fluxtop_lw" : num = 34         
        elif num == "fluxsurf_lw" : num = 35        
        elif num == "z_0" : num = 36          
        elif num == "thermal_inertia" : num = 37        
        elif num == "ground_albedo" : num = 38         
        elif num == "dod" : num = 39         
        elif num == "tauref" : num = 40          
        elif num == "dust_mmr" : num = 41         
        elif num == "dust_reff" : num = 42        
        elif num == "dust_dep" : num = 43         
        elif num == "co2ice" : num = 44         
        elif num == "surf_h2o_ice" : num = 45          
        elif num == "water_cap" : num = 46          
        elif num == "col_h2ovapor" : num = 47          
        elif num == "vmr_h2o" : num = 48         
        elif num == "col_h2oice" : num = 49          
        elif num == "vmr_h2oice" : num = 50          
        elif num == "h2oice_reff" : num = 51          
        elif num == "zmax" : num = 52          
        elif num == "wstar_up" : num = 53         
        elif num == "wstar_dn" : num = 54       
        elif num == "vvv" : num = 55        
        elif num == "vhf" : num = 56        
        elif num == "surfstress" : num = 57         
        elif num == "sensib_flux" : num = 58        
        elif num == "Cp" : num = 59         
        elif num == "gamma" : num = 60         
        elif num == "Rgas" : num = 61         
        elif num == "viscosity" : num = 62         
        elif num == "pscaleheight" : num = 63          
        elif num == "vmr_co2" : num = 64         
        elif num == "vmr_n2" : num = 65         
        elif num == "vmr_ar" : num = 66          
        elif num == "vmr_co" : num = 67         
        elif num == "vmr_o" : num = 68         
        elif num == "vmr_o2" : num = 69         
        elif num == "vmr_o3" : num = 70         
        elif num == "vmr_h" : num = 71         
        elif num == "vmr_h2" : num = 72         
        elif num == "vmr_he" : num = 73         
        elif num == "col_co2" : num = 74         
        elif num == "col_n2" : num = 75         
        elif num == "col_ar" : num = 76         
        elif num == "col_co" : num = 77         
        elif num == "col_o" : num = 78          
        elif num == "col_o2" : num = 79          
        elif num == "col_o3" : num = 80         
        elif num == "col_h" : num = 81         
        elif num == "col_h2" : num = 82         
        elif num == "col_he" : num = 83         
        elif num == "vmr_elec" : num = 84          
        elif num == "col_elec" : num = 85         
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
        ### now MCD request
        (self.pres, self.dens, self.temp, self.zonwind, self.merwind, \
         self.meanvar, self.extvar, self.seedout, self.ierr) \
         = \
         mcd.call_mcd(self.zkey,self.xz,self.lon,self.lat,self.hrkey, \
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
        strdate = str(self.datekey)+str(self.xdate)+str(self.xdates)+str(self.xdatee)
        name = str(self.zkey)+strxz+strlon+strlat+str(self.hrkey)+strdate+strloct+str(self.dust)
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
        limit=85 #mcd6 number of extravar
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
    ### --- prepare I/O arrays for 1d slices ---
    ### NB: use ininterv to define xcoord and ycoord, not done here
      if ndx is None:  print "No dimension in prepare. Exit. Set at least ndx." ; exit()
      #else:            self.xcoord = np.ones(ndx)
      if ndy is None:  dashape = (ndx)     ; dashapemean = (ndx,6)     ; dashapeext = (ndx,101)     #; self.ycoord = None
      else:            dashape = (ndx,ndy) ; dashapemean = (ndx,ndy,6) ; dashapeext = (ndx,ndy,101) #; self.ycoord = np.ones(ndy)
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

    def diurnal(self,nd=dflct):
    ### retrieve a local time slice
      self.fixedlt = True  ## local time is real local time
      save = self.loct
      self.xlabel = ltlab
      self.prepare(ndx=nd) ; self.ininterv(0.,24.,nd,start=self.locts,end=self.locte) 
      for i in range(nd): self.loct = self.xcoord[i] ; self.update() ; self.put1d(i)
      self.loct = save

    def zonal(self,nd=dfzon):
    ### retrieve a longitude slice
      save = self.lon
      save2 = self.loct
      self.xlabel = lonlab
      self.prepare(ndx=nd) ; self.ininterv(-180.,180.,nd,start=self.lons,end=self.lone)
      if not self.fixedlt: umst = self.loct
      for i in range(nd): 
          self.lon = self.xcoord[i]
          if not self.fixedlt: self.loct = (umst + self.lon/15.) % 24
          self.update() ; self.put1d(i)
      self.lon = save
      self.loct = save2

    def meridional(self,nd=dfmer):
    ### retrieve a latitude slice
      self.fixedlt = True  ## local time is real local time
      save = self.lat
      self.xlabel = latlab
      self.prepare(ndx=nd) ; self.ininterv(-90.,90.,nd,start=self.lats,end=self.late)
      for i in range(nd): self.lat = self.xcoord[i] ; self.update() ; self.put1d(i)
      self.lat = save

    def profile(self,nd=dfver,tabperso=None):
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

    def seasonal(self,nd=dfsea):
    ### retrieve a seasonal slice
      save = self.xdate
      self.xlabel = lslab
      self.prepare(ndx=nd) ; self.ininterv(0.,360.,nd,start=self.xdates,end=self.xdatee)
      for i in range(nd): self.xdate = self.xcoord[i] ; self.update() ; self.put1d(i)
      self.xdate = save

    def getascii(self,tabtodo,filename="output.txt",log=None):
    ### print out values in an ascii file
      import datetime
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      asciifile = open(filename, "w")
      if log is not None:
          logfile = open(log, "a")
      for i in range(len(tabtodo)):  
          txt = "##########################################################################################\n"
          ### print out header
          (field, fieldlab) = self.definefield(tabtodo[i])
          self.gettitle(oneline=True)
          txt = txt + "### " + self.title + "\n"
          txt = txt.replace("scenario.","scenario.\n###")
          txt = txt + "### --------------------------------------------------------------------------------------\n"
          txt = txt + "### Column 1 is " + self.xlabel + "\n"
          dim = field.ndim
          if (dim == 1):
            txt = txt + "### Column 2 is " + fieldlab + "\n"
            data = txt
          elif (dim == 2):
            txt = txt + "### Columns 2+ are " + fieldlab + "\n"
            txt = txt + "### Line 1 is " + self.ylabel + "\n"
          txt = txt + "### --------------------------------------------------------------------------------------\n"
          txt = txt + "### Retrieved on: " + datetime.datetime.today().isoformat() + "\n"
          txt = txt + "### " + self.ack + "\n"
          txt = txt + "##########################################################################################\n"
          ### print out data
          data = txt
          if (dim == 1):
            for ix in range(len(self.xcoord)):
              data = data + "%15.5e%15.5e\n" % ( self.xcoord[ix], field[ix] )
          elif (dim == 2):
            data = data + "---- ||"
            for iy in range(len(self.ycoord)):
              data = data + "%15.5e" % ( self.ycoord[iy] )
            data = data + "\n-----------------------------------\n"
            for ix in range(len(self.xcoord)):
             zestr = "%+.03d ||" % (self.xcoord[ix])
             for iy in range(len(self.ycoord)):
               zestr = zestr + "%15.5e" % (field[ix,iy])
             data = data + zestr+"\n"
          asciifile.write(data)
          if log is not None:
             logfile.write(txt)
      asciifile.close()
      if log is not None:
         logfile.close()
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
        elif self.xdates is not None:  ax.set_xticks(np.arange(0,361,30)) ; ax.set_xbound(lower=self.xdates, upper=self.xdatee)

        ## does not work
        #ax.ticklabel_format(useOffset=False,axis='x')
        #ax.ticklabel_format(useOffset=False,axis='y')

        ax.grid(True, linestyle=':', color='grey')

      self.gettitle()
      fig.text(0.5, 0.95, self.title, ha='center', transform=fig.gca().transAxes, fontweight='bold')
      fig.text(0.5, 0.01, self.ack, ha='center')
      canvas = FigureCanvasAgg(fig)
      # The size * the dpi gives the final image size
      #   a4"x4" image * 80 dpi ==> 320x320 pixel image
      canvas.print_figure(figname, dpi=self.dpi)

###################
### 2D analysis ###
###################

    def definetype(self,typex=None,typey=None):
      ### get the type of 2D plot
      # -- we keep a "forcing" option for retrocompatibility (e.g. for maps)
      # -- this is meant to ultimately disappear
      if typex is not None and typey is not None:
        self.typex = typex
        self.typey = typey
        return
      # -- retrieve kind of time axis
      if self.locts is not None: 
        zetypet = "loct"
        # case of a ls/lt map (all-time plot)
        if self.xdates is not None:
          self.typex = "ls" ; self.typey = zetypet
          return
      elif self.xdates is not None: 
        zetypet = "ls"
      else: 
        zetypet = None
      # -- retrieve kind of spatial axis
      if self.lons is not None: 
        zetypes = "lon"
        # case of a lat/lon map
        if self.lats is not None:
          self.typex = zetypes ; self.typey = "lat"
          return      
      elif self.lats is not None: 
        zetypes = "lat"
      # -- explore all possibilities
      if zetypet is None:
        # this is horizontal/vertical section
        self.typex = zetypes ; self.typey = "alt" # secalt
      else:
        # there is a time axis
        # -- first we consider that spatial could be altitude
        # -- others (lat & lon) were done above
        if self.xzs is not None: 
          zetypes = "alt"
        # -- second we organize coordinates as usual
        if zetypet == "loct":
          if zetypes == "alt":
            self.typex = zetypet ; self.typey = zetypes # hovmoller with alt
          else:
            self.typex = zetypes ; self.typey = zetypet # hovmoller
        elif zetypet == "ls":
          self.typex = zetypet ; self.typey = zetypes # seasonal
      return

    def query2d(self,typex=None,typey=None,ndmean=32):
    ### retrieve a 2D slice
      # save query
      save1 = self.lon ; save2 = self.xz ; save3 = self.loct ; save4 = self.lat ; save5 = self.xdate
      # initialize
      if not self.fixedlt: 
        umst = self.loct # local time is not fixed. user-defined local time is at longitude 0.
      # define the type of 2D plot
      # -- hard setting of typex and typey is meant to disappear
      self.definetype(typex=typex,typey=typey)
      # settings for averaging
      if self.averaging is not None:
        coordmean = self.meandim(ndmean=ndmean) 
        if self.averaging == "lon": 
          self.fixedlt = False
          umst = self.loct # shouldn't it be zero? does not matter anyways...
      # get coordinates
      self.fillcoord()
      # MAIN QUERY LOOP
      for i in range(self.xcoord.size):
       for j in range(self.ycoord.size):
         # fill in correct values for query
         self.filldim(self.xcoord[i],self.ycoord[j])  
         # emulate "natural" global local time if not free dim and not fixedlt
         if self.locts is None:
           if not self.fixedlt: 
             self.loct = (umst + self.lon/15.) % 24
         # get field (simple or average)
         if self.averaging is None:
           self.update()
         else:
           self.meanperform(coordmean,umst=umst) 
         # fill in 2D array
         self.put2d(i,j)
      # reinstall init state
      if not self.fixedlt: self.loct = umst
      self.lon = save1 ; self.xz = save2 ; self.loct = save3 ; self.lat = save4 ; self.xdate = save5

################################################
### retro-compatibility functions
### -- only latlon is used now in htmlmap2d query
### -- htmlplot2d uses now query2d() in which plot type is set
### -- but what is below could be useful in interactive mode
    def latlon(self,typex="lon",typey="lat"):
    ### retrieve a latitude/longitude slice
      self.query2d(typex=typex,typey=typey)
    def secalt(self,typex="lat"):
    ### retrieve a coordinate/altitude slice
      self.query2d(typex=typex,typey="alt")
    def hovmoller(self,typex="lat",typey="loct"):
    ### retrieve a time/other coordinate slice
      if typey == "ls":
         typey = typex ; typex = "ls" #usually in seasonal plots, Ls is in x-axis
      self.query2d(typex=typex,typey=typey)
    def zonalmean(self,ndmean=32,typey="alt",typex="lat"):
    ### retrieve a zonalmean lat/altitude or ls/lat slice
      self.averaging="lon"
      self.query2d(typex=typex,typey=typey)
################################################

    def fillcoord(self):
    ### define coordinates. using global default values for number of points.
    ### prepare output arrays once coordinates are defined.
      # fill in coord properties // x-axis
      if self.typex == "lat":
        ndx = dfmer
        self.xlabel = latlab
        self.ininterv(-90.,90.,ndx,start=self.lats,end=self.late) ## do we want to change ndx or just use dfzon?
      elif self.typex == "lon":
        ndx = dfzon
        self.xlabel = lonlab
        self.ininterv(-180.,180.,ndx,start=self.lons,end=self.lone)
      elif self.typex == "alt":
        ndx = dfver
        self.vertlabel()
        self.vertaxis(dfver)
      elif self.typex == "ls":
        ndx = dfsea
        self.xlabel = lslab
        self.ininterv(0.,360.,ndx,start=self.xdates,end=self.xdatee)
      elif self.typex == "loct":
        ndx = dflct
        self.xlabel = ltlab
        self.ininterv(0.,24.,ndx,start=self.locts,end=self.locte)
      # fill in coord properties // y-axis
      if self.typey == "lat":
        ndy = dfmer
        self.ylabel = latlab
        self.ininterv(-90.,90.,ndy,start=self.lats,end=self.late,yaxis=True) 
      elif self.typey == "lon":
        ndy = dfzon
        self.ylabel = lonlab
        self.ininterv(-180.,180.,ndy,start=self.lons,end=self.lone,yaxis=True)
      elif self.typey == "alt":
        ndy = dfver
        sav = self.xlabel # save because used below as intermediate (to be improved)
        self.vertlabel() ; self.ylabel = self.xlabel
        self.xlabel = sav
        self.vertaxis(ndy,yaxis=True)
      elif self.typey == "ls":
        ndy = dfsea
        self.ylabel = lslab
        self.ininterv(0.,360.,ndy,start=self.xdates,end=self.xdatee,yaxis=True)
      elif self.typey == "loct":
        ndy = dflct
        self.ylabel = ltlab
        self.ininterv(0.,24.,ndy,start=self.locts,end=self.locte,yaxis=True)
      # prepare arrays with correct dimensions
      self.prepare(ndx=ndx,ndy=ndy)

    def filldim(self,xx,yy):
      # fill in individual values in axis
      tab = [[self.typex,xx],[self.typey,yy]]
      for ttt in tab:
        if "lat"  in ttt[0]: self.lat   = ttt[1]
        if "lon"  in ttt[0]: self.lon   = ttt[1]
        if "alt"  in ttt[0]: self.xz    = ttt[1]
        if "loct" in ttt[0]: self.loct  = ttt[1]
        if "ls"   in ttt[0]: self.xdate = ttt[1]

    def meandim(self,ndmean=32):
    ### define averaging dimension
      sav = self.xcoord #using xcoord as an intermediate
      if self.averaging == "lon":
        self.ininterv(-180.,180.,ndmean)
      elif self.averaging == "loct":
        self.ininterv(0.,24.,ndmean)
      coordmean = self.xcoord
      self.xcoord = sav
      return coordmean

    def meanperform(self,coordmean,umst=None):
      ndmean = coordmean.size
      meanpres = 0. ; meandens = 0. ; meantemp = 0. ; meanzonwind = 0. ; meanmerwind = 0. ; meanmeanvar = np.zeros(5) ; meanextvar = np.zeros(100)        
      for ccc in coordmean:
        if self.averaging == "lon":
          # zonal averaging with forcing of local time
          self.lon = ccc
          self.loct = (umst + self.lon/15.) % 24 #fixedlt false for this case
        elif self.averaging == "loct":
          self.loct = ccc
        self.update() 
        meanpres = meanpres + self.pres/float(ndmean) ; meandens = meandens + self.dens/float(ndmean) ; meantemp = meantemp + self.temp/float(ndmean)
        meanzonwind = meanzonwind + self.zonwind/float(ndmean) ; meanmerwind = meanmerwind + self.merwind/float(ndmean)
        meanmeanvar = meanmeanvar + self.meanvar/float(ndmean) ; meanextvar = meanextvar + self.extvar/float(ndmean)
      self.pres=meanpres ; self.dens=meandens ; self.temp=meantemp ; self.zonwind=meanzonwind ; self.merwind=meanmerwind
      self.meanvar=meanmeanvar ; self.extvar=meanextvar

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
      if incwind:   mcdcomp.maplatlon(self.xcoord,self.ycoord,field,title=fieldlab,proj=proj,vecx=windx,vecy=windy,vmin=self.min2d,vmax=self.max2d) #,stride=1)
      else:         mcdcomp.maplatlon(self.xcoord,self.ycoord,field,title=fieldlab,proj=proj,vmin=self.min2d,vmax=self.max2d)
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

    def htmlmap2d(self,tabtodo,incwind=False,figname="temp.png"):
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
        from scipy.io import netcdf
        filename = dataloc+"/mola32.nc" ; back = "alt"
        zefile = netcdf.netcdf_file(filename, 'r') 
        fieldc = zefile.variables[back][::32,::32]/1000.
        yc = zefile.variables['latitude'][::32]
        xc = zefile.variables['longitude'][::32]
        zefile.close()
        havetopo = True
      except:
        print "Trouble with netCDF or surface.nc file. Continue without topo lines."
        havetopo = False

      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.

      fig = mcdcomp.setfig(len(tabtodo),proj=self.proj)
      subv,subh = mcdcomp.definesubplot( len(tabtodo) , fig )

      for i in range(len(tabtodo)):
        yeah = fig.add_subplot(subv,subh,i+1)
        choice = tabtodo[i]
        self.latlon() #ndx=64,ndy=48) 
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
          if self.palt is None: self.palt  = 99999999.
          
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
          elif self.proj == "nsper": yeah = Basemap(projection=self.proj,lat_0=self.plat,lon_0=self.plon,satellite_height=self.palt*1000.,resolution=None,ax=yeah)
          
          ## NB: resolution=None is here to avoid loading coastlines which caused problems with some (plat,plon) couples
          ### background
          if isback:
            img="/home/marshttp/www-mars/mcd_python/MarsMap_2500x1250.jpg"
            yeah.warpimage(img,scale=0.75)
          mertab,partab = np.r_[-180.:180.:30.],np.r_[-90.:90.:15.]
          merlab,parlab = [0,0,0,1],[1,0,0,0]
          #format = '%.1f'
          
	  if(self.proj=="nsper") : 
           yeah.drawmeridians(mertab,color="grey")
           yeah.drawparallels(partab,color="grey")
	  else :
           yeah.drawmeridians(mertab,color="grey",labels=merlab)
           yeah.drawparallels(partab,color="grey",labels=parlab)	    
	  
          [lon2d,lat2d] = np.meshgrid(lon,lat)
          x, y = yeah(lon2d, lat2d)
        else:
          x = lon ; y = lat

        ## define field. bound field.
        what_I_plot = np.transpose(field)
        zevmin,zevmax,limtype = mcdcomp.setbounds(what_I_plot,vmin=self.min2d,vmax=self.max2d)
        ## define contour field levels. define color palette
        ticks = ndiv + 1
        zelevels = np.linspace(zevmin,zevmax,ticks)
        palette = get_cmap(name=self.colorm)

        ## topography contours
        if havetopo:
          rcParams['contour.negative_linestyle'] = 'solid' # negative contours solid instead of dashed
          zelevc = np.linspace(-9.,20.,11.,0.)
          if isproj:
             [xc2,yc2] = np.meshgrid(xc,yc)
             xc3,yc3 = yeah(xc2,yc2)
             yeah.contour( xc3, yc3, fieldc, zelevc, colors='black',linewidths = 0.4 )
             [xc2,yc2] = np.meshgrid(np.array(xc)-360.,yc)
             xc3,yc3 = yeah(xc2,yc2)
             yeah.contour( xc3, yc3, fieldc, zelevc, colors='black',linewidths = 0.4 )
          else:
             yeah.contour( np.array(xc)       , yc, fieldc, zelevc, colors='black',linewidths = 0.4)
             yeah.contour( np.array(xc) - 360., yc, fieldc, zelevc, colors='black',linewidths = 0.4)

        # contour field
        if self.iscontour: c = yeah.contour( x, y, what_I_plot, zelevels, cmap = palette )
        else: c = yeah.contourf( x, y, what_I_plot, zelevels, cmap = palette, alpha = 1.-self.trans, extend=limtype )

        # colorbar
        if not isproj:               orientation='vertical'   ; frac = 0.15  ; pad = 0.04 ; lu = 0.5
        elif self.proj in ['moll']:  orientation="horizontal" ; frac = 0.08  ; pad = 0.03 ; lu = 1.0
        elif self.proj in ['robin']: orientation="horizontal" ; frac = 0.07  ; pad = 0.1  ; lu = 1.0
        elif self.proj in ['cyl']:   orientation="vertical"   ; frac = 0.023 ; pad = 0.03 ; lu = 0.5
        else:                        orientation='vertical'   ; frac = 0.05  ; pad = 0.03 ; lu = 0.5
        if np.abs(zevmax) == 1e-35: zevmax=0.
        if np.abs(zevmin) == 1e-35: zevmin=0.
        zelevpal = np.linspace(zevmin,zevmax,num=min([ticks/2+1,21]))
        clb = Figure.colorbar(fig,c,orientation=orientation,format=self.fmt,ticks=zelevpal,\
             fraction=frac,pad=pad,spacing='proportional')
        #clb.set_label(fieldlab)

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

        # operation on axis: title, labels, ticks, etc.
        ax = fig.gca()
        rcParams['axes.titlesize'] = rcParams['font.size']
        ax.set_title(fieldlab)
        if not isproj:
          ax.set_ylabel("Latitude") ; ax.set_xlabel("Longitude")
          # make intervals 
          self.makeinterv()
          ax.set_xticks(np.arange(-360,361,self.loninterv)) ; ax.set_xbound(lower=self.lons, upper=self.lone)
          ax.set_yticks(np.arange(-90,91,self.latinterv)) ; ax.set_ybound(lower=self.lats, upper=self.late)
          ax.set_xlim(xmin=-180.,xmax=180.)

      ## titles and final production
      self.gettitle()
      #ax.set_title(self.title,x=0.5,y=1.05)
      #ax.set_xlabel('\n'+self.ack,x=0.5,y=0.05)

      fig.text(0.5, 0.95, self.title, ha='center', transform=fig.gca().transAxes, fontweight='bold')
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
      from matplotlib import rcParams
      from matplotlib.figure import Figure
      from matplotlib.backends.backend_agg import FigureCanvasAgg
      from matplotlib.cm import get_cmap
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.

      fig = mcdcomp.setfig(len(tabtodo))
      subv,subh = mcdcomp.definesubplot( len(tabtodo) , fig )

      #######################
      self.query2d(self)
      #######################

      for i in range(len(tabtodo)):
        yeah = fig.add_subplot(subv,subh,i+1)
        choice = tabtodo[i]

        (field, fieldlab) = self.definefield(choice)

        colorb=self.colorm ; ndiv=20 

        ## define field. bound field.
        what_I_plot = np.transpose(field)
        zevmin,zevmax,limtype = mcdcomp.setbounds(what_I_plot,vmin=self.min2d,vmax=self.max2d)  
        ## define contour field levels. define color palette
        ticks = ndiv + 1
        zelevels = np.linspace(zevmin,zevmax,ticks)
        palette = get_cmap(name=colorb)
        # contour field
        c = yeah.contourf( self.xcoord, self.ycoord, what_I_plot, zelevels, cmap = palette, extend=limtype )
        if np.abs(zevmax) == 1e-35: zevmax=0.
        if np.abs(zevmin) == 1e-35: zevmin=0.
        clb = Figure.colorbar(fig,c,orientation='vertical',format=self.fmt,ticks=np.linspace(zevmin,zevmax,num=min([ticks/2+1,21])))
        #clb.set_label(fieldlab)
        ax = fig.gca() ; ax.set_ylabel(self.ylabel) ; ax.set_xlabel(self.xlabel)
        rcParams['axes.titlesize'] = rcParams['font.size']
        ax.set_title(fieldlab)

        self.makeinterv()
        if self.lons is not None:
          if self.xdates is not None:
            ax.set_yticks(np.arange(-360,361,self.loninterv)) ; ax.set_ybound(lower=self.lons, upper=self.lone)
          else:
            ax.set_xticks(np.arange(-360,361,self.loninterv)) ; ax.set_xbound(lower=self.lons, upper=self.lone)
        elif self.lats is not None: 
          if self.xdates is not None:
            ax.set_yticks(np.arange(-90,91,self.latinterv)) ; ax.set_ybound(lower=self.lats, upper=self.late)
          else:
            ax.set_xticks(np.arange(-90,91,self.latinterv)) ; ax.set_xbound(lower=self.lats, upper=self.late)

        if self.locts is not None: 
          if self.xzs is not None: 
            ax.set_xticks(np.arange(0,26,2)) ; ax.set_xbound(lower=self.locts, upper=self.locte)
          else:
            ax.set_yticks(np.arange(0,26,2)) ; ax.set_ybound(lower=self.locts, upper=self.locte)
        elif self.xdates is not None:
            ax.set_xticks(np.arange(0,361,30)) ; ax.set_xbound(lower=self.xdates, upper=self.xdatee)

        if self.zkey == 4 and self.xzs is not None: 
            ax.set_yscale('log') ; ax.set_ylim(ax.get_ylim()[::-1])
        else:
            #ax.set_yticks(np.arange(self.xzs,self.xze,10000.)) ; 
            ax.set_ybound(lower=self.xzs, upper=self.xze)

      self.gettitle()
      fig.text(0.5, 0.95, self.title, ha='center', transform=fig.gca().transAxes, fontweight='bold')
      fig.text(0.5, 0.01, self.ack, ha='center')
      canvas = FigureCanvasAgg(fig)
      # The size * the dpi gives the final image size
      #   a4"x4" image * 80 dpi ==> 320x320 pixel image
      canvas.print_figure(figname, dpi=self.dpi)

    ### TODO: makeplot2d, plot2d, passer plot settings

