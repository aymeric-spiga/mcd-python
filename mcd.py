#####################################################
### A Python Class for the Venus Climate Database ###
### ----------------------------------------------###
### Aymeric SPIGA 17-21/04/2012                   ###
### T. PIERRON & A. BIERJON 10/2021               ###
### ----------------------------------------------###
### (see quicktest.py for examples of use)        ###
#####################################################

###
from fvcd import vcd # VCD compiled with f2py
from fvcd import dataloc  # location of VCD data file
from fvcd import dataver  # compiled version of VCD 
###
import numpy as np
import matplotlib.pyplot as mpl
###

## default number of points for each dimension
dfzon = 96 #37 # zonal dimension
dfmer = 96 #19 # meridional dimension
dfver = 35 #20 # vertical dimension
dflct = 25 #13 # local time dimension
dfsea = 25 # solar long dimension
## default names
lslab = "Areocentric longitude (degrees)"
latlab = "North latitude (degrees)"
lonlab = "East longitude (degrees)"
ltlab = "Local time (Venusian hours)"
#NB: vertical labels are treated by .vertlabel()

def errormess(text,printvar=None):
    print text
    if printvar is not None: print printvar
    exit()
    return

class vcd_class():
 
    def __repr__(self):
    # print out a help string when help is invoked on the object
        whatprint = 'VCD object. \"help(mcd)\" for more information\n'
        return whatprint

########################
### Default settings ###
########################

    def __init__(self):
    # default settings
        ## 0. general stuff
        self.name      = "VCD_v"+dataver # this is coming from fvcd
        self.dset      = dataloc # this is coming from fvcd
        self.ack       = "Venus Climate Database (c) LMD/ESA"
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
                            # 0 = pressure level (Pa)
                            # 1 = radius from centre of planet (m)
                            # 2 = height above Venus ref sphere (m) (sphere of radius 6051848m)
                            # 3 = height above surface (m)
        self.datekey   = 1  # 0 = "Earth time": xdate is given in Julian days (localtime must be set to zero)
                            # 1 = "Venus localtime": xdate is the value of local time
        ## 2. climatological options
        self.scena     = 1 # climatological average scenario
        self.varE107   = 0 # value for the solar EUV (varE107=0 if scena!=7)
        self.hrkey     = 1 # set high resolution mode on (hrkey=0 to set high resolution off)
        ## 3. additional settings for advanced use
        self.extvarkey = np.ones(100) #extra output variables flag table (1: yes, 0: no)
        self.perturkey = 0  #integer perturkey ! perturbation type (0: none)
        self.seedin    = 1  #random number generator seed (unused if perturkey=0)
        self.gwlength  = 0. #gravity Wave wavelength (unused if perturkey=0)
        ## outputs. just to define attributes.
        ## --> in update
        self.zonwind = None ; self.merwind = None ; self.vertwind = None ; self.temp = None ; self.pres = None ; self.dens = None ; self.extvar = None
        self.seedout = None ; self.ierr = None
        ## --> in prepare
        self.xcoord = None ; self.ycoord = None
        self.zonwindtab = None ; self.merwindtab = None ; self.vertwindtab = None
        self.temptab = None ; self.prestab = None ; self.denstab = None ; self.extvartab = None
        ## plot stuff
        self.typex = None ; self.typey = None
        self.xlabel = None ; self.ylabel = None ; self.title = ""
        self.vertplot = False
        self.fmt = "%.1e" 
        self.colorm = "jet"
        self.fixedlt = False # can be used as Universal Venus Time (uvt), or as Local True Solar Time
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
        self.latpoint = None
        self.lonpoint = None

    def getscenalabel(self):
        if self.scena == 1: self.scenalabel = "Standard cloud albedo scenario, average solar EUV conditions"
        elif self.scena == 2: self.scenalabel = "Standard cloud albedo scenario, minimum solar EUV conditions"
        elif self.scena == 3: self.scenalabel = "Standard cloud albedo scenario, maximum solar EUV conditions"
        elif self.scena == 4: self.scenalabel = "Low cloud albedo scenario, average solar EUV conditions"
        elif self.scena == 5: self.scenalabel = "High cloud albedo scenario, average solar EUV conditions"
        elif self.scena == 6: self.scenalabel = "EUV input as deduced from the input Julian date" #+ " ("+str(self.varE107)+" s.f.u)"
        elif self.scena == 7: self.scenalabel = "EUV input specified by user ("+str(int(self.varE107))+"sfu)"

    def gettitle(self,oneline=False):
        self.getscenalabel()
        self.title = self.name + " with " + self.scenalabel + "."
        #if self.datekey == 1:    
          #if self.xdates is None:
          # self.title = self.title + " LT " + str(self.xdate) + "Vhrs"
          #if self.locts is None and self.averaging is None:
          #  if self.lons is not None: # if longitude is a free dimension
          #   if not self.fixedlt:  self.title = self.title + " (at longitude 0) "
          #     else: self.title = self.title + " (fixed at all longitudes) "
        if self.datekey == 0:  
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
        if self.datekey == 0:  
          self.title = self.title + " JD " + str(self.xdate) + "."
        elif self.datekey == 1:
          if self.averaging == "loct":
            self.title = self.title + " Diurnal mean over all local times."
          else:
            if self.locts is None and self.averaging != "lon":
              self.title = self.title + " Local time " + str(self.loct) + "Vhrs"
              if self.lons is not None: # if longitude is a free dimension
                if not self.fixedlt:  self.title = self.title + " (at longitude 0) "
                else: self.title = self.title + " (fixed at all longitudes) "

    def getextvarlab(self,num):
        ## we use the end of extvar (unused) to store the meanvar (see update)
        whichfield = { \
	93: "Horizontal wind speed (m/s)", \
        94: "W-E wind component (m/s)", \
        95: "S-N wind component (m/s)", \
        96: "Upward vertical wind component (m/s)", \
        97: "Temperature (K)", \
        98: "Pressure (Pa)", \
        99: "Density (kg/m3)", \
        1: "distance to planet centre (m)", \
        2: "altitude above the reference sphere (m)", \
        3: "altitude above local surface (m)", \
        4: "orographic height (m) (altitude of the surface above the reference sphere)", \
        5: "orographic height (m) in the GCM", \
        6: "Solar longitude of Venus (deg)", \
        7: "Sun-Venus distance (UA)", \
        8: "Local true solar time at longitude lon (Vhrs)", \
        9: "Universal solar time (i.e. local true solar time at longitude 0) (Vhrs)", \
        10: "Solar zenith angle (deg)", \
        11: "unperturbed (climatological) zonal wind (m/s)", \
        12: "unperturbed (climatological) meridional wind (m/s)", \
        13: "unperturbed (climatological) vertical wind (m/s)", \
        14: "unperturbed (climatological) atmospheric temperature (K) ", \
        15: "unperturbed (climatological) atmospheric pressure (Pa)", \
        16: "unperturbed (climatological) atmospheric density (kg/m3)", \
        17: "Surface temperature (K)", \
        18: "Surface pressure (Pa)", \
        20: "V-hourly variability (RMS) of the zonal wind (m/s)", \
        21: "V-hourly variability (RMS) of the meridional wind (m/s)", \
        22: "V-hourly variability (RMS) of the vertical wind (m/s)", \
        23: "V-hourly variability (RMS) of the atmospheric temperature (K)", \
        24: "V-hourly variability (RMS) of the atmospheric pressure (Pa)", \
        25: "V-hourly variability (RMS) of the atmospheric density (kg/m3)", \
        26: "V-hourly variability (RMS) of the surface temperature (K)", \
        27: "V-hourly variability (RMS) of the surface pressure (Pa)", \
        30: "Vday to Vday variability (RMS) of the zonal wind (m/s)", \
        31: "Vday to Vday variability (RMS) of the meridional wind (m/s)", \
        32: "Vday to Vday variability (RMS) of the vertical wind (m/s)", \
        33: "Vday to Vday variability (RMS) of the atmospheric temperature (K)", \
        34: "Vday to Vday variability (RMS) of the atmospheric pressure (Pa)", \
        35: "Vday to Vday variability (RMS) of the atmospheric density (kg/m3)", \
        36: "Vday to Vday variability (RMS) of the surface temperature (K)", \
        37: "Vday to Vday variability (RMS) of the surface pressure (Pa)", \
        40: "Net Solar Flux (SW) received at the top of the atmosphere (W/m2), positive downward", \
        41: "SW net flux at given altitude (W/m2), positive downward", \
        42: "LW net flux at given altitude (W/m2), positive upward", \
        43: "atmospheric scale height at given altitude (m)", \
        44: "atmospheric mean molar mass at given altitude (g/mol)", \
        45: "atmospheric speed of sound cs (m/s)", \
        46: "atmospheric reduced molecular gas constant r (J/K/kg)", \
        47: "atmospheric heat capacity Cp (J/kg/K) ", \
        48: "atmospheric specific heat ratio $\gamma$ (N/A) ", \
        49: "atmospheric viscosity estimation (N/s/m2)", \
        50: "CO2 volume mixing ratio (mol/mol)", \
        51: "CO volume mixing ratio (mol/mol)", \
        52: "O2 volume mixing ratio (mol/mol)", \
        53: "O volume mixing ratio (mol/mol)", \
        54: "H volume mixing ratio (mol/mol)", \
        55: "H2 volume mixing ratio (mol/mol)", \
        56: "H2O volume mixing ratio (mol/mol)", \
        57: "SO2 volume mixing ratio (mol/mol)", \
        58: "SO volume mixing ratio (mol/mol)", \
        59: "OCS volume mixing ratio (mol/mol)", \
        60: "O3 volume mixing ratio (mol/mol)", \
        61: "HCl volume mixing ratio (mol/mol)", \
        62: "N2 volume mixing ratio (mol/mol)", \
        63: "He volume mixing ratio (mol/mol)", \
        70: "CO2 column (kg/m2)", \
        71: "CO column (kg/m2)", \
        72: "O2 column (kg/m2)", \
        73: "O column (kg/m2)", \
        74: "H column (kg/m2)", \
        75: "H2 column (kg/m2)", \
        76: "H2O column (kg/m2)", \
        77: "SO2 column (kg/m2)", \
        78: "SO column (kg/m2)", \
        79: "OCS column (kg/m2)", \
        80: "O3 column (kg/m2)", \
        81: "HCl column (kg/m2)", \
        82: "N2 column (kg/m2)", \
        83: "He column (kg/m2)", \
        84: "Vertically integrated O2 nightglow ($\Delta$ emission) (MR)", \
        90: "VIRA temperature (K) at the same location", \
        91: "VIRA pressure (Pa) at the same location", \
        92: "VIRA density (kg/m3) at the same location"
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
        if num == "u": num = 94
        elif num == "v": num = 95
        elif num == "w": num = 96
        elif num == "wind": num = 93
        elif num == "t": num = 97
        elif num == "p": num = 98
        elif num == "rho": num = 99
        
        elif num == "zradius": num = 1
        elif num == "zsphere": num = 2
        elif num == "zsurface": num = 3
        elif num == "oroheight": num = 4
        elif num == "oro_gcm": num = 5
        elif num == "ls": num = 6
        elif num == "venusau": num = 7
        elif num == "LTST": num = 8
        elif num == "utime": num = 9
        elif num == "solzenangle": num = 10
        
        elif num == "u_np": num = 11
        elif num == "v_np": num = 12
        elif num == "w_np": num = 13
        elif num == "t_np": num = 14
        elif num == "p_np": num = 15
        elif num == "rho_np": num = 16
        
        elif num == "tsurf": num = 17
        elif num == "ps": num = 18
        
        elif num == "rms_u": num = 20
        elif num == "rms_v": num = 21
        elif num == "rms_w": num = 22
        elif num == "rms_t": num = 23
        elif num == "rms_p": num = 24
        elif num == "rms_rho": num = 25
        elif num == "rms_tsurf": num = 26
        elif num == "rms_ps": num = 27
        
        elif num == "dtd_u": num = 30
        elif num == "dtd_v": num = 31
        elif num == "dtd_w": num = 32
        elif num == "dtd_t": num = 33
        elif num == "dtd_p": num = 34
        elif num == "dtd_rho": num = 35
        elif num == "dtd_tsurf": num = 36
        elif num == "dtd_ps": num = 37

        elif num == "solar_flux_top": num = 40
        elif num == "solar_flux_sw": num = 41
        elif num == "solar_flux_lw": num = 42
        elif num == "pscaleheight": num = 43
        elif num == "mean_mol_mass": num = 44
        
        elif num == "cs": num = 45
        elif num == "r_gas": num = 46
        elif num == "cp": num = 47
        elif num == "gamma": num = 48
        elif num == "viscosity": num = 49
        
        elif num == "vmr_co2": num = 50
        elif num == "vmr_co": num = 51
        elif num == "vmr_o2": num = 52
        elif num == "vmr_o": num = 53
        elif num == "vmr_h": num = 54
        elif num == "vmr_h2": num = 55
        elif num == "vmr_h2o": num = 56
        elif num == "vmr_so2": num = 57
        elif num == "vmr_so": num = 58
        elif num == "vmr_ocs": num = 59
        elif num == "vmr_o3": num = 60
        elif num == "vmr_hcl": num = 61
        elif num == "vmr_n2": num = 62
        elif num == "vmr_he": num = 63
        
        elif num == "col_co2": num = 70
        elif num == "col_co": num = 71
        elif num == "col_o2": num = 72
        elif num == "col_o": num = 73
        elif num == "col_h": num = 74
        elif num == "col_h2": num = 75
        elif num == "col_h2o": num = 76
        elif num == "col_so2": num = 77
        elif num == "col_so": num = 78
        elif num == "col_ocs": num = 79
        elif num == "col_o3": num = 80
        elif num == "col_hcl": num = 81
        elif num == "col_n2": num = 82
        elif num == "col_he": num = 83
        elif num == "o2_ng": num = 84

        elif num == "vira_temp": num = 90
        elif num == "vira_pres": num = 91
        elif num == "vira_dens": num = 92
        
        elif not isinstance(num, np.int): errormess("field reference not found.")
        return num

###################
### One request ###
###################

    def update(self):
    # retrieve fields from VCD (call_vcd). more info in fvcd.call_vcd.__doc__
        ## sanity first
        self.loct = abs(self.loct)%24
        if self.locts is not None and self.locte is not None: 
            self.locts = abs(self.locts)%24
            self.locte = abs(self.locte)%24
            if self.locts == self.locte: self.locte = self.locts + 24
        ## ensure that local time = 0 when using Earth dates
        if self.datekey == 0: self.loct = 0.
        ### now VCD request
        (self.zonwind, self.merwind, self.vertwind, self.temp, self.pres, self.dens, \
         self.extvar, self.seedout, self.ierr) \
         = \
         vcd.call_vcd(self.zkey,self.xz,self.lon,self.lat,self.hrkey, \
             self.datekey,self.xdate,self.loct,self.dset,self.scena,self.varE107, \
             self.perturkey,self.seedin,self.gwlength,self.extvarkey )
        ## we use the end of extvar (unused) to store meanvar. this is convenient for getextvar(lab)
        self.extvar[93] = self.zonwind ; self.extvar[94] = self.merwind
        self.extvar[95] = self.vertwind ; self.extvar[96] = self.temp
        self.extvar[97] = self.pres ; self.extvar[98] = self.dens
	self.extvar[92] = np.sqrt(self.extvar[93]**2 + self.extvar[94]**2)

        ## treat missing values 
        if self.temp == -999: self.extvar[:] = np.NaN

    def printset(self):
    # print main settings
        if self.datekey == 0: #Earth date
          print "zkey",self.zkey,"xz",self.xz,"lon",self.lon,"lat",self.lat,"hrkey",self.hrkey, \
                "xdate",self.xdate,"scena",self.scena
        else: #Venusian hour
          print "zkey",self.zkey,"xz",self.xz,"lon",self.lon,"lat",self.lat,"hrkey",self.hrkey, \
                "loct",self.loct,"scena",self.scena

    def getnameset(self):
    # set a name referring to settings [convenient for databases]
        strlat = str(self.lat)+str(self.lats)+str(self.late)
        strlon = str(self.lon)+str(self.lons)+str(self.lone)
        strxz = str(self.xz)+str(self.xzs)+str(self.xze)
        if self.datekey == 0: #Earth date
          strdate = str(self.xdate)+str(self.xdates)+str(self.xdatee)
          name = str(self.zkey)+strxz+strlon+strlat+str(self.hrkey)+str(self.datekey)+strdate+str(self.scena)
        else: #Venusian hour
          strloct = str(self.loct)+str(self.locts)+str(self.locte)
          name = str(self.zkey)+strxz+strlon+strlat+str(self.hrkey)+str(self.datekey)+strloct+str(self.scena)
        return name

    def printcoord(self):
    # print requested space-time coordinates
        if self.datekey == 0: #Earth date
          print "LAT",self.lat,"LON",self.lon,"XDATE",self.xdate
        else: #Venusian hour
          print "LAT",self.lat,"LON",self.lon,"LOCT",self.loct

    def printmeanvar(self):
    # print mean VCD variables
        print "Zonal wind = %5.3f meters per second." % (self.zonwind)
        print "Meridional wind = %5.3f meters per second." % (self.merwind)
        print "Vertical wind = %5.3f meters per second." % (self.vertwind)
        print "Temperature = %3.0f kelvins (%4.0f degrees celsius)." % (self.temp,self.temp-273.15)
        print "Pressure = %5.3f pascals. " % (self.pres)
        print "Density = %5.3f kilograms per cubic meter. " % (self.dens)
        
    def printextvar(self,num):
    # print extra VCD variables
        num = self.convertlab(num)
        dastr = str(self.extvar[num-1])
        if dastr == "nan":   print "!!!! There is a problem, probably a value is requested below the surface !!!!"
        else:                print self.getextvarlab(num) + " ..... " + dastr

    def printallextvar(self):
    # print all extra VCD variables    
        limit=92
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
    # 1. call VCD 2. print settings 3. print mean vars
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
      self.zonwindtab = np.ones(dashape) ; self.merwindtab = np.ones(dashape) ; self.vertwindtab = np.ones(dashape)
      self.temptab = np.ones(dashape) ; self.prestab = np.ones(dashape) ; self.denstab = np.ones(dashape)
      self.extvartab = np.ones(dashapeext)

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
      if self.zkey != 0 or not vertcoord:   tabtab = np.linspace(first,second,nd)
      else:                                 tabtab = np.logspace(first,second,nd)
      if not yaxis:      self.xcoord = tabtab
      else:              self.ycoord = tabtab

    def correctbounds(self,start,end,vertcoord):
      if self.zkey != 0 or not vertcoord:
        # regular altitudes
        if start > end: first = end ; second = start
        else:           first = start ; second = end
      else:
        # pressure: reversed avis
        if start < end: first = np.log10(end) ; second = np.log10(start)
        else:           first = np.log10(start) ; second = np.log10(end)
      return first, second

    def vertlabel(self):
      if self.zkey == 0:   self.xlabel = "pressure (Pa)"
      elif self.zkey == 1: self.xlabel = "radius from centre of planet (m)"
      elif self.zkey == 2: self.xlabel = "altitude above Venus ref sphere (m)"
      elif self.zkey == 3: self.xlabel = "height above surface (m)"

    def vertunits(self):
      if self.zkey == 0:   self.vunits = "Pa"
      elif self.zkey == 1: self.vunits = "m CP"  #Centre of Planet
      elif self.zkey == 2: self.vunits = "m AVS" #Above Venus Sphere
      elif self.zkey == 3: self.vunits = "m ALS" #Above Local Surface

    def vertaxis(self,number,yaxis=False):
#------MARS VALUES-------------    
#      if self.zkey == 0:   self.ininterv(1000.,0.001,number,start=self.xzs,end=self.xze,yaxis=yaxis,vertcoord=True)
#      elif self.zkey == 1: self.ininterv(3396000,3596000,number,start=self.xzs,end=self.xze,yaxis=yaxis,vertcoord=True)
#      elif self.zkey == 2: self.ininterv(-5000.,100000.,number,start=self.xzs,end=self.xze,yaxis=yaxis,vertcoord=True)
#      elif self.zkey == 3: self.ininterv(0.,120000.,number,start=self.xzs,end=self.xze,yaxis=yaxis,vertcoord=True)
      if self.zkey == 0:   self.ininterv(1.e7,1.0,number,start=self.xzs,end=self.xze,yaxis=yaxis,vertcoord=True) #~0-150km, no thermosphere
      elif self.zkey == 1: self.ininterv(6051848,6301848,number,start=self.xzs,end=self.xze,yaxis=yaxis,vertcoord=True) #~0-250km, with thermosphere
      elif self.zkey == 2: self.ininterv(-2000.,145000.,number,start=self.xzs,end=self.xze,yaxis=yaxis,vertcoord=True) #-2 - 145km, no thermosphere
      elif self.zkey == 3: self.ininterv(0.,150000.,number,start=self.xzs,end=self.xze,yaxis=yaxis,vertcoord=True) # 0-150km, no thermosphere

###################
### 1D analysis ###
###################

    def put1d(self,i):
    ## fill in subscript i in output arrays
    ## (arrays must have been correctly defined through prepare)
      if self.prestab is None:  errormess("arrays must be prepared first through self.prepare")
      self.zonwindtab[i] = self.zonwind ; self.merwindtab[i] = self.merwind ; self.vertwindtab[i] = self.vertwind
      self.temptab[i] = self.temp ; self.prestab[i] = self.pres ; self.denstab[i] = self.dens
      self.extvartab[i,1:100] = self.extvar[0:99] ## note: var numbering according to VCD manual is kept

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
      if not self.fixedlt: uvt = self.loct
      for i in range(nd): 
          self.lon = self.xcoord[i]
          if not self.fixedlt: self.loct = (uvt - self.lon/15.) % 24  # Venus has retrograde rotation
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
      if self.zkey == 0: mpl.semilogy() ; ax = mpl.gca() ; ax.set_ylim(ax.get_ylim()[::-1])
      if not self.vertplot and self.islog: mpl.semilogy()
      if self.vertplot and self.islog: mpl.semilogx()
      mpl.figtext(0.5, 0.01, self.ack, ha='center')

    def plot1d(self,tabtodo):
    ### complete 1D figure with possible multiplots
      import vcdcomp as vcdcomp
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      fig = mpl.figure() ; subv,subh = vcdcomp.definesubplot( len(tabtodo) , fig ) 
      for i in range(len(tabtodo)): mpl.subplot(subv,subh,i+1).grid(True, linestyle=':', color='grey') ; self.makeplot1d(tabtodo[i])
      mpl.show()

    def htmlplot1d(self,tabtodo,figname="temp.png",title=""):
    ### complete 1D figure with possible multiplots
    ### added in 09/2012 for online MCD
    ### see http://www.dalkescientific.com/writings/diary/archive/2005/04/23/matplotlib_without_gui.html
      import vcdcomp as vcdcomp
      from matplotlib.figure import Figure
      from matplotlib.backends.backend_agg import FigureCanvasAgg
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.

      howmanyplots = len(tabtodo)
      if howmanyplots == 1: fig = Figure(figsize=(16,8))
      elif howmanyplots == 2: fig = Figure(figsize=(8,8))
      elif howmanyplots == 3: fig = Figure(figsize=(8,16))
      elif howmanyplots == 4: fig = Figure(figsize=(16,8))

      subv,subh = vcdcomp.definesubplot( len(tabtodo) , fig )
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

        if self.xzs is not None and self.zkey == 0: ax.set_yscale('log') ; ax.set_ylim(ax.get_ylim()[::-1])
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
            self.typex = zetypet ; self.typey = zetypes # hovmoller
        elif zetypet == "ls":
          self.typex = zetypet ; self.typey = zetypes # seasonal
      return

    def query2d(self,typex=None,typey=None,ndmean=32):
    ### retrieve a 2D slice
      # save query
      save1 = self.lon ; save2 = self.xz ; save3 = self.loct ; save4 = self.lat ; save5 = self.xdate
      # initialize
      if not self.fixedlt: 
        uvt = self.loct # local time is not fixed. user-defined local time is at longitude 0.
      # define the type of 2D plot
      # -- hard setting of typex and typey is meant to disappear
      self.definetype(typex=typex,typey=typey)
      # settings for averaging
      if self.averaging is not None:
        coordmean = self.meandim(ndmean=ndmean) 
        if self.averaging == "lon": 
          self.fixedlt = False
          uvt = self.loct # shouldn't it be zero? does not matter anyways...
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
             self.loct = (uvt - self.lon/15.) % 24  # Venus has retrograde rotation
         # get field (simple or average)
         if self.averaging is None:
           self.update()
         else:
           self.meanperform(coordmean,uvt=uvt) 
         # fill in 2D array
         self.put2d(i,j)
      # reinstall init state
      if not self.fixedlt: self.loct = uvt # AB: why do this since it is changed just below?
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

    def meanperform(self,coordmean,uvt=None):
      ndmean = coordmean.size
      meanzonwind = 0. ; meanmerwind = 0. ; meanvertwind = 0. ; meantemp = 0. ; meanpres = 0. ; meandens = 0. ; meanextvar = np.zeros(100)        
      for ccc in coordmean:
        if self.averaging == "lon":
          # zonal averaging with forcing of local time
          self.lon = ccc
          self.loct = (uvt - self.lon/15.) % 24 #fixedlt false for this case  # Venus has retrograde rotation
        elif self.averaging == "loct":
          self.loct = ccc
        self.update() 
        meanzonwind = meanzonwind + self.zonwind/float(ndmean) ; meanmerwind = meanmerwind + self.merwind/float(ndmean) ; meanvertwind = meanvertwind + self.vertwind/float(ndmean)
        meantemp = meantemp + self.temp/float(ndmean) ; meanpres = meanpres + self.pres/float(ndmean) ; meandens = meandens + self.dens/float(ndmean)
        meanextvar = meanextvar + self.extvar/float(ndmean)
      self.zonwind=meanzonwind ; self.merwind=meanmerwind ; self.vertwind=meanvertwind
      self.temp=meantemp ; self.pres=meanpres ; self.dens=meandens
      self.extvar=meanextvar

    def put2d(self,i,j):
    ## fill in subscript i,j in output arrays
    ## (arrays must have been correctly defined through prepare)
      if self.prestab is None:  errormess("arrays must be prepared first through self.prepare")
      self.zonwindtab[i,j] = self.zonwind ; self.merwindtab[i,j] = self.merwind ; self.vertwindtab[i,j] = self.vertwind
      self.temptab[i,j] = self.temp ; self.prestab[i,j] = self.pres ; self.denstab[i,j] = self.dens
      self.extvartab[i,j,1:100] = self.extvar[0:99] ## note: var numbering according to VCD manual is kept

    def makemap2d(self,choice,incwind=False,proj="cyl"):
    ### one 2D map is created for the user-defined variable in choice.
      import vcdcomp as vcdcomp
      self.latlon() ## a map is implicitely a lat-lon plot. otherwise it is a plot (cf. makeplot2d)
      if choice == "wind" or incwind:
          (windx, fieldlabwx) = self.definefield("u")
          (windy, fieldlabwy) = self.definefield("v")
      if choice == "wind":
          field = np.sqrt(windx*windx + windy*windy)
	  print windx
          fieldlab = "Horizontal wind speed (m/s)"
      else:    
          (field, fieldlab) = self.definefield(choice)
      if incwind:   vcdcomp.maplatlon(self.xcoord,self.ycoord,field,title=fieldlab,proj=proj,vecx=windx,vecy=windy,vmin=self.min2d,vmax=self.max2d) #,stride=1)
      else:         vcdcomp.maplatlon(self.xcoord,self.ycoord,field,title=fieldlab,proj=proj,vmin=self.min2d,vmax=self.max2d)
      mpl.figtext(0.5, 0.0, self.ack, ha='center')

    def map2d(self,tabtodo,incwind=False,proj="cyl"):
    ### complete 2D figure with possible multiplots
      import vcdcomp as vcdcomp
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      fig = mpl.figure()
      subv,subh = vcdcomp.definesubplot( len(tabtodo) , fig ) 
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
      import vcdcomp as vcdcomp
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
        filename = dataloc+"/relief.nc" ; back = "relief"
        zefile = netcdf.netcdf_file(filename, 'r') 
        fieldc = zefile.variables[back][::32,::32]/1000.
        yc = zefile.variables['latitude'][::32]
        xc = zefile.variables['longitude'][::32]
        zefile.close()
        havetopo = True
      except:
        print "Trouble with netCDF or relief.nc file. Continue without topo lines."
        havetopo = False

      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.

      fig = vcdcomp.setfig(len(tabtodo),proj=self.proj)
      subv,subh = vcdcomp.definesubplot( len(tabtodo) , fig )

      for i in range(len(tabtodo)):
        yeah = fig.add_subplot(subv,subh,i+1)
        choice = tabtodo[i]
        self.latlon() #ndx=96,ndy=96) 
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
            img="/home/marshttp/www-venus/vcd_python/mapvenus.jpg"
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
        zevmin,zevmax,limtype = vcdcomp.setbounds(what_I_plot,vmin=self.min2d,vmax=self.max2d)
        ## define contour field levels. define color palette
        ticks = ndiv + 1
        zelevels = np.linspace(zevmin,zevmax,ticks)
        palette = get_cmap(name=self.colorm)

        ## topography contours
        if havetopo:
          rcParams['contour.negative_linestyle'] = 'solid' # negative contours solid instead of dashed
          zelevc = np.linspace(-9.,20.,11.,0.)
          if isproj:
             [xc2,yc2] = np.meshgrid(np.array(xc)+180.,yc)
             xc3,yc3 = yeah(xc2,yc2)
             yeah.contour( xc3, yc3, fieldc, zelevc, colors='black',linewidths = 0.4 )
             [xc2,yc2] = np.meshgrid(np.array(xc)-180.,yc)
             xc3,yc3 = yeah(xc2,yc2)
             yeah.contour( xc3, yc3, fieldc, zelevc, colors='black',linewidths = 0.4 )
          else:
             yeah.contour( np.array(xc) + 180 , yc, fieldc, zelevc, colors='black',linewidths = 0.4)
             yeah.contour( np.array(xc) - 180., yc, fieldc, zelevc, colors='black',linewidths = 0.4)

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
      import vcdcomp as vcdcomp
      from matplotlib import rcParams
      from matplotlib.figure import Figure
      from matplotlib.backends.backend_agg import FigureCanvasAgg
      from matplotlib.cm import get_cmap
      if isinstance(tabtodo,np.str): tabtodo=[tabtodo] ## so that asking one element without [] is possible.
      if isinstance(tabtodo,np.int): tabtodo=[tabtodo] ## so that asking one element without [] is possible.

      fig = vcdcomp.setfig(len(tabtodo))
      subv,subh = vcdcomp.definesubplot( len(tabtodo) , fig )

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
        zevmin,zevmax,limtype = vcdcomp.setbounds(what_I_plot,vmin=self.min2d,vmax=self.max2d)  
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
          if self.locts is not None:
            ax.set_yticks(np.arange(-360,361,self.loninterv)) ; ax.set_ybound(lower=self.lons, upper=self.lone)
          else:
            ax.set_xticks(np.arange(-360,361,self.loninterv)) ; ax.set_xbound(lower=self.lons, upper=self.lone)
        elif self.lats is not None: 
          if self.locts is not None:
            ax.set_yticks(np.arange(-90,91,self.latinterv)) ; ax.set_ybound(lower=self.lats, upper=self.late)
          else:
            ax.set_xticks(np.arange(-90,91,self.latinterv)) ; ax.set_xbound(lower=self.lats, upper=self.late)

        if self.locts is not None: 
            ax.set_xticks(np.arange(0,26,2)) ; ax.set_xbound(lower=self.locts, upper=self.locte)
        elif self.xdates is not None:
            ax.set_xticks(np.arange(0,361,30)) ; ax.set_xbound(lower=self.xdates, upper=self.xdatee)

        if self.zkey == 0 and self.xzs is not None: 
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

