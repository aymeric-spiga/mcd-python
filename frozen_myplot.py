## Author: AS
def errormess(text,printvar=None):
    print text
    if printvar is not None: print printvar 
    exit()
    return

## Author: AS
def adjust_length (tab, zelen):
    from numpy import ones
    if tab is None:
        outtab = ones(zelen) * -999999
    else:
        if zelen != len(tab):
            print "not enough or too much values... setting same values all variables"
            outtab = ones(zelen) * tab[0]
        else:
            outtab = tab
    return outtab

## Author: AS
def getname(var=False,var2=False,winds=False,anomaly=False):
    if var and winds:     basename = var + '_UV'
    elif var:             basename = var
    elif winds:           basename = 'UV'
    else:                 errormess("please set at least winds or var",printvar=nc.variables)
    if anomaly:           basename = 'd' + basename
    if var2:              basename = basename + '_' + var2
    return basename

## Author: AS + AC
def localtime(time,lon,namefile): # lon is the mean longitude of the plot, not of the domain. central lon of domain is taken from cen_lon
    import numpy as np
    from netCDF4 import Dataset
    ## THIS IS FOR MESOSCALE
    nc = Dataset(namefile)
    ## get start date and intervals
    dt_hour=1. ; start=0.
    if hasattr(nc,'TITLE'):
        title=getattr(nc, 'TITLE')
        if hasattr(nc,'DT') and hasattr(nc,'START_DATE') and 'MRAMS' in title: 
            ## we must adapt what is done in getlschar to MRAMS (outputs from ic.py)
            dt_hour=getattr(nc, 'DT')/60. 
            start_date=getattr(nc, 'START_DATE')
            start_hour=np.float(start_date[11:13])
            start_minute=np.float(start_date[14:16])/60.
            start=start_hour+start_minute # start is the local time of simu at longitude 0
            #LMD MMM is 1 output/hour (and not 1 output/timestep)
            #MRAMS is 1 output/timestep, unless stride is added in ic.py
        elif 'WRF' in title: 
            [dummy,start,dt_hour] = getlschar ( namefile ) # get start hour and interval hour
    ## get longitude 
    if lon is not None:
       if lon[0,1]!=lon[0,0]: mean_lon_plot = 0.5*(lon[0,1]-lon[0,0])
       else: mean_lon_plot=lon[0,0]
    elif hasattr(nc, 'CEN_LON'): mean_lon_plot=getattr(nc, 'CEN_LON')
    else: mean_lon_plot=0.
    ## calculate local time
    ltst = start + time*dt_hour + mean_lon_plot / 15.
    ltst = int (ltst * 10) / 10.
    ltst = ltst % 24
    return ltst

## Author: AC
def check_localtime(time):
    a=-1
    for i in range(len(time)-1):
       if (time[i] > time[i+1]): a=i
    if a >= 0 and a < (len(time)-1)/2.:
       print "Sorry, time axis is not regular."
       print "Contourf needs regular axis... recasting"
       for i in range(a+1):
          time[i]=time[i]-24.
    if a >= 0 and a >= (len(time)-1)/2.:
       print "Sorry, time axis is not regular."
       print "Contourf needs regular axis... recasting"
       for i in range((len(time)-1) - a):
          time[a+1+i]=time[a+1+i]+24.
    return time

## Author: AS, AC, JL
def whatkindfile (nc):
    typefile = 'gcm' # default
    if 'controle' in nc.variables:             typefile = 'gcm'
    elif 'phisinit' in nc.variables:           typefile = 'gcm'
    elif 'phis' in nc.variables:               typefile = 'gcm'
    elif 'time_counter' in nc.variables:       typefile = 'earthgcm'
    elif hasattr(nc,'START_DATE'):             typefile = 'meso' 
    elif 'HGT_M' in nc.variables:              typefile = 'geo'
    elif hasattr(nc,'institution'):
      if "European Centre" in getattr(nc,'institution'):  typefile = 'ecmwf'
    return typefile

## Author: AS
def getfield (nc,var):
    ## this allows to get much faster (than simply referring to nc.variables[var])
    import numpy as np
    dimension = len(nc.variables[var].dimensions)
    #print "   Opening variable with", dimension, "dimensions ..."
    if dimension == 2:    field = nc.variables[var][:,:]
    elif dimension == 3:  field = nc.variables[var][:,:,:]
    elif dimension == 4:  field = nc.variables[var][:,:,:,:]
    elif dimension == 1:  field = nc.variables[var][:]
    # if there are NaNs in the ncdf, they should be loaded as a masked array which will be
    # recasted as a regular array later in reducefield
    if (np.isnan(np.sum(field)) and (type(field).__name__ not in 'MaskedArray')):
       print "Warning: netcdf as nan values but is not loaded as a Masked Array."
       print "recasting array type"
       out=np.ma.masked_invalid(field)
       out.set_fill_value([np.NaN])
    else:
    # missing values from zrecast or hrecast are -1e-33
       masked=np.ma.masked_where(field < -1e30,field)
       masked2=np.ma.masked_where(field > 1e35,field)
       masked.set_fill_value([np.NaN]) ; masked2.set_fill_value([np.NaN])
       mask = np.ma.getmask(masked) ; mask2 = np.ma.getmask(masked2)
       if (True in np.array(mask)):
          out=masked
          print "Masked array... Missing value is NaN"
       elif (True in np.array(mask2)):
          out=masked2
          print "Masked array... Missing value is NaN"
#       else:
#       # missing values from api are 1e36
#          masked=np.ma.masked_where(field > 1e35,field)
#          masked.set_fill_value([np.NaN])
#          mask = np.ma.getmask(masked)
#          if (True in np.array(mask)):out=masked
#          else:out=field
       else:
#       # missing values from MRAMS files are 0.100E+32
          masked=np.ma.masked_where(field > 1e30,field)
          masked.set_fill_value([np.NaN])
          mask = np.ma.getmask(masked)
          if (True in np.array(mask)):out=masked
          else:out=field
#       else:out=field
    return out

## Author: AC
# Compute the norm of the winds or return an hodograph
# The corresponding variable to call is UV or uvmet (to use api)
def windamplitude (nc,mode):
    import numpy as np
    varinfile = nc.variables.keys()
    if "U" in varinfile: zu=getfield(nc,'U')
    elif "Um" in varinfile: zu=getfield(nc,'Um')
    else: errormess("you need slopex or U or Um in your file.")
    if "V" in varinfile: zv=getfield(nc,'V')
    elif "Vm" in varinfile: zv=getfield(nc,'Vm')
    else: errormess("you need V or Vm in your file.")
    znt,znz,zny,znx = np.array(zu).shape
    if hasattr(nc,'WEST-EAST_PATCH_END_UNSTAG'):znx=getattr(nc, 'WEST-EAST_PATCH_END_UNSTAG')
    zuint = np.zeros([znt,znz,zny,znx])
    zvint = np.zeros([znt,znz,zny,znx])
    if "U" in varinfile:
       if hasattr(nc,'SOUTH-NORTH_PATCH_END_STAG'): zny_stag=getattr(nc, 'SOUTH-NORTH_PATCH_END_STAG')
       if hasattr(nc,'WEST-EAST_PATCH_END_STAG'): znx_stag=getattr(nc, 'WEST-EAST_PATCH_END_STAG')
       if zny_stag == zny: zvint=zv
       else:
          for yy in np.arange(zny):      zvint[:,:,yy,:] = (zv[:,:,yy,:] + zv[:,:,yy+1,:])/2.
       if znx_stag == znx: zuint=zu
       else:
          for xx in np.arange(znx):      zuint[:,:,:,xx] = (zu[:,:,:,xx] + zu[:,:,:,xx+1])/2.
    else:
       zuint=zu
       zvint=zv
    if mode=='amplitude': return np.sqrt(zuint**2 + zvint**2)
    if mode=='hodograph': return zuint,zvint
    if mode=='hodograph_2': return None, 360.*np.arctan(zvint/zuint)/(2.*np.pi)

## Author: AC
# Compute the enrichment factor of non condensible gases
# The corresponding variable to call is enfact
# enrichment factor is computed as in Yuan Lian et al. 2012
# i.e. you need to have VL2 site at LS 135 in your data
# this only requires co2col so that you can concat.nc at low cost
def enrichment_factor(nc,lon,lat,time):
    import numpy as np
    from myplot import reducefield
    varinfile = nc.variables.keys()
    if "co2col" in varinfile: co2col=getfield(nc,'co2col')
    else: print "error, you need co2col var in your file"
    if "ps" in varinfile: ps=getfield(nc,'ps')
    else: print "error, you need ps var in your file"
    dimension = len(nc.variables['co2col'].dimensions)
    if dimension == 2: 
      zny,znx = np.array(co2col).shape
      znt=1
    elif dimension == 3: znt,zny,znx = np.array(co2col).shape
    mmrarcol = np.zeros([znt,zny,znx])
    enfact = np.zeros([znt,zny,znx])
    grav=3.72
    mmrarcol[:,:,:] = 1. - grav*co2col[:,:,:]/ps[:,:,:]
# Computation with reference argon mmr at VL2 Ls 135 (as in Yuan Lian et al 2012)
    lonvl2=np.zeros([1,2])
    latvl2=np.zeros([1,2])
    timevl2=np.zeros([1,2])
    lonvl2[0,0]=-180
    lonvl2[0,1]=180
    latvl2[:,:]=48.16
    timevl2[:,:]=135.
    indexlon  = getsindex(lonvl2,0,lon)
    indexlat  = getsindex(latvl2,0,lat)
    indextime = getsindex(timevl2,0,time)
    mmrvl2, error = reducefield( mmrarcol, d4=indextime, d1=indexlon, d2=indexlat)
    print "VL2 Ls 135 mmr arcol:", mmrvl2
    enfact[:,:,:] = mmrarcol[:,:,:]/mmrvl2
    return enfact

## Author: AC
# Compute the norm of the slope angles
# The corresponding variable to call is SLOPEXY
def slopeamplitude (nc):
    import numpy as np
    varinfile = nc.variables.keys()
    if "slopex" in varinfile: zu=getfield(nc,'slopex')
    elif "SLOPEX" in varinfile: zu=getfield(nc,'SLOPEX')
    else: errormess("you need slopex or SLOPEX in your file.") 
    if "slopey" in varinfile: zv=getfield(nc,'slopey')
    elif "SLOPEY" in varinfile: zv=getfield(nc,'SLOPEY')
    else: errormess("you need slopey or SLOPEY in your file.")
    znt,zny,znx = np.array(zu).shape
    zuint = np.zeros([znt,zny,znx])
    zvint = np.zeros([znt,zny,znx])
    zuint=zu
    zvint=zv
    return np.sqrt(zuint**2 + zvint**2)

## Author: AC
# Compute the temperature difference between surface and first level.
# API is automatically called to get TSURF and TK.
# The corresponding variable to call is DELTAT 
def deltat0t1 (nc):
    import numpy as np
    varinfile = nc.variables.keys()
    if "tsurf" in varinfile: zu=getfield(nc,'tsurf')
    elif "TSURF" in varinfile: zu=getfield(nc,'TSURF')
    else: errormess("You need tsurf or TSURF in your file")
    if "tk" in varinfile: zv=getfield(nc,'tk')
    elif "TK" in varinfile: zv=getfield(nc,'TK')
    else: errormess("You need tk or TK in your file. (might need to use API. try to add -i 4 -l XXX)")
    znt,zny,znx = np.array(zu).shape
    zuint = np.zeros([znt,zny,znx])
    zuint=zu - zv[:,0,:,:]
    return zuint

## Author: AS + TN + AC
def reducefield (input,d4=None,d3=None,d2=None,d1=None,yint=False,alt=None,anomaly=False,redope=None,mesharea=None,unidim=999):
    ### we do it the reverse way to be compliant with netcdf "t z y x" or "t y x" or "y x"
    ### it would be actually better to name d4 d3 d2 d1 as t z y x
    ### ... note, anomaly is only computed over d1 and d2 for the moment
    import numpy as np
    from frozen_mymath import max,mean,min,sum,getmask
    csmooth = 12 ## a fair amount of grid points (too high results in high computation time)
    if redope is not None:
       if   redope == "mint":     input = min(input,axis=0) ; d1 = None
       elif redope == "maxt":     input = max(input,axis=0) ; d1 = None
       elif redope == "edge_y1":  input = input[:,:,0,:]    ; d2 = None
       elif redope == "edge_y2":  input = input[:,:,-1,:]   ; d2 = None
       elif redope == "edge_x1":  input = input[:,:,:,0]    ; d1 = None
       elif redope == "edge_x2":  input = input[:,:,:,-1]   ; d1 = None
       else:                      errormess("not supported. but try lines in reducefield beforehand.")
       #elif redope == "minz":     input = min(input,axis=1) ; d2 = None
       #elif redope == "maxz":     input = max(input,axis=1) ; d2 = None
       #elif redope == "miny":     input = min(input,axis=2) ; d3 = None
       #elif redope == "maxy":     input = max(input,axis=2) ; d3 = None
       #elif redope == "minx":     input = min(input,axis=3) ; d4 = None
       #elif redope == "maxx":     input = max(input,axis=3) ; d4 = None
    dimension = np.array(input).ndim
    shape = np.array(np.array(input).shape)
    #print 'd1,d2,d3,d4: ',d1,d2,d3,d4
    if anomaly: print 'ANOMALY ANOMALY'
    output = input
    error = False
    #### this is needed to cope the case where d4,d3,d2,d1 are single integers and not arrays
    if d4 is not None and not isinstance(d4, np.ndarray): d4=[d4]
    if d3 is not None and not isinstance(d3, np.ndarray): d3=[d3]
    if d2 is not None and not isinstance(d2, np.ndarray): d2=[d2]
    if d1 is not None and not isinstance(d1, np.ndarray): d1=[d1]
    ### now the main part
    if dimension == 2:
	#### this is needed for 1d-type files (where dim=2 but axes are time-vert and not lat-lon)
        if unidim==1: d2=d4 ; d1=d3 ; d4 = None ; d3 = None
        if mesharea is None: mesharea=np.ones(shape)
        if   max(d2) >= shape[0]: error = True
        elif max(d1) >= shape[1]: error = True
        elif d1 is not None and d2 is not None:
          try:
            totalarea = np.ma.masked_where(getmask(output),mesharea)
            totalarea = mean(totalarea[d2,:],axis=0);totalarea = mean(totalarea[d1])
          except: print "(problem with areas. I skip this)" ; mesharea = 1. ; totalarea = 1.
          output = output*mesharea; output = mean(output[d2,:],axis=0); output = mean(output[d1])/totalarea
        elif d1 is not None:         output = mean(input[:,d1],axis=1)
        elif d2 is not None:
          try:
            totalarea = np.ma.masked_where(getmask(output),mesharea)
            totalarea = mean(totalarea[d2,:],axis=0)
          except: print "(problem with areas. I skip this)" ; mesharea = 1. ; totalarea = 1.
          output = output*mesharea; output = mean(output[d2,:],axis=0)/totalarea
    elif dimension == 3:
        if mesharea is None: mesharea=np.ones(shape[[1,2]])
        if   max(d4) >= shape[0]: error = True
        elif max(d2) >= shape[1]: error = True
        elif max(d1) >= shape[2]: error = True
        elif d4 is not None and d2 is not None and d1 is not None:
          output = mean(input[d4,:,:],axis=0)
          try:
            totalarea = np.ma.masked_where(getmask(output),mesharea)
            totalarea = mean(totalarea[d2,:],axis=0);totalarea = mean(totalarea[d1])
          except: print "(problem with areas. I skip this)" ; mesharea = 1. ; totalarea = 1.
          output = output*mesharea; output = mean(output[d2,:],axis=0); output = mean(output[d1])/totalarea
        elif d4 is not None and d2 is not None:
	  output = mean(input[d4,:,:],axis=0)
          try:
            totalarea = np.ma.masked_where(getmask(output),mesharea)
            totalarea = mean(totalarea[d2,:],axis=0)
          except: print "(problem with areas. I skip this)" ; mesharea = 1. ; totalarea = 1.
          output = output*mesharea; output = mean(output[d2,:],axis=0)/totalarea
        elif d4 is not None and d1 is not None:    output = mean(input[d4,:,:],axis=0); output=mean(output[:,d1],axis=1)
        elif d2 is not None and d1 is not None:
          try:
            totalarea = np.tile(mesharea,(output.shape[0],1,1))
            totalarea = np.ma.masked_where(getmask(output),totalarea)
            totalarea = mean(totalarea[:,d2,:],axis=1);totalarea = mean(totalarea[:,d1],axis=1)
          except: print "(problem with areas. I skip this)" ; mesharea = 1. ; totalarea = 1.
          output = output*mesharea; output = mean(output[:,d2,:],axis=1); output = mean(output[:,d1],axis=1)/totalarea
        elif d1 is not None:   output = mean(input[:,:,d1],axis=2)
        elif d2 is not None:   
          try:
            totalarea = np.tile(mesharea,(output.shape[0],1,1))
            totalarea = np.ma.masked_where(getmask(output),totalarea)
            totalarea = mean(totalarea[:,d2,:],axis=1)
          except: print "(problem with areas. I skip this)" ; mesharea = 1. ; totalarea = 1.
          output = output*mesharea; output = mean(output[:,d2,:],axis=1)/totalarea
        elif d4 is not None:   output = mean(input[d4,:,:],axis=0)
    elif dimension == 4:
        if mesharea is None: mesharea=np.ones(shape[[2,3]]) # mesharea=np.random.random_sample(shape[[2,3]])*5. + 2. # pour tester
        if   max(d4) >= shape[0]: error = True
        elif max(d3) >= shape[1]: error = True
        elif max(d2) >= shape[2]: error = True
        elif max(d1) >= shape[3]: error = True
        elif d4 is not None and d3 is not None and d2 is not None and d1 is not None:
             output = mean(input[d4,:,:,:],axis=0)
             output = reduce_zaxis(output[d3,:,:],ax=0,yint=yint,vert=alt,indice=d3)
             if anomaly: output = 100. * ((output / smooth(output,csmooth)) - 1.)
             try:
               totalarea = np.ma.masked_where(np.isnan(output),mesharea)
               totalarea = mean(totalarea[d2,:],axis=0); totalarea = mean(totalarea[d1])
             except: print "(problem with areas. I skip this)" ; mesharea = 1. ; totalarea = 1.
             output = output*mesharea; output = mean(output[d2,:],axis=0); output = mean(output[d1])/totalarea
        elif d4 is not None and d3 is not None and d2 is not None: 
             output = mean(input[d4,:,:,:],axis=0)
             output = reduce_zaxis(output[d3,:,:],ax=0,yint=yint,vert=alt,indice=d3)
             if anomaly: output = 100. * ((output / smooth(output,csmooth)) - 1.)
             try:
               totalarea = np.ma.masked_where(np.isnan(output),mesharea)
               totalarea = mean(totalarea[d2,:],axis=0)
             except: print "(problem with areas. I skip this)" ; mesharea = 1. ; totalarea = 1.
             output = output*mesharea; output = mean(output[d2,:],axis=0)/totalarea
        elif d4 is not None and d3 is not None and d1 is not None: 
             output = mean(input[d4,:,:,:],axis=0)
             output = reduce_zaxis(output[d3,:,:],ax=0,yint=yint,vert=alt,indice=d3)
             if anomaly: output = 100. * ((output / smooth(output,csmooth)) - 1.) 
             output = mean(output[:,d1],axis=1)
        elif d4 is not None and d2 is not None and d1 is not None: 
             output = mean(input[d4,:,:,:],axis=0)
             if anomaly:
                 for k in range(output.shape[0]):  output[k,:,:] = 100. * ((output[k,:,:] / smooth(output[k,:,:],csmooth)) - 1.)
             try:
               totalarea = np.tile(mesharea,(output.shape[0],1,1))
               totalarea = np.ma.masked_where(getmask(output),mesharea)
               totalarea = mean(totalarea[:,d2,:],axis=1); totalarea = mean(totalarea[:,d1],axis=1)
             except: print "(problem with areas. I skip this)" ; mesharea = 1. ; totalarea = 1.
             output = output*mesharea; output = mean(output[:,d2,:],axis=1); output = mean(output[:,d1],axis=1)/totalarea
             #noperturb = smooth1d(output,window_len=7)
             #lenlen = len(output) ; output = output[1:lenlen-7] ; yeye = noperturb[4:lenlen-4]
             #plot(output) ; plot(yeye) ; show() ; plot(output-yeye) ; show()
        elif d3 is not None and d2 is not None and d1 is not None:
             output = reduce_zaxis(input[:,d3,:,:],ax=1,yint=yint,vert=alt,indice=d3)
             if anomaly:
                 for k in range(output.shape[0]):  output[k,:,:] = 100. * ((output[k,:,:] / smooth(output[k,:,:],csmooth)) - 1.)
             try:
               totalarea = np.tile(mesharea,(output.shape[0],1,1))
               totalarea = np.ma.masked_where(getmask(output),totalarea)
               totalarea = mean(totalarea[:,d2,:],axis=1); totalarea = mean(totalarea[:,d1],axis=1)
             except: print "(problem with areas. I skip this)" ; mesharea = 1. ; totalarea = 1.
             output = output*mesharea; output = mean(output[:,d2,:],axis=1); output = mean(output[:,d1],axis=1)/totalarea
        elif d4 is not None and d3 is not None:  
             output = mean(input[d4,:,:,:],axis=0) 
             output = reduce_zaxis(output[d3,:,:],ax=0,yint=yint,vert=alt,indice=d3)
             if anomaly: output = 100. * ((output / smooth(output,csmooth)) - 1.) 
        elif d4 is not None and d2 is not None:  
             output = mean(input[d4,:,:,:],axis=0)
             try:
               totalarea = np.tile(mesharea,(output.shape[0],1,1))
               totalarea = np.ma.masked_where(getmask(output),mesharea)
               totalarea = mean(totalarea[:,d2,:],axis=1)
             except: print "(problem with areas. I skip this)" ; mesharea = 1. ; totalarea = 1.
             output = output*mesharea; output = mean(output[:,d2,:],axis=1)/totalarea
        elif d4 is not None and d1 is not None:  
             output = mean(input[d4,:,:,:],axis=0) 
             output = mean(output[:,:,d1],axis=2)
        elif d3 is not None and d2 is not None:
             output = reduce_zaxis(input[:,d3,:,:],ax=1,yint=yint,vert=alt,indice=d3)
             try:
               totalarea = np.tile(mesharea,(output.shape[0],1,1))
               totalarea = np.ma.masked_where(getmask(output),mesharea)
               totalarea = mean(totalarea[:,d2,:],axis=1)
             except: print "(problem with areas. I skip this)" ; mesharea = 1. ; totalarea = 1.
             output = output*mesharea; output = mean(output[:,d2,:],axis=1)/totalarea
        elif d3 is not None and d1 is not None:  
             output = reduce_zaxis(input[:,d3,:,:],ax=1,yint=yint,vert=alt,indice=d3)
             output = mean(output[:,:,d1],axis=2)
        elif d2 is not None and d1 is not None:
             try:
               totalarea = np.tile(mesharea,(output.shape[0],output.shape[1],1,1))
               totalarea = np.ma.masked_where(getmask(output),totalarea)
               totalarea = mean(totalarea[:,:,d2,:],axis=2); totalarea = mean(totalarea[:,:,d1],axis=1)
             except: print "(problem with areas. I skip this)" ; mesharea = 1. ; totalarea = 1.
             output = output*mesharea; output = mean(output[:,:,d2,:],axis=2); output = mean(output[:,:,d1],axis=2)/totalarea
        elif d1 is not None:        output = mean(input[:,:,:,d1],axis=3)
        elif d2 is not None:
             try:
               totalarea = np.tile(mesharea,(output.shape[0],output.shape[1],1,output.shape[3]))
               totalarea = np.ma.masked_where(getmask(output),totalarea)
               totalarea = mean(totalarea[:,:,d2,:],axis=2)
             except: print "(problem with areas. I skip this)" ; mesharea = 1. ; totalarea = 1.
             output = output*mesharea; output = mean(output[:,:,d2,:],axis=2)/totalarea
        elif d3 is not None:        output = reduce_zaxis(input[:,d3,:,:],ax=1,yint=yint,vert=alt,indice=d3)
        elif d4 is not None:        output = mean(input[d4,:,:,:],axis=0)
    dimension2 = np.array(output).ndim
    shape2 = np.array(output).shape
    print 'REDUCEFIELD dim,shape: ',dimension,shape,' >>> ',dimension2,shape2 
    return output, error

## Author: AC + AS
def reduce_zaxis (input,ax=None,yint=False,vert=None,indice=None):
    from frozen_mymath import max,mean
    from scipy import integrate
    import numpy as np
    if yint and vert is not None and indice is not None:
       if type(input).__name__=='MaskedArray':
          input.set_fill_value([np.NaN])
          output = integrate.trapz(input.filled(),x=vert[indice],axis=ax)
       else:
          output = integrate.trapz(input,x=vert[indice],axis=ax)
    else:
       output = mean(input,axis=ax)
    return output

## Author: AS + TN
def definesubplot ( numplot, fig, ipreferline=False):
    from matplotlib.pyplot import rcParams
    rcParams['font.size'] = 12. ## default (important for multiple calls)
    if numplot <= 0:
        subv = 99999
        subh = 99999
    elif numplot == 1: 
        subv = 1
        subh = 1
    elif numplot == 2:
        subv = 2 #1 #2
        subh = 1 #2 #1
        fig.subplots_adjust(wspace = 0.35)
        fig.subplots_adjust(hspace = 0.3)
        #rcParams['font.size'] = int( rcParams['font.size'] * 3. / 4. )
    elif numplot == 3:
        subv = 3
        subh = 1
        fig.subplots_adjust(hspace = 0.25)
        fig.subplots_adjust(wspace = 0.35)
        if ipreferline: subv = 1 ; subh = 3 ; fig.subplots_adjust(wspace = 0.35)
        #rcParams['font.size'] = int( rcParams['font.size'] * 3. / 4. )
    elif numplot == 4:
        subv = 2
        subh = 2
        #fig.subplots_adjust(wspace = 0.4, hspace = 0.6)
        fig.subplots_adjust(wspace = 0.25, hspace = 0.3)
        #rcParams['font.size'] = int( rcParams['font.size'] * 3. / 4. )
    elif numplot <= 6:
        subv = 2
        subh = 3
        #fig.subplots_adjust(wspace = 0.4, hspace = 0.0)
        fig.subplots_adjust(wspace = 0.5, hspace = 0.3)
        rcParams['font.size'] = int( rcParams['font.size'] * 1. / 2. )
    elif numplot <= 8:
        subv = 2
        subh = 4
        fig.subplots_adjust(wspace = 0.3, hspace = 0.3)
        rcParams['font.size'] = int( rcParams['font.size'] * 1. / 2. )
    elif numplot <= 9:
        subv = 3
        subh = 3
        fig.subplots_adjust(wspace = 0.3, hspace = 0.3)
        rcParams['font.size'] = int( rcParams['font.size'] * 1. / 2. )
    elif numplot <= 12:
        subv = 3
        subh = 4
        fig.subplots_adjust(wspace = 0, hspace = 0.1)
        rcParams['font.size'] = int( rcParams['font.size'] * 1. / 2. )
    elif numplot <= 16:
        subv = 4
        subh = 4
        fig.subplots_adjust(wspace = 0.3, hspace = 0.3)
        rcParams['font.size'] = int( rcParams['font.size'] * 1. / 2. )
    else:
        print "number of plot supported: 1 to 16"
        exit()
    return subv,subh

## Author: AS
def getstralt(nc,nvert):
    varinfile = nc.variables.keys()
    if 'vert' not in varinfile:
        stralt = "_lvl" + str(nvert)
    else:
        zelevel = int(nc.variables['vert'][nvert])
        if abs(zelevel) < 10000.:   strheight=str(zelevel)+"m"
        else:                       strheight=str(int(zelevel/1000.))+"km"
        if 'altitude'       in nc.dimensions:   stralt = "_"+strheight+"-AMR"
        elif 'altitude_abg' in nc.dimensions:   stralt = "_"+strheight+"-ALS"
        elif 'bottom_top'   in nc.dimensions:   stralt = "_"+strheight
        elif 'pressure'     in nc.dimensions:   stralt = "_"+str(zelevel)+"Pa"
        else:                                   stralt = ""
    return stralt

## Author: AS
def getlschar ( namefile, getaxis=False ):
    from netCDF4 import Dataset
    from timestuff import sol2ls
    from numpy import array
    from string import rstrip
    import os as daos
    namefiletest = rstrip( rstrip( rstrip( namefile, chars="_z"), chars="_zabg"), chars="_p")
    testexist = daos.path.isfile(namefiletest)
    zetime = None
    if testexist:  
      namefile = namefiletest
      #### we assume that wrfout is next to wrfout_z and wrfout_zabg
      nc  = Dataset(namefile)
      zetime = None
      days_in_month = [61, 66, 66, 65, 60, 54, 50, 46, 47, 47, 51, 56]
      plus_in_month = [ 0, 61,127,193,258,318,372,422,468,515,562,613]
      if 'Times' in nc.variables:
        zetime = nc.variables['Times'][0]
        shape = array(nc.variables['Times']).shape
        if shape[0] < 2: zetime = None
    if zetime is not None \
       and 'vert' not in nc.variables:
        ##### strangely enough this does not work for api or ncrcat results!
        zesol = plus_in_month[ int(zetime[5]+zetime[6])-1 ] + int(zetime[8]+zetime[9]) - 1 ##les sols GCM commencent a 0
        dals = int( 10. * sol2ls ( zesol ) ) / 10.
        ###
        zetime2 = nc.variables['Times'][1]
        one  = int(zetime[11]+zetime[12]) + int(zetime[14]+zetime[15])/37.
        next = int(zetime2[11]+zetime2[12]) + int(zetime2[14]+zetime2[15])/37. 
        zehour    = one
        zehourin  = abs ( next - one )
        if not getaxis:
            lschar = "_Ls"+str(dals)
        else:
            zelen = len(nc.variables['Times'][:])
            yeye = range(zelen) ; lsaxis = range(zelen) ; solaxis = range(zelen) ; ltaxis = range(zelen)
            for iii in yeye:
               zetime = nc.variables['Times'][iii] 
               ltaxis[iii] = int(zetime[11]+zetime[12]) + int(zetime[14]+zetime[15])/37.
               solaxis[iii] = ltaxis[iii] / 24. + plus_in_month[ int(zetime[5]+zetime[6])-1 ] + int(zetime[8]+zetime[9]) - 1 ##les sols GCM commencent a 0
               lsaxis[iii] = sol2ls ( solaxis[iii] )
               if ltaxis[iii] < ltaxis[iii-1]: ltaxis[iii] = ltaxis[iii] + 24.
               #print ltaxis[iii], solaxis[iii], lsaxis[iii], getattr( nc, 'JULDAY' )
            lschar = lsaxis ; zehour = solaxis ; zehourin = ltaxis
    else:
        lschar=""
        zehour = 0
        zehourin = 1  
    return lschar, zehour, zehourin

## Author: AS
def getprefix (nc):
    prefix = 'LMD_MMM_'
    prefix = prefix + 'd'+str(getattr(nc,'GRID_ID'))+'_'
    prefix = prefix + str(int(getattr(nc,'DX')/1000.))+'km_'
    return prefix

## Author: AS
def getproj (nc):
    typefile = whatkindfile(nc)
    if typefile in ['meso','geo']:
        ### (il faudrait passer CEN_LON dans la projection ?)
        map_proj = getattr(nc, 'MAP_PROJ')
        cen_lat  = getattr(nc, 'CEN_LAT')
        if map_proj == 2:
            if cen_lat > 10.:    
                proj="npstere"
                #print "NP stereographic polar domain" 
            else:            
                proj="spstere"
                #print "SP stereographic polar domain"
        elif map_proj == 1: 
            #print "lambert projection domain" 
            proj="lcc"
        elif map_proj == 3: 
            #print "mercator projection"
            proj="merc"
        else:
            proj="merc"
    elif typefile in ['gcm']:        proj="cyl"    ## pb avec les autres (de trace derriere la sphere ?)
    else:                            proj="ortho"
    return proj    

## Author: AS
def ptitle (name):
    from matplotlib.pyplot import title
    title(name)
    print name

## Author: AS
def polarinterv (lon2d,lat2d):
    import numpy as np
    wlon = [np.min(lon2d),np.max(lon2d)]
    ind = np.array(lat2d).shape[0] / 2  ## to get a good boundlat and to get the pole
    wlat = [np.min(lat2d[ind,:]),np.max(lat2d[ind,:])]
    return [wlon,wlat]

## Author: AS
def simplinterv (lon2d,lat2d):
    import numpy as np
    return [[np.min(lon2d),np.max(lon2d)],[np.min(lat2d),np.max(lat2d)]]

## Author: AS
def wrfinterv (lon2d,lat2d):
    nx = len(lon2d[0,:])-1
    ny = len(lon2d[:,0])-1
    lon1 = lon2d[0,0] 
    lon2 = lon2d[nx,ny] 
    lat1 = lat2d[0,0] 
    lat2 = lat2d[nx,ny] 
    if abs(0.5*(lat1+lat2)) > 60.:  wider = 0.5 * (abs(lon1)+abs(lon2)) * 0.1
    else:                           wider = 0.
    if lon1 < lon2:  wlon = [lon1, lon2 + wider]  
    else:            wlon = [lon2, lon1 + wider]
    if lat1 < lat2:  wlat = [lat1, lat2]
    else:            wlat = [lat2, lat1]
    return [wlon,wlat]

## Author: AS
def makeplotres (filename,res=None,pad_inches_value=0.25,folder='',disp=True,ext='png',erase=False):
    import  matplotlib.pyplot as plt
    from os import system 
    addstr = ""
    if res is not None:
        res = int(res)
        addstr = "_"+str(res)
    name = filename+addstr+"."+ext
    if folder != '':      name = folder+'/'+name
    plt.savefig(name,dpi=res,bbox_inches='tight',pad_inches=pad_inches_value)
    if disp:              display(name)
    if ext in ['eps','ps','svg']:   system("tar czvf "+name+".tar.gz "+name+" ; rm -f "+name)
    if erase:   system("mv "+name+" to_be_erased")		
    return

## Author: AS + AC
def dumpbdy (field,n,stag=None,condition=False,onlyx=False,onlyy=False):
    nx = len(field[0,:])-1
    ny = len(field[:,0])-1
    if condition:
      if stag == 'U': nx = nx-1
      if stag == 'V': ny = ny-1
      if stag == 'W': nx = nx+1 #special les case when we dump stag on W
    if onlyx:    result = field[:,n:nx-n]
    elif onlyy:  result = field[n:ny-n,:]
    else:        result = field[n:ny-n,n:nx-n]
    return result

## Author: AS + AC
def getcoorddef ( nc ):   
    import numpy as np
    ## getcoord2d for predefined types
    typefile = whatkindfile(nc)
    if typefile in ['meso']:
        if '9999' not in getattr(nc,'START_DATE') :   
            ## regular mesoscale 
            [lon2d,lat2d] = getcoord2d(nc)  
        else:                     
            ## idealized mesoscale                      
            nx=getattr(nc,'WEST-EAST_GRID_DIMENSION')
            ny=getattr(nc,'SOUTH-NORTH_GRID_DIMENSION')
            dlat=getattr(nc,'DX')
            ## this is dirty because Martian-specific
            # ... but this just intended to get "lat-lon" like info
            falselon = np.arange(-dlat*(nx-1)/2.,dlat*(nx-1)/2.,dlat)/60000.
            falselat = np.arange(-dlat*(ny-1)/2.,dlat*(ny-1)/2.,dlat)/60000.
            [lon2d,lat2d] = np.meshgrid(falselon,falselat) ## dummy coordinates
            print "WARNING: domain plot artificially centered on lat,lon 0,0"
    elif typefile in ['gcm','earthgcm','ecmwf']: 
        #### n est ce pas nc.variables ?
        if "longitude" in nc.dimensions:  dalon = "longitude"
        elif "lon" in nc.dimensions:      dalon = "lon"
        else:                             dalon = "nothing"
        if "latitude" in nc.dimensions:   dalat = "latitude"
        elif "lat" in nc.dimensions:      dalat = "lat"
        else:                             dalat = "nothing"
        [lon2d,lat2d] = getcoord2d(nc,nlat=dalat,nlon=dalon,is1d=True)
    elif typefile in ['geo']:
        [lon2d,lat2d] = getcoord2d(nc,nlat='XLAT_M',nlon='XLONG_M')
    return lon2d,lat2d    

## Author: AS
def getcoord2d (nc,nlat='XLAT',nlon='XLONG',is1d=False):
    import numpy as np
    if nlon == "nothing" or nlat == "nothing":
        print "NO LAT LON FIELDS. I AM TRYING MY BEST. I ASSUME GLOBAL FIELD."
        lon = np.linspace(-180.,180.,getdimfromvar(nc)[-1])
        lat = np.linspace(-90.,90.,getdimfromvar(nc)[-2])
        [lon2d,lat2d] = np.meshgrid(lon,lat)
    else:
        if is1d:
            lat = nc.variables[nlat][:]
            lon = nc.variables[nlon][:]
            [lon2d,lat2d] = np.meshgrid(lon,lat)
        else:
            lat = nc.variables[nlat][0,:,:]
            lon = nc.variables[nlon][0,:,:]
            [lon2d,lat2d] = [lon,lat]
    return lon2d,lat2d

## Author: AS
def getdimfromvar (nc):
    varinfile = nc.variables.keys()
    dim = nc.variables[varinfile[-1]].shape ## usually the last variable is 4D or 3D
    return dim

## FROM COOKBOOK http://www.scipy.org/Cookbook/SignalSmooth
def smooth1d(x,window_len=11,window='hanning'):
    import numpy
    """smooth the data using a window with requested size.
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
    output:
        the smoothed signal
    example:
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    see also: 
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
    TODO: the window parameter could be the window itself if an array instead of a string   
    """
    x = numpy.array(x)
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')
    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y

## Author: AS
def smooth (field, coeff):
	## actually blur_image could work with different coeff on x and y
	if coeff > 1:	result = blur_image(field,int(coeff))
	else:		result = field
	return result

## FROM COOKBOOK http://www.scipy.org/Cookbook/SignalSmooth
def gauss_kern(size, sizey=None):
	import numpy as np
    	# Returns a normalized 2D gauss kernel array for convolutions
    	size = int(size)
    	if not sizey:
	        sizey = size
	else:
	        sizey = int(sizey)
	x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
	g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
	return g / g.sum()

## FROM COOKBOOK http://www.scipy.org/Cookbook/SignalSmooth
def blur_image(im, n, ny=None) :
	from scipy.signal import convolve
	# blurs the image by convolving with a gaussian kernel of typical size n. 
	# The optional keyword argument ny allows for a different size in the y direction.
    	g = gauss_kern(n, sizey=ny)
    	improc = convolve(im, g, mode='same')
    	return improc

## Author: AS
def getwinddef (nc):    
    ###
    varinfile = nc.variables.keys()
    if 'Um' in varinfile:   [uchar,vchar] = ['Um','Vm'] #; print "this is API meso file"
    elif 'U' in varinfile:  [uchar,vchar] = ['U','V']   #; print "this is RAW meso file"
    elif 'u' in varinfile:  [uchar,vchar] = ['u','v']   #; print "this is GCM file"
    elif 'vitu' in varinfile:  [uchar,vchar] = ['vitu','vitv']   #; print "this is GCM v5 file"
    ### you can add choices here !
    else:                   [uchar,vchar] = ['not found','not found']
    ###
    if uchar in ['U']:         metwind = False ## geometrical (wrt grid) 
    else:                      metwind = True  ## meteorological (zon/mer)
    if metwind is False:       print "Not using meteorological winds. You trust numerical grid as being (x,y)"
    ###
    return uchar,vchar,metwind

## Author: AS
def vectorfield (u, v, x, y, stride=3, scale=15., factor=250., color='black', csmooth=1, key=True):
    ## scale regle la reference du vecteur
    ## factor regle toutes les longueurs (dont la reference). l'AUGMENTER pour raccourcir les vecteurs.
    import  matplotlib.pyplot               as plt
    import  numpy                           as np
    #posx = np.min(x) - np.std(x) / 10.
    #posy = np.min(y) - np.std(y) / 10.
    posx = np.min(x) 
    posy = np.min(y) - 4.*np.std(y) / 10.
    u = smooth(u,csmooth)
    v = smooth(v,csmooth)
    widthvec = 0.003 #0.005 #0.003
    q = plt.quiver( x[::stride,::stride],\
                    y[::stride,::stride],\
                    u[::stride,::stride],\
                    v[::stride,::stride],\
                    angles='xy',color=color,pivot='middle',\
                    scale=factor,width=widthvec )
    if color in ['white','yellow']:     kcolor='black'
    else:                               kcolor=color
    if key: p = plt.quiverkey(q,posx,posy,scale,\
                   str(int(scale)),coordinates='data',color=kcolor,labelpos='S',labelsep = 0.03)
    return 

## Author: AS
def display (name):
    from os import system
    system("display "+name+" > /dev/null 2> /dev/null &")
    return name

## Author: AS
def findstep (wlon):
    steplon = int((wlon[1]-wlon[0])/4.)  #3
    step = 120.
    while step > steplon and step > 15. :       step = step / 2.
    if step <= 15.:
        while step > steplon and step > 5.  :   step = step - 5.
    if step <= 5.:
        while step > steplon and step > 1.  :   step = step - 1.
    if step <= 1.:
        step = 1. 
    return step

## Author: AS
def define_proj (char,wlon,wlat,back=None,blat=None,blon=None):
    from    mpl_toolkits.basemap            import Basemap
    import  numpy                           as np
    import  matplotlib                      as mpl
    from frozen_mymath import max
    meanlon = 0.5*(wlon[0]+wlon[1])
    meanlat = 0.5*(wlat[0]+wlat[1])
    zewidth = np.abs(wlon[0]-wlon[1])*60000.*np.cos(3.14*meanlat/180.)
    zeheight = np.abs(wlat[0]-wlat[1])*60000.
    if blat is None:
        ortholat=meanlat
        if   wlat[0] >= 80.:   blat =  -40. 
        elif wlat[1] <= -80.:  blat = -40.
        elif wlat[1] >= 0.:    blat = wlat[0] 
        elif wlat[0] <= 0.:    blat = wlat[1]
    else:  ortholat=blat
    if blon is None:  ortholon=meanlon
    else:             ortholon=blon
    h = 50.  ## en km
    radius = 3397200.
    if   char == "cyl":     m = Basemap(rsphere=radius,projection='cyl',\
                              llcrnrlat=wlat[0],urcrnrlat=wlat[1],llcrnrlon=wlon[0],urcrnrlon=wlon[1])
    elif char == "moll":    m = Basemap(rsphere=radius,projection='moll',lon_0=meanlon)
    elif char == "ortho":   m = Basemap(rsphere=radius,projection='ortho',lon_0=ortholon,lat_0=ortholat)
    elif char == "lcc":     m = Basemap(rsphere=radius,projection='lcc',lat_1=meanlat,lat_0=meanlat,lon_0=meanlon,\
                              width=zewidth,height=zeheight)
                              #llcrnrlat=wlat[0],urcrnrlat=wlat[1],llcrnrlon=wlon[0],urcrnrlon=wlon[1])
    elif char == "npstere": m = Basemap(rsphere=radius,projection='npstere', boundinglat=blat, lon_0=0.)
    elif char == "spstere": m = Basemap(rsphere=radius,projection='spstere', boundinglat=blat, lon_0=180.)
    elif char == "nplaea":  m = Basemap(rsphere=radius,projection='nplaea', boundinglat=wlat[0], lon_0=meanlon)
    elif char == "laea":    m = Basemap(rsphere=radius,projection='laea',lon_0=meanlon,lat_0=meanlat,lat_ts=meanlat,\
                              width=zewidth,height=zeheight)
                              #llcrnrlat=wlat[0],urcrnrlat=wlat[1],llcrnrlon=wlon[0],urcrnrlon=wlon[1])
    elif char == "nsper":   m = Basemap(rsphere=radius,projection='nsper',lon_0=meanlon,lat_0=meanlat,satellite_height=h*1000.)
    elif char == "merc":    m = Basemap(rsphere=radius,projection='merc',lat_ts=0.,\
                              llcrnrlat=wlat[0],urcrnrlat=wlat[1],llcrnrlon=wlon[0],urcrnrlon=wlon[1])
    elif char == "geos":    m = Basemap(rsphere=radius,projection='geos',lon_0=meanlon)
    elif char == "robin":   m = Basemap(rsphere=radius,projection='robin',lon_0=0)
    elif char == "cass":    
             if zewidth > 60000.:  ## approx. more than one degree
                m = Basemap(rsphere=radius,projection='cass',\
                              lon_0=meanlon,lat_0=meanlat,\
                              width=zewidth,height=zeheight)
             else:
                m = Basemap(rsphere=radius,projection='cass',\
                              lon_0=meanlon,lat_0=meanlat,\
                              llcrnrlat=wlat[0],urcrnrlat=wlat[1],llcrnrlon=wlon[0],urcrnrlon=wlon[1])
    else:                   errormess("projection not supported.")
    fontsizemer = int(mpl.rcParams['font.size']*3./4.)
    if zewidth > 60000.:
        if char in ["cyl","lcc","merc","nsper","laea"]:   step = findstep(wlon)
        else:                                             step = 10.
        steplon = step*2.
    else:
        print "very small domain !"
        steplon = 0.5
        step = 0.5
    zecolor ='grey'
    zelinewidth = 1
    zelatmax = 80.
    if meanlat > 75.:  zelatmax = 90. ; step = step/2. ; steplon = steplon*2.
#    # to show gcm grid:
#    #zecolor = 'r'
#    #zelinewidth = 1
#    #step = 180./48.
#    #steplon = 360./64.
#    #zelatmax = 90. - step/3
#    if char not in ["moll","robin"]:
#        if wlon[1]-wlon[0] < 2.:  ## LOCAL MODE
#            m.drawmeridians(np.r_[-1.:1.:0.05], labels=[0,0,0,1], color=zecolor, linewidth=zelinewidth, fontsize=fontsizemer, fmt='%5.2f')
#            m.drawparallels(np.r_[-1.:1.:0.05], labels=[1,0,0,0], color=zecolor, linewidth=zelinewidth, fontsize=fontsizemer, fmt='%5.2f')
#        else:  ## GLOBAL OR REGIONAL MODE
#            m.drawmeridians(np.r_[-180.:180.:steplon], labels=[0,0,0,1], color=zecolor, linewidth=zelinewidth, fontsize=fontsizemer, latmax=zelatmax)
#            m.drawparallels(np.r_[-90.:90.:step], labels=[1,0,0,0], color=zecolor, linewidth=zelinewidth, fontsize=fontsizemer, latmax=zelatmax)
    if back: 
      if back not in ["coast","sea"]:   m.warpimage(marsmap(back),scale=0.75)
      elif back == "coast":             m.drawcoastlines()
      elif back == "sea":               m.drawlsmask(land_color='white',ocean_color='aqua')
            #if not back:
            #    if not var:                                        back = "mola"    ## if no var:         draw mola
            #    elif typefile in ['mesoapi','meso','geo'] \
            #       and proj not in ['merc','lcc','nsper','laea']:  back = "molabw"  ## if var but meso:   draw molabw
            #    else:                                              pass             ## else:              draw None
    return m

## Author: AS
#### test temporaire
def putpoints (map,plot):
    #### from http://www.scipy.org/Cookbook/Matplotlib/Maps
    # lat/lon coordinates of five cities.
    lats = [18.4]
    lons = [-134.0]
    points=['Olympus Mons']
    # compute the native map projection coordinates for cities.
    x,y = map(lons,lats)
    # plot filled circles at the locations of the cities.
    map.plot(x,y,'bo')
    # plot the names of those five cities.
    wherept = 0 #1000 #50000
    for name,xpt,ypt in zip(points,x,y):
       plot.text(xpt+wherept,ypt+wherept,name)
    ## le nom ne s'affiche pas...
    return

## Author: AS
def calculate_bounds(field,vmin=None,vmax=None):
    import numpy as np
    from frozen_mymath import max,min,mean
    ind = np.where(field < 9e+35)
    fieldcalc = field[ ind ] # la syntaxe compacte ne marche si field est un tuple
    ###
    dev = np.std(fieldcalc)*3.0
    ###
    if vmin is None:  zevmin = mean(fieldcalc) - dev
    else:             zevmin = vmin
    ###
    if vmax is None:  zevmax = mean(fieldcalc) + dev
    else:             zevmax = vmax
    if vmin == vmax:
                      zevmin = mean(fieldcalc) - dev  ### for continuity
                      zevmax = mean(fieldcalc) + dev  ### for continuity            
    ###
    if zevmin < min(fieldcalc): zevmin = min(fieldcalc)
    if zevmax > max(fieldcalc): zevmax = max(fieldcalc)
    #if zevmin < 0. and min(fieldcalc) > 0.: zevmin = 0.
    #print "BOUNDS field ", min(fieldcalc), max(fieldcalc), " //// adopted", zevmin, zevmax
    return zevmin, zevmax

## Author: AS
def bounds(what_I_plot,zevmin,zevmax):
    from frozen_mymath import max,min,mean
    ### might be convenient to add the missing value in arguments
    #what_I_plot[ what_I_plot < zevmin ] = zevmin#*(1. + 1.e-7)
    if zevmin < 0: what_I_plot[ what_I_plot < zevmin*(1. - 1.e-7) ] = zevmin*(1. - 1.e-7)
    else:          what_I_plot[ what_I_plot < zevmin*(1. + 1.e-7) ] = zevmin*(1. + 1.e-7)
    #print "NEW MIN ", min(what_I_plot)
    what_I_plot[ what_I_plot > 9e+35  ] = -9e+35
    what_I_plot[ what_I_plot > zevmax ] = zevmax*(1. - 1.e-7)
    #print "NEW MAX ", max(what_I_plot)
    return what_I_plot

## Author: AS
def nolow(what_I_plot):
    from frozen_mymath import max,min
    lim = 0.15*0.5*(abs(max(what_I_plot))+abs(min(what_I_plot)))
    print "NO PLOT BELOW VALUE ", lim
    what_I_plot [ abs(what_I_plot) < lim ] = 1.e40 
    return what_I_plot


## Author : AC
def hole_bounds(what_I_plot,zevmin,zevmax):
    import numpy as np
    zi=0
    for i in what_I_plot:
        zj=0
        for j in i:
            if ((j < zevmin) or (j > zevmax)):what_I_plot[zi,zj]=np.NaN
            zj=zj+1
        zi=zi+1

    return what_I_plot

## Author: AS
def zoomset (wlon,wlat,zoom):
    dlon = abs(wlon[1]-wlon[0])/2.
    dlat = abs(wlat[1]-wlat[0])/2.
    [wlon,wlat] = [ [wlon[0]+zoom*dlon/100.,wlon[1]-zoom*dlon/100.],\
                    [wlat[0]+zoom*dlat/100.,wlat[1]-zoom*dlat/100.] ]
    print "ZOOM %",zoom,wlon,wlat
    return wlon,wlat

## Author: AS
def fmtvar (whichvar="def"):
    fmtvar    =     { \
             "MIXED":        "%.0f",\
             "UPDRAFT":      "%.0f",\
             "DOWNDRAFT":    "%.0f",\
             "TK":           "%.0f",\
             "T":            "%.0f",\
             "MARS_TI":      "%.0f",\
          "THERMAL_INERTIA": "%.0f",\
             #"ZMAX_TH":      "%.0f",\
             #"WSTAR":        "%.0f",\
             # Variables from TES ncdf format
             "T_NADIR_DAY":  "%.0f",\
             "T_NADIR_NIT":  "%.0f",\
             # Variables from tes.py ncdf format
             "TEMP_DAY":     "%.0f",\
             "TEMP_NIGHT":   "%.0f",\
             # Variables from MCS and mcs.py ncdf format
             "DTEMP":        "%.0f",\
             "NTEMP":        "%.0f",\
             "DNUMBINTEMP":  "%.0f",\
             "NNUMBINTEMP":  "%.0f",\
             # other stuff
             "TPOT":         "%.0f",\
             "TSURF":        "%.0f",\
             "TSK":          "%.0f",\
             "U_OUT1":       "%.0f",\
             "T_OUT1":       "%.0f",\
             "def":          "%.1e",\
             "PTOT":         "%.0f",\
             "PSFC":         "%.1f",\
             "HGT":          "%.1e",\
             "USTM":         "%.2f",\
             "HFX":          "%.0f",\
             "ICETOT":       "%.1e",\
             "TAU_ICE":      "%.2f",\
             "TAUICE":       "%.2f",\
             "VMR_ICE":      "%.1e",\
             "MTOT":         "%.1f",\
             "ANOMALY":      "%.1f",\
             "W":            "%.2f",\
             "WMAX_TH":      "%.1f",\
             "WSTAR":        "%.1f",\
             "QSURFICE":     "%.0f",\
             "UM":           "%.0f",\
             "WIND":         "%.0f",\
             "UVMET":         "%.0f",\
             "UV":         "%.0f",\
             "ALBBARE":      "%.2f",\
             "TAU":          "%.1f",\
             "CO2":          "%.2f",\
             "ENFACT":       "%.1f",\
             "QDUST":        "%.6f",\
             #### T.N.
             "TEMP":         "%.0f",\
             "VMR_H2OICE":   "%.0f",\
             "VMR_H2OVAP":   "%.0f",\
             "TAUTES":       "%.2f",\
             "TAUTESAP":     "%.2f",\

                    }
    if "TSURF" in whichvar: whichvar = "TSURF"
    if whichvar not in fmtvar:
        whichvar = "def"
    return fmtvar[whichvar]

## Author: AS
####################################################################################################################
### Colorbars http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps?action=AttachFile&do=get&target=colormaps3.png
def defcolorb (whichone="def"):
    whichcolorb =    { \
             "def":          "spectral",\
             "HGT":          "spectral",\
             "HGT_M":        "spectral",\
             "TK":           "gist_heat",\
             "TPOT":         "Paired",\
             "TSURF":        "RdBu_r",\
             "TSK":          "RdBu_r",\
             "QH2O":         "PuBu",\
             "USTM":         "YlOrRd",\
             "WIND":         "YlOrRd",\
             #"T_nadir_nit":  "RdBu_r",\
             #"T_nadir_day":  "RdBu_r",\
             "HFX":          "RdYlBu",\
             "ICETOT":       "YlGnBu_r",\
             #"MTOT":         "PuBu",\
             "CCNQ":         "YlOrBr",\
             "CCNN":         "YlOrBr",\
             "TEMP":         "Jet",\
             "TAU_ICE":      "Blues",\
             "TAUICE":       "Blues",\
             "VMR_ICE":      "Blues",\
             "W":            "jet",\
             "WMAX_TH":      "spectral",\
             "ANOMALY":      "RdBu_r",\
             "QSURFICE":     "hot_r",\
             "ALBBARE":      "spectral",\
             "TAU":          "YlOrBr_r",\
             "CO2":          "YlOrBr_r",\
             "MIXED":        "GnBu",\
             #### T.N.
             "MTOT":         "spectral",\
             "H2O_ICE_S":    "RdBu",\
             "VMR_H2OICE":   "PuBu",\
             "VMR_H2OVAP":   "PuBu",\
             "WATERCAPTAG":  "Blues",\
                     }
#W --> spectral ou jet
#spectral BrBG RdBu_r
    #print "predefined colorbars"
    if "TSURF" in whichone: whichone = "TSURF"
    if whichone not in whichcolorb:
        whichone = "def"
    return whichcolorb[whichone]

## Author: AS
def definecolorvec (whichone="def"):
        whichcolor =    { \
                "def":          "black",\
                "vis":          "yellow",\
                "vishires":     "green",\
                "molabw":       "yellow",\
                "mola":         "black",\
                "gist_heat":    "white",\
                "hot":          "tk",\
                "gist_rainbow": "black",\
                "spectral":     "black",\
                "gray":         "red",\
                "PuBu":         "black",\
                "titan":        "red",\
                        }
        if whichone not in whichcolor:
                whichone = "def"
        return whichcolor[whichone]

## Author: AS
def marsmap (whichone="vishires"):
        from os import uname
        mymachine = uname()[1]
        ### not sure about speed-up with this method... looks the same
        if "lmd.jussieu.fr" in mymachine:   domain = "/u/aslmd/WWW/maps/"
        elif "aymeric-laptop" in mymachine: domain = "/home/aymeric/Dropbox/Public/"
        else:                               domain = "http://www.lmd.jussieu.fr/~aslmd/maps/"
	whichlink = 	{ \
		#"vis":		"http://maps.jpl.nasa.gov/pix/mar0kuu2.jpg",\
		#"vishires":	"http://www.lmd.jussieu.fr/~aslmd/maps/MarsMap_2500x1250.jpg",\
                #"geolocal":    "http://dl.dropbox.com/u/11078310/geolocal.jpg",\
		#"mola":	"http://www.lns.cornell.edu/~seb/celestia/mars-mola-2k.jpg",\
		#"molabw":	"http://dl.dropbox.com/u/11078310/MarsElevation_2500x1250.jpg",\
                "thermalday":   domain+"thermalday.jpg",\
                "thermalnight": domain+"thermalnight.jpg",\
                "tesalbedo":    domain+"tesalbedo.jpg",\
                "vis":         domain+"mar0kuu2.jpg",\
                "vishires":    domain+"MarsMap_2500x1250.jpg",\
                "geolocal":    domain+"geolocal.jpg",\
                "mola":        domain+"mars-mola-2k.jpg",\
                "molabw":      domain+"MarsElevation_2500x1250.jpg",\
                "clouds":      "http://www.johnstonsarchive.net/spaceart/marswcloudmap.jpg",\
                "jupiter":     "http://www.mmedia.is/~bjj/data/jupiter_css/jupiter_css.jpg",\
                "jupiter_voy": "http://www.mmedia.is/~bjj/data/jupiter/jupiter_vgr2.jpg",\
                #"bw":          domain+"EarthElevation_2500x1250.jpg",\
                "bw":          "http://users.info.unicaen.fr/~karczma/TEACH/InfoGeo/Images/Planets/EarthElevation_2500x1250.jpg",\
                "contrast":    "http://users.info.unicaen.fr/~karczma/TEACH/InfoGeo/Images/Planets/EarthMapAtmos_2500x1250.jpg",\
                "nice":        "http://users.info.unicaen.fr/~karczma/TEACH/InfoGeo/Images/Planets/earthmap1k.jpg",\
                "blue":        "http://eoimages.gsfc.nasa.gov/ve/2430/land_ocean_ice_2048.jpg",\
                "blueclouds":  "http://eoimages.gsfc.nasa.gov/ve/2431/land_ocean_ice_cloud_2048.jpg",\
                "justclouds":  "http://eoimages.gsfc.nasa.gov/ve/2432/cloud_combined_2048.jpg",\
                "pluto":       "http://www.boulder.swri.edu/~buie/pluto/pluto_all.png",\
                "triton":      "http://laps.noaa.gov/albers/sos/neptune/triton/triton_rgb_cyl_www.jpg",\
                "titan":       "http://laps.noaa.gov/albers/sos/saturn/titan/titan_rgb_cyl_www.jpg",\
                #"titan":       "http://laps.noaa.gov/albers/sos/celestia/titan_50.jpg",\
                "titanuni":    "http://maps.jpl.nasa.gov/pix/sat6fss1.jpg",\
                "venus":       "http://laps.noaa.gov/albers/sos/venus/venus4/venus4_rgb_cyl_www.jpg",\
                "cosmic":      "http://laps.noaa.gov/albers/sos/universe/wmap/wmap_rgb_cyl_www.jpg",\
			}
        ### see http://www.mmedia.is/~bjj/planetary_maps.html
	if whichone not in whichlink: 
		print "marsmap: choice not defined... you'll get the default one... "
		whichone = "vishires"  
        return whichlink[whichone]

#def earthmap (whichone):
#	if   whichone == "contrast":	whichlink="http://users.info.unicaen.fr/~karczma/TEACH/InfoGeo/Images/Planets/EarthMapAtmos_2500x1250.jpg"
#	elif whichone == "bw":		whichlink="http://users.info.unicaen.fr/~karczma/TEACH/InfoGeo/Images/Planets/EarthElevation_2500x1250.jpg"
#	elif whichone == "nice":	whichlink="http://users.info.unicaen.fr/~karczma/TEACH/InfoGeo/Images/Planets/earthmap1k.jpg"
#	return whichlink

## Author: AS
def latinterv (area="Whole"):
    list =    { \
        "Europe":                [[ 20., 80.],[- 50.,  50.]],\
        "Central_America":       [[-10., 40.],[ 230., 300.]],\
        "Africa":                [[-20., 50.],[- 50.,  50.]],\
        "Whole":                 [[-90., 90.],[-180., 180.]],\
        "Southern_Hemisphere":   [[-90., 60.],[-180., 180.]],\
        "Northern_Hemisphere":   [[-60., 90.],[-180., 180.]],\
        "Tharsis":               [[-30., 60.],[-170.,- 10.]],\
        "Whole_No_High":         [[-60., 60.],[-180., 180.]],\
        "Chryse":                [[-60., 60.],[- 60.,  60.]],\
        "North_Pole":            [[ 50., 90.],[-180., 180.]],\
        "Close_North_Pole":      [[ 75., 90.],[-180., 180.]],\
        "Far_South_Pole":        [[-90.,-40.],[-180., 180.]],\
        "South_Pole":            [[-90.,-50.],[-180., 180.]],\
        "Close_South_Pole":      [[-90.,-75.],[-180., 180.]],\
        "Sirenum_Crater_large":  [[-46.,-34.],[-166.,-151.]],\
        "Sirenum_Crater_small":  [[-36.,-26.],[-168.,-156.]],\
        "Rupes":                 [[ 72., 90.],[-120.,- 20.]],\
        "Xanadu":                [[-40., 20.],[  40., 120.]],\
        "Hyperboreae":           [[ 80., 87.],[- 70.,- 10.]],\
              }
    if area not in list:   area = "Whole"
    [olat,olon] = list[area]
    return olon,olat

## Author: TN
def separatenames (name):
  from numpy import concatenate
  # look for comas in the input name to separate different names (files, variables,etc ..)
  if name is None:
     names = None
  else:
    names = []
    stop = 0
    currentname = name
    while stop == 0:
      indexvir = currentname.find(',')
      if indexvir == -1:
        stop = 1
        name1 = currentname
      else:
        name1 = currentname[0:indexvir]
      names = concatenate((names,[name1]))
      currentname = currentname[indexvir+1:len(currentname)]
  return names


## Author: TN
def readslices(saxis):
  from numpy import empty
  if saxis == None:
     zesaxis = None
  else:
     zesaxis = empty((len(saxis),2))
     for i in range(len(saxis)):
        a = separatenames(saxis[i])
        if len(a) == 1:
           zesaxis[i,:] = float(a[0])
        else:
           zesaxis[i,0] = float(a[0])
           zesaxis[i,1] = float(a[1])
           
  return zesaxis

## Author: TN
def readdata(data,datatype,coord1,coord2):
## Read sparse data
  from numpy import empty
  if datatype == 'txt':
     if len(data[coord1].shape) == 1:
         return data[coord1][:]
     elif len(data[coord1].shape) == 2:
         return data[coord1][:,int(coord2)-1]
     else:
         errormess('error in readdata')
  elif datatype == 'sav':
     return data[coord1][coord2]
  else:
     errormess(datatype+' type is not supported!')


## Author: AS
def bidimfind(lon2d,lat2d,vlon,vlat,file=None):
   import numpy as np
   import matplotlib.pyplot as mpl
   if vlat is None:    array = (lon2d - vlon)**2
   elif vlon is None:  array = (lat2d - vlat)**2
   else:               array = (lon2d - vlon)**2 + (lat2d - vlat)**2
   idy,idx = np.unravel_index( np.argmin(array), lon2d.shape )
   if vlon is not None:
      if (np.abs(lon2d[idy,idx]-vlon)) > 5: errormess("longitude not found ",printvar=lon2d)
   if vlat is not None:
      if (np.abs(lat2d[idy,idx]-vlat)) > 5: errormess("latitude not found ",printvar=lat2d)
   if file is not None:
      print idx,idy,lon2d[idy,idx],vlon
      print idx,idy,lat2d[idy,idx],vlat
      var = file.variables["HGT"][:,:,:]
      mpl.contourf(var[0,:,:],30,cmap = mpl.get_cmap(name="Greys_r") ) ; mpl.axis('off') ; mpl.plot(idx,idy,'mx',mew=4.0,ms=20.0)
      mpl.show()
   return idy,idx

## Author: TN
def getsindex(saxis,index,axis):
# input  : all the desired slices and the good index
# output : all indexes to be taken into account for reducing field
  import numpy as np
  if ( np.array(axis).ndim == 2):
      axis = axis[:,0]
  if saxis is None:
      zeindex = None
  else:
      aaa = int(np.argmin(abs(saxis[index,0] - axis)))
      bbb = int(np.argmin(abs(saxis[index,1] - axis)))
      [imin,imax] = np.sort(np.array([aaa,bbb]))
      zeindex = np.array(range(imax-imin+1))+imin
      # because -180 and 180 are the same point in longitude,
      # we get rid of one for averaging purposes.
      if axis[imin] == -180 and axis[imax] == 180:
         zeindex = zeindex[0:len(zeindex)-1]
         print "INFO: whole longitude averaging asked, so last point is not taken into account."
  return zeindex
     
## Author: TN
def define_axis(lon,lat,vert,time,indexlon,indexlat,indexvert,indextime,what_I_plot,dim0,vertmode,redope):
# Purpose of define_axis is to find x and y axis scales in a smart way
# x axis priority: 1/time 2/lon 3/lat 4/vertical
# To be improved !!!...
   from numpy import array,swapaxes
   x = None
   y = None
   count = 0
   what_I_plot = array(what_I_plot)
   shape = what_I_plot.shape
   if indextime is None and len(time) > 1:
      print "AXIS is time"
      x = time
      count = count+1
   if indexlon is None and len(lon) > 1 and redope not in ['edge_x1','edge_x2']:
      print "AXIS is lon"
      if count == 0: x = lon
      else: y = lon
      count = count+1
   if indexlat is None and len(lat) > 1 and redope not in ['edge_y1','edge_y2']:
      print "AXIS is lat"
      if count == 0: x = lat
      else: y = lat
      count = count+1
   if indexvert is None and len(vert) > 1 and ((dim0 == 4) or (y is None)):
      print "AXIS is vert"
      if vertmode == 0: # vertical axis is as is (GCM grid)
         if count == 0: x=range(len(vert))
         else: y=range(len(vert))
         count = count+1
      else: # vertical axis is in kms
         if count == 0: x = vert
         else: y = vert
         count = count+1
   x = array(x)
   y = array(y)
   print "CHECK SHAPE: what_I_plot, x, y", what_I_plot.shape, x.shape, y.shape
   if len(shape) == 1:
       if shape[0] != len(x):           print "WARNING: shape[0] != len(x). Correcting." ; what_I_plot = what_I_plot[0:len(x)]
       if len(y.shape) > 0:             y = ()
   elif len(shape) == 2:
       if shape[1] == len(y) and shape[0] == len(x) and shape[0] != shape[1]:
           print "INFO: swapaxes: ",what_I_plot.shape,shape ; what_I_plot = swapaxes(what_I_plot,0,1)
       else:
           if shape[0] != len(y):       print "WARNING: shape[0] != len(y). Correcting." ; what_I_plot = what_I_plot[0:len(y),:]
           elif shape[1] != len(x):     print "WARNING: shape[1] != len(x). Correcting." ; what_I_plot = what_I_plot[:,0:len(x)]
   elif len(shape) == 3:
       if vertmode < 0: print "not supported. must check array dimensions at some point. not difficult to implement though."
   return what_I_plot,x,y

# Author: TN + AS + AC
def determineplot(slon, slat, svert, stime, redope):
    nlon = 1 # number of longitudinal slices -- 1 is None
    nlat = 1
    nvert = 1
    ntime = 1
    nslices = 1
    if slon is not None:
        length=len(slon[:,0])
        nslices = nslices*length
        nlon = len(slon)
    if slat is not None:
        length=len(slat[:,0])
        nslices = nslices*length
        nlat = len(slat)
    if svert is not None:
        length=len(svert[:,0])
        nslices = nslices*length
        nvert = len(svert)
    if stime is not None:
        length=len(stime[:,0])
        nslices = nslices*length
        ntime = len(stime)
    #else:
    #    nslices = 2  
    mapmode = 0
    if slon is None and slat is None and redope not in ['edge_x1','edge_x2','edge_y1','edge_y2']:
       mapmode = 1 # in this case we plot a map, with the given projection
    return nlon, nlat, nvert, ntime, mapmode, nslices

## Author : AS
def maplatlon( lon,lat,field,\
               proj="cyl",colorb="jet",ndiv=10,zeback="molabw",trans=0.6,title="",\
               vecx=None,vecy=None,stride=2 ):
    ### an easy way to map a field over lat/lon grid
    import numpy as np
    import matplotlib.pyplot as mpl
    from matplotlib.cm import get_cmap
    ## get lon and lat in 2D version. get lat/lon intervals
    numdim = len(np.array(lon).shape)
    if numdim == 2: 	[lon2d,lat2d] = [lon,lat]
    elif numdim == 1:	[lon2d,lat2d] = np.meshgrid(lon,lat)
    else:               errormess("lon and lat arrays must be 1D or 2D")
    #[wlon,wlat] = latinterv()
    [wlon,wlat] = simplinterv(lon2d,lat2d)
    ## define projection and background. define x and y given the projection
    m = define_proj(proj,wlon,wlat,back=zeback,blat=None,blon=None)
    x, y = m(lon2d, lat2d)
    ## define field. bound field.
    what_I_plot = np.transpose(field)
    zevmin, zevmax = calculate_bounds(what_I_plot)  ## vmin=min(what_I_plot_frame), vmax=max(what_I_plot_frame))
    what_I_plot = bounds(what_I_plot,zevmin,zevmax)
    ## define contour field levels. define color palette
    ticks = ndiv + 1
    zelevels = np.linspace(zevmin,zevmax,ticks)
    palette = get_cmap(name=colorb)
    ## contour field
    m.contourf( x, y, what_I_plot, zelevels, cmap = palette, alpha = trans )
    ## draw colorbar
    if proj in ['moll','cyl']:        zeorientation="horizontal" ; zepad = 0.07
    else:                             zeorientation="vertical" ; zepad = 0.03
    #daformat = fmtvar(fvar.upper())
    daformat = "%.0f"
    zecb = mpl.colorbar( fraction=0.05,pad=zepad,format=daformat,orientation=zeorientation,\
                 ticks=np.linspace(zevmin,zevmax,num=min([ticks/2+1,21])),extend='neither',spacing='proportional' ) 
    ## give a title
    if zeorientation == "horizontal": zecb.ax.set_xlabel(title)
    else:                             ptitle(title)
    ## draw vector
    if vecx is not None and vecy is not None:
       [vecx_frame,vecy_frame] = m.rotate_vector( np.transpose(vecx), np.transpose(vecy), lon2d, lat2d ) ## for metwinds
       vectorfield(vecx_frame, vecy_frame, x, y, stride=stride, csmooth=2,\
                                             scale=30., factor=500., color=definecolorvec(colorb), key=True)
    ## scale regle la reference du vecteur. factor regle toutes les longueurs (dont la reference). l'AUGMENTER pour raccourcir les vecteurs.
    return
## Author : AC
## Handles calls to specific computations (e.g. wind norm, enrichment factor...)
def select_getfield(zvarname=None,znc=None,ztypefile=None,mode=None,ztsat=None,ylon=None,ylat=None,yalt=None,ytime=None,analysis=None):
      from frozen_mymath import get_tsat
 
      ## Specific variables are described here:
      # for the mesoscale:
      specificname_meso = ['UV','uv','uvmet','slopexy','SLOPEXY','deltat','DELTAT','hodograph','tk','hodograph_2']
      # for the gcm:
      specificname_gcm = ['enfact']

      ## Check for variable in file:
      if mode == 'check':      
          varname = zvarname
          varinfile=znc.variables.keys()
	  logical_novarname = zvarname not in znc.variables
          logical_nospecificname_meso = not ((ztypefile in ['meso']) and (zvarname in specificname_meso))
          logical_nospecificname_gcm = not ((ztypefile in ['gcm']) and (zvarname in specificname_gcm))
          if ( logical_novarname and logical_nospecificname_meso and logical_nospecificname_gcm ):
              if len(varinfile) == 1:   varname = varinfile[0]
              else:                     varname = False
          ## Return the variable name:
          return varname

      ## Get the corresponding variable:
      if mode == 'getvar':
          plot_x = None ; plot_y = None ;
          ### ----------- 1. saturation temperature
          if zvarname in ["temp","t","T_nadir_nit","T_nadir_day","temp_day","temp_night"] and ztsat:
              tt=getfield(znc,zvarname) ; print "computing Tsat-T, I ASSUME Z-AXIS IS PRESSURE"
              if type(tt).__name__=='MaskedArray':  tt.set_fill_value([np.NaN]) ; tinput=tt.filled()
              else:                                 tinput=tt
              all_var=get_tsat(yalt,tinput,zlon=ylon,zlat=ylat,zalt=yalt,ztime=ytime)
          ### ----------- 2. wind amplitude
          elif ((zvarname in ['UV','uv','uvmet']) and (ztypefile in ['meso']) and (zvarname not in znc.variables)):
              all_var=windamplitude(znc,'amplitude')
          elif ((zvarname in ['hodograph','hodograph_2']) and (ztypefile in ['meso']) and (zvarname not in znc.variables)):
              plot_x, plot_y = windamplitude(znc,zvarname)
              if plot_x is not None: all_var=plot_x # dummy
              else: all_var=plot_y ; plot_x = None ; plot_y = None # Hodograph type 2 is not 'xy' mode
          elif ((zvarname in ['slopexy','SLOPEXY']) and (ztypefile in ['meso']) and (zvarname not in znc.variables)):
              all_var=slopeamplitude(znc)
          ### ------------ 3. Near surface instability
          elif ((zvarname in ['DELTAT','deltat']) and (ztypefile in ['meso']) and (zvarname not in znc.variables)):
              all_var=deltat0t1(znc)
          ### ------------ 4. Enrichment factor
          elif ((ztypefile in ['gcm']) and (zvarname in ['enfact'])):
              all_var=enrichment_factor(znc,ylon,ylat,ytime)
          ### ------------ 5. teta -> temp
          elif ((ztypefile in ['meso']) and (zvarname in ['tk']) and ('tk' not in znc.variables.keys())):
              all_var=teta_to_tk(znc)
          else:
          ### -----------  999. Normal case
              all_var = getfield(znc,zvarname)
          if analysis is not None:
             if analysis in ['histo','density','histodensity']: plot_y=all_var ; plot_x = plot_y
             elif analysis == 'fft': plot_y, plot_x = spectrum(all_var,ytime,yalt,ylat,ylon) ; all_var = plot_y
          return all_var, plot_x, plot_y

# Author : A.C
# FFT is computed before reducefield voluntarily, because we dont want to compute
# ffts on averaged fields (which would kill all waves). Instead, we take the fft everywhere
# (which is not efficient but it is still ok) and then, make the average (if the user wants to)
def spectrum(var,time,vert,lat,lon):
    import numpy as np
    fft=np.fft.fft(var,axis=1)
    N=len(vert)
    step=(vert[1]-vert[0])*1000.
    print "step is: ",step
    fftfreq=np.fft.fftfreq(N,d=step)
    fftfreq=np.fft.fftshift(fftfreq) # spatial FFT => this is the wavenumber
    fft=np.fft.fftshift(fft)
    fftfreq = 1./fftfreq # => wavelength (div by 0 expected, don't panic)
    fft=np.abs(fft) # => amplitude spectrum
#    fft=np.abs(fft)**2 # => power spectrum
    return fft,fftfreq

# Author : A.C.
# Computes temperature from potential temperature for mesoscale files, without the need to use API, i.e. using natural vertical grid
def teta_to_tk(nc):
    import numpy as np
    varinfile = nc.variables.keys() 
    p0=610.
    t0=220.
    r_cp=1./3.89419
    if "T" in varinfile: zteta=getfield(nc,'T')
    else: errormess("you need T in your file.")
    if "PTOT" in varinfile: zptot=getfield(nc,'PTOT')
    else: errormess("you need PTOT in your file.")
    zt=(zteta+220.)*(zptot/p0)**(r_cp)
    return zt

# Author : A.C.
# Find the lon and lat index of the dust devil with the largest pressure gradient
# Steps :
# 1/ convert the chosen PSFC frame to an image of the PSFC anomaly with respect to the mean
# 2/ apply the Sobel operator
# (The Sobel operator performs a 2-D spatial gradient measurement on an image and so emphasizes regions of high spatial frequency that correspond to edges.)
# 3/ find the maximum of the resulting field 
# 4/ find the points in a 5 pixel radius around the maximum for which the value of the Sobel transform is greater than half the maximum
# 5/ define a slab of points encompassing the above selected points, including the potential points 'inside' them (if the above points are a hollow circle for example)
# 6/ in this slab, find the point at which the surface pressure is minimum
def find_devil(nc,indextime):
    import numpy as np
    from scipy import ndimage
    from frozen_mymath import array2image,image2array

    varinfile = nc.variables.keys()
    if "PSFC" not in varinfile: errormess("You need PSFC in your file to find dust devils")
    else: psfc_full=getfield(nc,'PSFC')
    psfc,error=reducefield( psfc_full, d4=indextime)
    psfcim=array2image(1000.*(psfc-psfc.mean()))
    sx = ndimage.sobel(psfcim, axis=0, mode='constant') ; sy = ndimage.sobel(psfcim, axis=1, mode='constant')
    sob = np.hypot(sx, sy)
    zemax=np.max(sob)
    goodvalues = sob[sob >= zemax/2]
    ix = np.in1d(sob.ravel(), goodvalues).reshape(sob.shape)
    idxs,idys=np.where(ix)
    maxvalue = sob[sob == zemax]
    ixmax = np.in1d(sob.ravel(), maxvalue[0]).reshape(sob.shape)
    idxmax,idymax=np.where(ixmax)
    valok=[]
    for i in np.arange(len(idxs)):
        a=np.sqrt((idxmax-idxs[i])**2 + (idymax-idys[i])**2)
        if 0 < a <= 5.*np.sqrt(2.): valok.append(goodvalues[i])
    ix = np.in1d(sob.ravel(), valok).reshape(sob.shape)
    idxs,idys=np.where(ix)
    hyperslab=psfc[np.min(idxs):np.max(idxs),np.min(idys):np.max(idys)]
    idxsub,idysub=np.where(hyperslab==hyperslab.min())
    idx=idxsub[0]+np.min(idxs) ; idy=idysub[0]+np.min(idys)
    return np.int(idx),np.int(idy)
