
## Those are additional functions
## Useful only for plots with mcd.py
## This is not intended to be improved
## -- instead planetoplot will be used one day

def min (field,axis=None): 
        import numpy as np
        if field is None: return None
        if type(field).__name__=='MaskedArray':
              field.set_fill_value(np.NaN)
              return np.ma.array(field).min(axis=axis)
        elif (np.isnan(np.sum(field)) and (type(field).__name__ not in 'MaskedArray')):
              return np.ma.masked_invalid(field).min(axis=axis)
        else: return np.array(field).min(axis=axis)

def max (field,axis=None):
        import numpy as np
        if field is None: return None
        if type(field).__name__=='MaskedArray':
              field.set_fill_value(np.NaN)
              return np.ma.array(field).max(axis=axis)
        elif (np.isnan(np.sum(field)) and (type(field).__name__ not in 'MaskedArray')):
              return np.ma.masked_invalid(field).max(axis=axis) 
        else: return np.array(field).max(axis=axis)

def mean (field,axis=None):
        import numpy as np
        if field is None: return None
        else: 
           if type(field).__name__=='MaskedArray':
              field.set_fill_value(np.NaN)
              zout=np.ma.array(field).mean(axis=axis)
              if axis is not None:
                 zout.set_fill_value(np.NaN)
                 return zout.filled()
              else:return zout
           elif (np.isnan(np.sum(field)) and (type(field).__name__ not in 'MaskedArray')):
              zout=np.ma.masked_invalid(field).mean(axis=axis)
              if axis is not None:
                 zout.set_fill_value([np.NaN])
                 return zout.filled()
              else:return zout
           else: 
              return np.array(field).mean(axis=axis)

def sum (field,axis=None):
        import numpy as np
        if field is None: return None
        else:
           if type(field).__name__=='MaskedArray':
              field.set_fill_value(np.NaN)
              zout=np.ma.array(field).sum(axis=axis)
              if axis is not None:
                 zout.set_fill_value(np.NaN)
                 return zout.filled()
              else:return zout
           elif (np.isnan(np.sum(field)) and (type(field).__name__ not in 'MaskedArray')):
              zout=np.ma.masked_invalid(field).sum(axis=axis)
              if axis is not None:
                 zout.set_fill_value([np.NaN])
                 return zout.filled()
              else:return zout
           else:
              return np.array(field).sum(axis=axis)

#################################################################
#################################################################
#################################################################

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
    if back is not None: 
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
def calculate_bounds(field,vmin=None,vmax=None):
    import numpy as np
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
    lim = 0.15*0.5*(abs(max(what_I_plot))+abs(min(what_I_plot)))
    print "NO PLOT BELOW VALUE ", lim
    what_I_plot [ abs(what_I_plot) < lim ] = 1.e40 
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
     
## Author : AS
def maplatlon( lon,lat,field,\
               proj="cyl",colorb="jet",ndiv=10,zeback=None,trans=0.6,title="",\
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
