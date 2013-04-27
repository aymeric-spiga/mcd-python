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
              
def getmask (field):
        import numpy as np
        if field is None: return None
        if type(field).__name__=='MaskedArray':
              return np.ma.getmask(field)
        else:
              return np.isnan(field)
        

def deg ():
        return u'\u00b0'

def writeascii ( tab, filename ):
    mydata = tab
    myfile = open(filename, 'w')
    for line in mydata:
        zeline = str(line)
        zeline = zeline.replace('[','')
        zeline = zeline.replace(']','')
        myfile.write(zeline + '\n')
    myfile.close()
    return


# A.C. routine to compute saturation temperature
# Be Carefull, when asking for tsat-t, this routine outputs a masked array.
# To be correctly handled, this call to tsat must be done before the call to
# reduce_field, which handles correctly masked array with the new mean() function.
def get_tsat(pressure,temp=None,zlon=None,zlat=None,zalt=None,ztime=None):
    import math as mt
    import numpy as np
    acond=3.2403751E-04
    bcond=7.3383721E-03
    # if temp is not in input, the routine simply outputs the vertical profile
    # of Tsat
    if temp is None:
      # Identify dimensions in temperature field
      output=np.zeros(np.array(pressure).shape)
      if len(np.array(pressure).shape) is 1:
         #pressure field is a 1d column, (i.e. the altitude coordinate)
         #temperature has to have a z-axis
         i=0
         for pp in pressure:
            output[i]=1./(bcond-acond*mt.log(.0095*pp))          
            i=i+1
      else:
         #pressure field is a field present in the file. Unhandled
         #by this routine for now, which only loads unique variables.
         print "3D pressure field not handled for now, exiting in tsat"
         print "Use a vertical pressure coordinate if you want to compute Tsat"
         exit()
    # if temp is in input, the routine computes Tsat-T by detecting where the 
    # vertical axis is in temp
    else:
      output=np.zeros(np.array(temp).shape)
      vardim=get_dim(zlon,zlat,zalt,ztime,temp)
      if 'altitude' not in vardim.keys():
         print 'no altitude coordinate in temperature field for Tsat computation'
         exit()
      zdim=vardim['altitude']
      ndim=len(np.array(temp).shape)
      print '--- in tsat(). vardim,zdim,ndim: ',vardim,zdim,ndim
      i=0
      for pp in pressure:
        if ndim is 1:
           output[i]=1./(bcond-acond*mt.log(.0095*pp))-temp[i]
        elif ndim is 2:
           if zdim is 0:
              output[i,:]=1./(bcond-acond*mt.log(.0095*pp))-temp[i,:]
           elif zdim is 1:
              output[:,i]=1./(bcond-acond*mt.log(.0095*pp))-temp[:,i]
           else:
              print "stop in get_tsat: zdim: ",zdim
              exit()
        elif ndim is 3:
           if zdim is 0:
              output[i,:,:]=1./(bcond-acond*mt.log(.0095*pp))-temp[i,:,:]
           elif zdim is 1:
              output[:,i,:]=1./(bcond-acond*mt.log(.0095*pp))-temp[:,i,:]
           elif zdim is 2:
              output[:,:,i]=1./(bcond-acond*mt.log(.0095*pp))-temp[:,:,i]
           else:
              print "stop in get_tsat: zdim: ",zdim
              exit()
        elif ndim is 4:
           if zdim is 0:
              output[i,:,:,:]=1./(bcond-acond*mt.log(.0095*pp))-temp[i,:,:,:]
           elif zdim is 1:
              output[:,i,:,:]=1./(bcond-acond*mt.log(.0095*pp))-temp[:,i,:,:]
           elif zdim is 2:
              output[:,:,i,:]=1./(bcond-acond*mt.log(.0095*pp))-temp[:,:,i,:]
           elif zdim is 3:
              output[:,:,:,i]=1./(bcond-acond*mt.log(.0095*pp))-temp[:,:,:,i]
           else:
              print "stop in get_tsat: zdim: ", zdim
              exit()
        else:
           print "stop in get_tsat: ndim: ",ndim
           exit()
        i=i+1
    m=np.ma.masked_invalid(temp,copy=False)
    zoutput=np.ma.array(output,mask=m.mask,fill_value=np.NaN)
    return zoutput

# A.C. Dirty routine to determine where are the axis of a variable
def get_dim(zlon,zlat,zalt,ztime,zvar):
   import numpy as np
   nx,ny,nz,nt=0,0,0,0
   if zlon is not None:
      nx=len(zlon)
   if zlat is not None:
      ny=len(zlat)
   if zalt is not None:
      nz=len(zalt)
   if ztime is not None:
      nt=len(ztime)
   zdims={}
   zdims['longitude']=nx
   zdims['latitude']=ny
   zdims['altitude']=nz
   zdims['Time']=nt
   zvardim=np.array(zvar).shape
   ndim=len(zvardim)
   zzvardim=[[]]*ndim
   j=0
   output={}
   for dim in zvardim:
       if dim not in zdims.values():
          print "WARNING -----------------------------"
          print "Dimensions given to subroutine do not match variables dimensions :"
          exit()
       else:
          a=get_key(zdims,dim)
          if len(a) is not 1:
             if j is 0:                ##this should solve most conflicts with Time
                zzvardim[j]=a[1]
             else:
                zzvardim[j]=a[0]
          else:
              zzvardim[j]=a[0]
          output[zzvardim[j]]=j
          j=j+1
   return output

# A.C. routine that gets keys from a dictionnary value
def get_key(self, value):
    """find the key(s) as a list given a value"""
    return [item[0] for item in self.items() if item[1] == value]

# A.C. routine that gets the nearest value index of array and value
def find_nearest(arr,value,axis=None,strict=False):
    import numpy as np
    # Special case when the value is nan
    if value*0 != 0: return np.NaN
    # Check that the value we search is inside the array for the strict mode
    if strict:
       min=arr.min()
       max=arr.max()
       if ((value > max) or (value < min)): return np.NaN

    if type(arr).__name__=='MaskedArray':
       mask=np.ma.getmask(arr)
       idx=np.ma.argmin(np.abs(arr-value),axis=axis)
    # Special case when there are only missing values on the axis
       if mask[idx]:
          idx=np.NaN
    else:
       idx=(np.abs(arr-value)).argmin(axis=axis)
    return idx

# Author: A.C.
def fig2data ( fig ):
    import numpy
    """
    @brief Convert a Matplotlib figure to a 4D numpy array with RGBA channels and return it
    @param fig a matplotlib figure
    @return a numpy 3D array of RGBA values
    """
    # draw the renderer
    fig.canvas.draw ( )
 
    # Get the RGBA buffer from the figure
    w,h = fig.canvas.get_width_height()
    buf = numpy.fromstring ( fig.canvas.tostring_argb(), dtype=numpy.uint8 )
    buf.shape = ( w, h,4 )
 
    # canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
    buf = numpy.roll ( buf, 3, axis = 2 )
    return buf

# Author: A.C.
def fig2img ( fig ):
    import Image
    import numpy
    """
    @brief Convert a Matplotlib figure to a PIL Image in RGBA format and return it
    @param fig a matplotlib figure
    @return a Python Imaging Library ( PIL ) image
    """
    # put the figure pixmap into a numpy array
    buf = fig2data ( fig )
    w, h, d = buf.shape
    return Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )

# Author: A.C.
# Convert a single layer image object (greyscale) to an array
def image2array(im):
    import numpy as np
    if im.mode not in ("L", "F"):
        raise ValueError, ("can only convert single-layer images", im.mode)
    if im.mode == "L":
        a = np.fromstring(im.tostring(), np.uint8)
    else:
        a = np.fromstring(im.tostring(), np.float32)
    a.shape = im.size[1], im.size[0]
    return a

# Author: A.C.
# Convert a 2D array to a single layer image object (greyscale)
def array2image(a):
    import numpy as np
    import Image
    if a.dtype == np.uint8:
        mode = "L"
    elif a.dtype == np.float32:
        mode = "F"
    else:
        raise ValueError, "unsupported image mode"
    return Image.fromstring(mode, (a.shape[1], a.shape[0]), a.tostring())

