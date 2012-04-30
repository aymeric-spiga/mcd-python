#! /bin/bash

###
### fmcd_g95.sh
###
### ---> This should be used to make MCD fortran stuff directly accessible in python 
### ---> This script is for g95, but it is easy to adapt to your other compilers
### ---> A file fmcd.so should be created
### ---> See mcd.py for use in python. Very easy!
###
### AS. 17/04/2012. 
###

### LINKS
NETCDF=/home/aymeric/Science/MODELES/MESOSCALE_DEV/NETCDF/netcdf64-4.0.1_gfortran/
wheremcd=/home/aymeric/Science/MCD_v4.3/mcd/

### LOG FILE
touch fmcd.log
\rm fmcd.log

### COPY/PREPARE SOURCES
### perform changes that makes f2py not to fail
sed s/"\!\!'"/"'"/g $wheremcd/call_mcd.F         | sed s/"\!'"/"'"/g | sed s/"\!"/"\n\!"/g > tmp.call_mcd.F
sed s/"\!\!'"/"'"/g $wheremcd/julian.F           | sed s/"\!'"/"'"/g | sed s/"\!"/"\n\!"/g > tmp.julian.F
sed s/"\!\!'"/"'"/g $wheremcd/heights.F          | sed s/"\!'"/"'"/g | sed s/"\!"/"\n\!"/g > tmp.heights.F
sed s/"\!\!'"/"'"/g $wheremcd/constants_mcd.inc  | sed s/"\!'"/"'"/g | sed s/"\!"/"\n\!"/g > constants_mcd.inc

### BUILD THROUGH f2py WHAT IS NECESSARY TO CREATE THE PYTHON FUNCTIONS
touch fmcd.pyf
\rm fmcd.pyf
f2py -h fmcd.pyf -m fmcd tmp.call_mcd.F tmp.julian.F tmp.heights.F > fmcd.log 2>&1

### IMPORTANT: we teach f2py about variables in the call_mcd subroutines which are intended to be out
sed s/"real :: pres"/"real, intent(out) :: pres"/g fmcd.pyf | \
sed s/"real :: dens"/"real, intent(out) :: dens"/g | \
sed s/"real :: temp"/"real, intent(out) :: temp"/g | \
sed s/"real :: zonwind"/"real, intent(out) :: zonwind"/g | \
sed s/"real :: merwind"/"real, intent(out) :: merwind"/g | \
sed s/"real dimension(5) :: meanvar"/"real dimension(5),intent(out) :: meanvar"/g | \
sed s/"real dimension(100) :: extvar"/"real dimension(100),intent(out) :: extvar"/g | \
sed s/"real :: seedout"/"real, intent(out) :: seedout"/g | \
sed s/"integer :: ier"/"integer, intent(out) :: ier"/g > fmcd.pyf.modif
mv fmcd.pyf.modif fmcd.pyf

### BUILD
f2py -c fmcd.pyf tmp.call_mcd.F tmp.julian.F tmp.heights.F --fcompiler=g95 \
  -L$NETCDF/lib -lnetcdf \
  -lm -I$NETCDF/include \
  --f90flags="-Wall -Wno=112,141,137,155,110 -fno-second-underscore" \
  --verbose \
  > fmcd.log 2>&1
#  --include-paths $NETCDF/include:$NETCDF/lib \ ---> makes it fail

### CLEAN THE PLACE
\rm tmp.call_mcd.F tmp.julian.F tmp.heights.F constants_mcd.inc
