#! /bin/bash

###
### fmcd_gfortran.sh
###
### ---> This should be used to make MCD fortran stuff directly accessible in python 
### ---> This script is for gfortran, but it is easy to adapt to your other compilers
### ---> A file fmcd.so should be created
### ---> See mcd.py for use in python. Very easy!
###
### AS. 17/04/2012. 
###

### LINKS
NETCDF=/home/marshttp/NETCDF/netcdf64-4.0.1_gfortran_fPIC/
wheremcd=/home/marshttp/MCD_v5.0//mcd/

### LOG FILE
touch fmcd5.log
\rm fmcd5.log

### COPY/PREPARE SOURCES
### perform changes that makes f2py not to fail
sed s/"\!\!'"/"'"/g $wheremcd/call_mcd.F         | sed s/"\!'"/"'"/g | sed s/"\!"/"\n\!"/g > tmp.call_mcd.F
sed s/"\!\!'"/"'"/g $wheremcd/julian.F           | sed s/"\!'"/"'"/g | sed s/"\!"/"\n\!"/g > tmp.julian.F
sed s/"\!\!'"/"'"/g $wheremcd/heights.F          | sed s/"\!'"/"'"/g | sed s/"\!"/"\n\!"/g > tmp.heights.F
sed s/"\!\!'"/"'"/g $wheremcd/constants_mcd.inc  | sed s/"\!'"/"'"/g | sed s/"\!"/"\n\!"/g > constants_mcd.inc

### BUILD THROUGH f2py WHAT IS NECESSARY TO CREATE THE PYTHON FUNCTIONS
touch fmcd5.pyf
\rm fmcd5.pyf
f2py -h fmcd5.pyf -m fmcd5 tmp.call_mcd.F tmp.julian.F tmp.heights.F > fmcd5.log 2>&1

#### IMPORTANT: we teach f2py about variables in the call_mcd subroutines which are intended to be out
sed s/"real :: pres"/"real, intent(out) :: pres"/g fmcd5.pyf | \
sed s/"real :: dens"/"real, intent(out) :: dens"/g | \
sed s/"real :: temp"/"real, intent(out) :: temp"/g | \
sed s/"real :: zonwind"/"real, intent(out) :: zonwind"/g | \
sed s/"real :: merwind"/"real, intent(out) :: merwind"/g | \
sed s/"real dimension(5) :: meanvar"/"real dimension(5),intent(out) :: meanvar"/g | \
sed s/"real dimension(100) :: extvar"/"real dimension(100),intent(out) :: extvar"/g | \
sed s/"real :: seedout"/"real, intent(out) :: seedout"/g | \
sed s/"integer :: ier"/"integer, intent(out) :: ier"/g > fmcd5.pyf.modif
mv fmcd5.pyf.modif fmcd5.pyf

### BUILD
f2py -c fmcd5.pyf tmp.call_mcd.F tmp.julian.F tmp.heights.F --fcompiler=gnu95 \
  -L$NETCDF/lib -lnetcdf \
  -lm -I$NETCDF/include \
  --f90flags="-fPIC" \
  --f77flags="-fPIC" \
  --verbose \
  > fmcd5.log 2>&1
#  --include-paths $NETCDF/include:$NETCDF/lib \ ---> makes it fail

### CLEAN THE PLACE
\rm tmp.call_mcd.F tmp.julian.F tmp.heights.F constants_mcd.inc
