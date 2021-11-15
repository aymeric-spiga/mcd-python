#! /bin/bash

######################################################################################
###
### compile_gfortran.sh
###
### - requirements: f2py tool + netCDF librairies compiled in Fortran
###
### ---> This should be used to make MCD fortran stuff directly accessible in python 
### ---> This script is for gfortran, but it is easy to adapt to your other compilers
### ---> A file fmcd.so should be created
### ---> See mcd.py for use in python. Very easy!
###
### AS. 17/04/2012. 
###
######################################################################################

NETCDF=../mcd-python/netcdf/gfortran_netcdf-4.0.1/
wheremcd=../MCD_v5.2/
version="5.2"

NETCDF=/home/marshttp/NETCDF/netcdf64-4.0.1_gfortran_fPIC/
wheremcd="/home/marshttp/MCD_5.3/"
version="5.3"

######################################################################################
######################################################################################
######################################################################################

### LOG FILE
num=""
touch fmcd$num.log
\rm fmcd$num.log

### COPY/PREPARE SOURCES
### perform changes that makes f2py not to fail
sed s/"\!\!'"/"'"/g $wheremcd/mcd/call_mcd.F         | sed s/"\!'"/"'"/g | sed -e 's/!/\'$'\n!/g' > tmp.call_mcd.F
sed s/"\!\!'"/"'"/g $wheremcd/mcd/julian.F           | sed s/"\!'"/"'"/g | sed -e 's/!/\'$'\n!/g' > tmp.julian.F
sed s/"\!\!'"/"'"/g $wheremcd/mcd/heights.F          | sed s/"\!'"/"'"/g | sed -e 's/!/\'$'\n!/g' > tmp.heights.F
sed s/"\!\!'"/"'"/g $wheremcd/mcd/constants_mcd.inc  | sed s/"\!'"/"'"/g | sed -e 's/!/\'$'\n!/g' > constants_mcd.inc

### BUILD THROUGH f2py WHAT IS NECESSARY TO CREATE THE PYTHON FUNCTIONS
touch fmcd$num.pyf
\rm fmcd$num.pyf
echo -e "  First f2py call : \n " > fmcd$num.log
f2py -h fmcd$num.pyf -m fmcd$num tmp.call_mcd.F tmp.julian.F tmp.heights.F >> fmcd$num.log 2>&1

### IMPORTANT: we teach f2py about variables in the call_mcd subroutines which are intended to be out
sed s/"real :: pres"/"real, intent(out) :: pres"/g fmcd$num.pyf | \
sed s/"real :: dens"/"real, intent(out) :: dens"/g | \
sed s/"real :: temp"/"real, intent(out) :: temp"/g | \
sed s/"real :: zonwind"/"real, intent(out) :: zonwind"/g | \
sed s/"real :: merwind"/"real, intent(out) :: merwind"/g | \
sed s/"real dimension(5) :: meanvar"/"real dimension(5),intent(out) :: meanvar"/g | \
sed s/"real dimension(100) :: extvar"/"real dimension(100),intent(out) :: extvar"/g | \
sed s/"real :: seedout"/"real, intent(out) :: seedout"/g | \
sed s/"integer :: ier"/"integer, intent(out) :: ier"/g > fmcd$num.pyf.modif
mv fmcd$num.pyf.modif fmcd$num.pyf

### CUSTOMIZE fmcd.pyf TO ADD BUILT-IN INFO ABOUT VERSION and DATA LINKS
datalink=`cd $wheremcd; pwd`"/data/"
echo "    usercode '''" > patchtmp.txt
echo '      char dataloc[] = "'$datalink'";' >> patchtmp.txt
echo '      char dataver[] = "'$version'";' >> patchtmp.txt
echo "    '''" >> patchtmp.txt
sed '/python module fmcd ! in/r patchtmp.txt' fmcd$num.pyf > tmp ; mv tmp fmcd$num.pyf
sed '/interface  ! in :fmcd/r patch.txt' fmcd$num.pyf > tmp ; mv tmp fmcd$num.pyf

### BUILD
echo -e " \n   Second f2py call : \n " >> fmcd$num.log
f2py -c fmcd$num.pyf tmp.call_mcd.F tmp.julian.F tmp.heights.F --fcompiler=gnu95 \
  -L$NETCDF/lib -lnetcdf \
  -lm -I$NETCDF/include \
  --f90flags="-fPIC" \
  --f77flags="-fPIC" \
  --verbose \
  >> fmcd$num.log 2>&1

#### CLEAN THE PLACE
\rm tmp.call_mcd.F tmp.julian.F tmp.heights.F constants_mcd.inc
\rm patchtmp.txt
