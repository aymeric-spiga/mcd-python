#! /bin/bash

######################################################################################
###
### compile_gfortran.sh
###
### - requirements: f2py tool + netCDF librairies compiled in Fortran
###
### ---> This should be used to make VCD fortran stuff directly accessible in python 
### ---> This script is for gfortran, but it is easy to adapt to your other compilers
### ---> A file fvcd.so should be created
### ---> See mcd.py for use in python. Very easy! 
###
######################################################################################

NETCDF=/home/marshttp/NETCDF/netcdf64-4.0.1_gfortran_fPIC
wherevcd="/home/marshttp/VCD_1.1/"
version="1.1"

######################################################################################
######################################################################################
######################################################################################

### LOG FILE
num=""
touch fvcd$num.log
\rm fvcd$num.log fvcd*.so

### COPY/PREPARE SOURCES
### perform changes that makes f2py not to fail
sed s/"\!\!'"/"'"/g $wherevcd/vcd/VCD_var.F90 | sed s/"\!'"/"'"/g | sed -e 's/!/\'$'\n!/g' > tmp.VCD_var.F90
sed s/"\!\!'"/"'"/g $wherevcd/vcd/VCD.F90     | sed s/"\!'"/"'"/g | sed -e 's/!/\'$'\n!/g' > tmp.VCD.F90
sed s/"\!\!'"/"'"/g $wherevcd/vcd/julian.F90  | sed s/"\!'"/"'"/g | sed -e 's/!/\'$'\n!/g' > tmp.julian.F90

### BUILD THROUGH f2py WHAT IS NECESSARY TO CREATE THE PYTHON FUNCTIONS
touch fvcd$num.pyf
\rm fvcd$num.pyf
echo -e "  First f2py call : \n " > fvcd$num.log
f2py -h fvcd$num.pyf -m fvcd$num tmp.VCD_var.F90 tmp.VCD.F90 tmp.julian.F90 >> fvcd$num.log 2>&1

### IMPORTANT: we teach f2py about variables in the call_mcd subroutines which are intended to be out
sed s/"real :: pres"/"real, intent(out) :: pres"/g fvcd$num.pyf | \
sed s/"real :: dens"/"real, intent(out) :: dens"/g | \
sed s/"real :: temp"/"real, intent(out) :: temp"/g | \
sed s/"real :: zonwind"/"real, intent(out) :: zonwind"/g | \
sed s/"real :: merwind"/"real, intent(out) :: merwind"/g | \
sed s/"real dimension(5) :: meanvar"/"real dimension(5),intent(out) :: meanvar"/g | \
sed s/"real dimension(100) :: extvar"/"real dimension(100),intent(out) :: extvar"/g | \
sed s/"real :: seedout"/"real, intent(out) :: seedout"/g | \
sed s/"integer :: ier"/"integer, intent(out) :: ier"/g > fvcd$num.pyf.modif
mv fvcd$num.pyf.modif fvcd$num.pyf

### CUSTOMIZE fvcd.pyf TO ADD BUILT-IN INFO ABOUT VERSION and DATA LINKS
datalink=`cd $wherevcd; pwd`"/data/"
echo "    usercode '''" > patchtmp.txt
echo '      char dataloc[] = "'$datalink'";' >> patchtmp.txt
echo '      char dataver[] = "'$version'";' >> patchtmp.txt
echo "    '''" >> patchtmp.txt
sed '/python module fvcd ! in/r patchtmp.txt' fvcd$num.pyf > tmp ; mv tmp fvcd$num.pyf
sed '/interface  ! in :fvcd/r patch.txt' fvcd$num.pyf > tmp ; mv tmp fvcd$num.pyf

### BUILD
echo -e " \n   Second f2py call : \n " >> fvcd$num.log
f2py -c fvcd$num.pyf tmp.VCD_var.F90 tmp.VCD.F90 tmp.julian.F90 --fcompiler=gnu95 \
  -L$NETCDF/lib -lnetcdf \
  -lm -I$NETCDF/include \
  --f90flags="-fPIC -ffree-form -ffree-line-length-none" \
  --f77flags="-fPIC" \
  --verbose \
  >> fvcd$num.log 2>&1

#### CLEAN THE PLACE
\rm tmp.VCD_var.F90 tmp.VCD.F90 tmp.julian.F90
\rm patchtmp.txt
