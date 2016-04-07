**MCD-PYTHON: python-based interface to the Mars Climate Database**

Open source code and contact information [available on github](https://github.com/aymeric-spiga) [no registration needed]

* To get sources through git 
~~~
git clone https://github.com/aymeric-spiga/mcd-python
~~~

* To get sources through SVN 
~~~
svn co https://github.com/aymeric-spiga/mcd-python/trunk mcd-python
~~~

* To get a static ZIP file of the current version of the code, 
[click here](https://github.com/aymeric-spiga/mcd-python/archive/master.zip)

----

**How to install?**

It is assumed you were able to compile successfully the `mcd` sources with `gfortran`. 

The `netCDF` library should be installed on your system. 
Moreover, from our experience, it seems that it must have been build
using the `-fPIC` (for `gfortran`; the name of the option changes with compilers) 
option which generates position independent code suitable for use in a shared library.
An example script is given in the `netcdf` folder.

The installation below relies on `f2py` utility, which is part of the `numpy` package.

 1. Getting the environment variables right: add the `mcd-python` folder to `PYTHONPATH` in your environment file (e.g. `.bashrc`)

        export PYTHONPATH=$PYTHONPATH:adapt_to_your_own/mcd-python

 2. Modify the compile script `compile_fmcd.sh` to link your local `netCDF` libraries and `mcd` distribution (Fortran sources)

 3. Check that `f2py` is included in your `python` library suite.

 4. Run `compile_fmcd.sh` and check for the created `.so` file (its size should be about 1 Mo)

