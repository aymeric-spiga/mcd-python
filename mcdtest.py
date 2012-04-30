#! /usr/bin/env python

from mcd import mcd

query = mcd()
query.update()
query.printcoord()
query.printmeanvar()
query.printextvar(22)
query.printallextvar()

query.xz = 30000.
query.printmcd()

query.viking1()
query.xz = 1.
query.xdate = 150.
query.loct = 12.
query.xdate = 90.

import matplotlib.pyplot as mpl

query.diurnal()
query.plot1d(15)
mpl.show()
query.plot1d(["t","p","u","v"])
mpl.show()

query.latlon()
query.map2d("tsurf")
mpl.show()
query.map2d(["ps","mtot","u","olr"])
mpl.show()

import myplot as cuisine

figname = query.getnameset()+'.png'
query.map2d("tsurf")
mpl.savefig(figname,dpi=85,bbox_inches='tight',pad_inches=0.25)
