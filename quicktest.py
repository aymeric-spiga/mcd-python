#! /usr/bin/env python

from mcd import mcd
req = mcd()
req.lat = -4.6
req.lon = 137.4
req.loct = 15
req.xz = 1
req.xdate = 150.6
req.update()
req.printmeanvar()

