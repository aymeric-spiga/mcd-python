#! /usr/bin/env python

from mcd import vcd_class
req = vcd_class()
req.lat = -4.6
req.lon = 137.4
req.loct = 15
req.xz = 1
req.update()
req.printmeanvar()
req.printextvar(94)

