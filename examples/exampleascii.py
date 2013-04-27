#! /usr/bin/python

from mcd import mcd


yeah = mcd()

yeah.diurnal()
yeah.getascii(["u","t"],filename="diurnal.txt")


yeah.zonal()
yeah.getascii(["u","t"],filename="zonal.txt")

