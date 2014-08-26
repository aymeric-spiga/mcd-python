#! /usr/bin/env python
from mcd import mcd
yeah = mcd()


yeah.xzs = -5000.
yeah.xze = 15000.

yeah.zkey = 2


yeah.lat = 25.
yeah.lon = 195.
yeah.loct = 4.2
#yeah.loct = 15.
yeah.xdate = 140.






yeah.profile(nd=50)


print yeah.temptab
print yeah.prestab

tpot = yeah.temptab*((610./yeah.prestab)**(1.0/3.9))

print tpot

#yeah.getascii(["t"],filename="profile.txt")

import ppplot as pp
p = pp.plot1d()
p.f = tpot
p.swap = True
p.makeshow()

