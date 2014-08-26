#! /usr/bin/env python
from mcd import mcd
import matplotlib.pyplot as mpl
test = mcd()
test.loct = 15.
test.xz = 100.
test.map2d("wind",proj="robin")
mpl.show()



