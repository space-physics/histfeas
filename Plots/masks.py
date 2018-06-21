#!/usr/bin/env python
"""
We use Bresenham rasterizers for an example showing a non-uniqueness issue with obliquely-viewed aurora
"""
from numpy import zeros
from skimage.morphology import disk
from skimage.draw import circle

rad = 4

print(disk(rad))
