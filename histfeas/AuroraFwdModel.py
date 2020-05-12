#!/usr/bin/env python3
"""
generates modeled aurora volume emission rate
Michael Hirsch
GPLv3+
"""
import logging
from numpy import isnan, zeros
from warnings import warn

#
from .arcclass import getver
from .transcararc import getColumnVER


def getSimVER(Phi0, Mp, Fwd, sim, arc, tInd):
    if sim.realdata:
        return

    Pfwd = zeros((Fwd["z"].size, Phi0.shape[1]), float)  # TODO for non transcar cases
    for k, a in arc.items():
        if a.zshape in ("flat", "impulse", "transcar"):
            """
            recall that phi(z,E) is 2-D matrix, and other matrices involved
            must be updated accordingly from original April 2014 try code
            """
            # NOTE: The line below IS **NOT** += !
            return getColumnVER(sim.useztranscar, Mp["ztc"], Mp["Mp"], Phi0)
        else:
            warn("these cases not yet updated for upsampled time!")
            Pfwd += getver(
                Fwd["x"],
                Fwd["z"],
                Mp,
                Phi0,
                a.Wkm,
                a.Hkm,
                a.X0km,
                a.Z0km,
                a.Xshape,
                a.Zshape,
                a.Pnorm,
            )

    if isnan(Pfwd).any():
        logging.critical("NaN detected in VER")
        return None
    else:
        return Pfwd
