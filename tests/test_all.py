#!/usr/bin/env python
"""
Registration case for HiST program
"""
from pathlib import Path
from numpy.testing import assert_allclose
import h5py

#
from histfeas import userinput, hist_figure


def readCheck(Phi0, Phifit, reffn):
    with h5py.File(str(reffn), "r", libver="latest") as f:
        assert_allclose(f["/phifwd/phi"], Phi0[..., 0])
        # noise makes inversion result differ uniquely each run
        xerrpct = (f["/phifwd/x0"] - Phifit[0]["gx0"]) / f["/phifwd/x0"] * 100
        xmsg = "x0 estimation error [km] {:.1f} %".format(xerrpct[0])
        assert abs(xerrpct[0]) < 30, "B_\perp location error out of tolerance"

        Eerrpct = (f["/phifwd/E0"] - Phifit[0]["gE0"]) / f["/phifwd/E0"] * 100
        Emsg = "E0 estimation error [eV] {:.1f} %".format(Eerrpct[0])
        assert abs(Eerrpct[0]) < 30, "E_0 error out of tolerance"

        print(xmsg)
        print(Emsg)


def writeout(regh5):
    with h5py.File(str(regh5), "a", libver="latest") as f:
        f["/phifwd/E0"] = 7500.0
        f["/phifwd/x0"] = 1.0


if __name__ == "__main__":
    #%% path hack
    # allows running from top or test directory
    from os import chdir

    root = Path(__file__).parent
    chdir(str(root))
    #%%
    from histfeas.loadAnalyze import readresults, findxlsh5  # here for matplotlib import

    P = userinput(ini="registration.ini", outdir="out/reg")

    if not P["load"]:
        hist_figure(P)
    #%% find result HDF5
    h5list, P = findxlsh5(P)
    #%% load result
    Phi0, Phifit = readresults(h5list, P)
    #%% check vs known result
    readCheck(Phi0, Phifit, "registration.h5")
    print("\nOK:  simulation registration case")
#%% real data
