"""
Michael Hirsch
GPLv3+
"""
import logging
from numpy import absolute,asfortranarray,diff,ones,inf,empty_like,isfinite
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from numpy.linalg import norm
from time import time
from warnings import warn
#
from .transcararc import getColumnVER
from .plotsnew import getx0E0


def FitVERopt(L,bn,Phi0,MpDict,sim,cam,Fwd,tInd,makeplot,verbose):
    assert L.ndim==2
    assert bn.ndim==1   and bn.flags['F_CONTIGUOUS']==True
    assert Phi0.ndim==1 and Phi0.flags['F_CONTIGUOUS']==True

    vfit = {}; bfit = {}; Phifit = {'x':None} #in case optim not run
    minverbose=bool(verbose)
#%% scaling brightness
    """
    We could repeatedly downscale simulted brightness in loop, but that consumes a lot of CPU.
    It is equivalent to temporarily upscale observed brightness once before minimization
    Then downscale once after minimization
    """
    bscale = [C.dn2intens for C in cam]
    cInd   = [C.ind       for C in cam]
    bnu = empty_like(bn)
    for s,c in zip(bscale,cInd):
        bnu[c] = bn[c] * s #DONT use 1/intens2dn --that's wrong for real data case!
#%%
    Mp,zTranscar,EK,EKpcolor = MpDict['Mp'],MpDict['ztc'],MpDict['Ek'],MpDict['EKpcolor']

    if sim.useztranscar:
        Tm = Mp
    else: #interpolate A to be on the same altitude grid as b
        warn("using interpolated VER, use caution that peaks aren't missed")
        fint = interp1d(zTranscar, Mp, kind='linear',axis=0) #faster than loop
        Tm = asfortranarray( fint(Fwd['z']) )

    sz,nEnergy = Tm.shape
    assert sz == Fwd['sz']
    assert Tm.flags['F_CONTIGUOUS'] is True

    '''in case optim not run -- don't remove'''
    Phifit['x'] = None #nans((nEnergy,Fwd['sx'])) #don't remove this
    Phifit['EK'] = EK
    Phifit['EKpcolor'] = EKpcolor
#%% optimization
    '''
    Note: Only SLSQP and COBYA allow constraints (Not L-BFGS-B)
    http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html#constrained-minimization-of-multivariate-scalar-functions-minimize
    http://stackoverflow.com/questions/20075714/scipy-minimize-with-constraints
    http://stackoverflow.com/questions/23476152/dynamically-writing-the-objective-function-and-constraints-for-scipy-optimize-mi
    '''

    if 'gaussian' in makeplot or 'optim' in makeplot:
        optimmeth= sim.optimfitmeth.lower()
        maxiter = sim.optimmaxiter #it's already int
        sx = Fwd['sx']

        cons = None
        optimbound = sim.minflux*ones((nEnergy*sx,2))     # lower bound
        optimbound[:,1] = inf  #None seems to give error  # upper bound
        if optimmeth=='nelder-mead':
            optimopt = {'maxiter':maxiter,'disp':minverbose} #100
        elif optimmeth=='bfgs':
            optimopt = {'maxiter':maxiter,'disp':minverbose,'norm':2} #20
        elif optimmeth=='tnc':
            optimopt = {'maxiter':maxiter,'disp':minverbose} #20
        elif optimmeth=='l-bfgs-b':
            # defaults: maxfun=5*nEnergy*sx, maxiter=10
            optimopt = {'maxfun':maxiter*nEnergy*sx,'maxiter':maxiter,
                        'disp':minverbose} #100 maxiter works well
        elif optimmeth=='slsqp':
            optimopt = {'maxiter':maxiter,'disp':minverbose} #2
            cons = {'type': 'ineq',
                    'fun': difffun}
        elif optimmeth=='cobyla':
            optimopt = {'maxiter':maxiter,'disp':minverbose,'rhobeg':1e1,'tol':1} #10
        else:
            raise TypeError('unknown minimization method: {}'.format(optimmeth))

        tic = time()
        #
        Phifit = minimize(optfun,
                        x0=Phi0, #Phi0 is a vector b/c that's what minimize() needs
                        args=(L.tocsr(),Tm,
                              bnu, #scaled version of bn (do once instead of in loop)
                              nEnergy,sx),
                        method=sim.optimfitmeth,
                        bounds=optimbound, #non-negativity
                        constraints=cons,
                        options=optimopt,
                        )
        #
        logging.info('{:0.1f} seconds to fit.'.format(time()-tic))

        logging.info('Minimizer says: {}'.format(Phifit.message))

        Phifit.x = Phifit.x.reshape(nEnergy,sx,order='F')

        logging.info('residual={:.1e} after {} func evaluations.'.format(Phifit.fun,
                                                                  Phifit.nfev))

        # we do this here so that we don't have to carry so many variables around
        vfit['optim'] = getColumnVER(sim.useztranscar,zTranscar, Mp, Phifit.x)
#%% downscale result to complement upscaling
        bfitu = L @ vfit['optim'].ravel(order='F')

        for s,c in zip(bscale,cInd):
            bfitu[c] /= s

        bfit['optim'] = bfitu
#%%
        # this is repeated because the assignment overwrites from minimize()
        Phifit['EK'] = EK
        Phifit['EKpcolor'] = EKpcolor
        # don't remove the two lines above (ek,ekpcolor)
#%% gaussian fit
        #print('max |diff(phi)| = ' + str(np.abs(np.diff(fitp.x, n=1, axis=0)).max()))
        gx0,gE0 = getx0E0(None,Phifit['x'],Phifit['EK'],Fwd['x'],tInd,None,[None])
 #       print(gx0)
#        print(gE0)
        if isfinite([gx0[0],gE0[0]]).all():
            print('Model input: (B_\perp,E_0) = ({:.2f}, {:.0f})'.format(gx0[0],gE0[0]))
        print('Estimated (B_\perp, E_0) = ({:0.2f}, {:0.0f})'.format(gx0[1],gE0[1]))
        Phifit['gx0'] = gx0[1]
        Phifit['gE0'] = gE0[1]

    return vfit,Phifit,Tm,bfit

def optfun(phiinv,L,Tm,b_obs,nEnergy,sx):
    """this provides the quantity to minimize
    Phi0 is a vector b/c that's what minimize needs, reshape is low cost (but this many times?)
    """
    pinv = Tm @ phiinv.reshape(nEnergy,sx,order='F')
    binv = L  @ pinv.ravel(order='F')
    return norm(binv - b_obs, ord=2)

def difffun(jfit,nEnergy=33,sx=109):
    '''used only for slsqp method'''
    # computes difference down columns (top to bottom)
    return 1e5-absolute(diff(jfit.reshape(nEnergy,sx,order='F'), n=1, axis=0)).max()
