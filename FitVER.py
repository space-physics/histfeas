"""
Michael Hirsch
GPLv3+
"""
from __future__ import print_function, division
from numpy import absolute,asfortranarray,diff,zeros,inf,empty_like
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from numpy.linalg import norm
from time import time
from warnings import warn
#
from transcararc import getColumnVER

def FitVERopt(L,bn,Phi0,MpDict,sim,cam,Fwd,makeplot,dbglvl):
    vfit = {}; bfit = {}; jfit = {'x':None} #in case optim not run
#%% scaling brightness
    """
    We could repeatedly downscale brightness in loop, but that consumes a lot of CPU.
    It is equivalent to temporarily upscale observed brightness once before minimization
    Then downscale once after minimization
    """
    bscale = [cam[c].intens2dn for c in cam]
    cInd = [cam[c].ind for c in cam]
    bnu = empty_like(bn)
    for s,c in zip(bscale,cInd):
        bnu[c] = bn[c] / s
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
    jfit['x'] = None #nans((nEnergy,Fwd['sx'])) #don't remove this
    jfit['EK'] = EK
    jfit['EKpcolor'] = EKpcolor
#%% optimization
    '''
    Note: Only SLSQP and COBYA allow constraints (Not L-BFGS-B)
    http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html#constrained-minimization-of-multivariate-scalar-functions-minimize
    http://stackoverflow.com/questions/20075714/scipy-minimize-with-constraints
    http://stackoverflow.com/questions/23476152/dynamically-writing-the-objective-function-and-constraints-for-scipy-optimize-mi
    '''

    if 'gaussian' in makeplot or 'optim' in makeplot:
        optimmeth= sim.optimfitmeth
        maxiter = sim.optimmaxiter #it's already int
        sx = Fwd['sx']

        cons = None
        nonnegbound = zeros((nEnergy*sx,2))
        nonnegbound[:,1] = inf  #None seems to give error
        if optimmeth.lower()=='nelder-mead':
            optimopt = {'maxiter':maxiter,'disp':True} #100
        elif optimmeth.lower()=='bfgs':
            optimopt = {'maxiter':maxiter,'disp':True,'norm':2} #20
        elif optimmeth.lower()=='tnc':
            optimopt = {'maxiter':maxiter,'disp':True} #20
        elif optimmeth.lower()=='l-bfgs-b':
            # defaults: maxfun=5*nEnergy*sx, maxiter=10
            optimopt = {'maxfun':maxiter*nEnergy*sx,'maxiter':maxiter,'disp':True} #100 maxiter works well
        elif optimmeth.lower()=='slsqp':
            optimopt = {'maxiter':maxiter,'disp':True} #2
            cons = {'type': 'ineq',
                    'fun': difffun}
        elif optimmeth.lower()=='cobyla':
            optimopt = {'maxiter':maxiter,'disp':True,'rhobeg':1e1,'tol':1} #10
        else:
            raise TypeError('unknown minimization method: {}'.format(optimmeth))

        tic = time()
        #
        jfit = minimize(optfun,
                        x0=Phi0, #Phi0 is a vector b/c that's what minimize() needs
                        args=(L.tocsr(),Tm,
                              bnu, #scaled version of bn (do once instead of in loop)
                              nEnergy,sx),
                        method=optimmeth,
                        bounds=nonnegbound, #non-negativity
                        constraints=cons,
                        options=optimopt,
                        )
        #
        print('{:0.1f} seconds to fit.'.format(time()-tic))

        jfit.x = jfit.x.reshape(nEnergy,sx,order='F')
        #jfit['optimresidual'] = jfit.fun
        if dbglvl>0:
            print('residual={:0.1f} after {} func evaluations.'.format(jfit.fun,jfit.nfev))

        # we do this here so that we don't have to carry so many variables around
        vfit['optim'] = getColumnVER(sim.useztranscar,zTranscar, Mp, jfit.x, Fwd['z'])
#%% downscale result to complement upscaling
        bfitu = L.dot( vfit['optim'].ravel(order='F') )
        for s,c in zip(bscale,cInd):
            bfitu[c] *= s
        bfit['optim'] = bfitu
#%%
        # this is repeated because the assignment overwrites from minimize()
        jfit['EK'] = EK
        jfit['EKpcolor'] = EKpcolor
        # don't remove the two lines above (ek,ekpcolor)

    return vfit,jfit,Tm,bfit

def optfun(phiinv,L,Tm,b_obs,nEnergy,sx):
    """this provides the quantity to minimize
    Phi0 is a vector b/c that's what minimize needs, reshape is low cost (but this many times?)
    """
    pinv = Tm.dot(phiinv.reshape(nEnergy,sx,order='F'))
    binv = L.dot(pinv.ravel(order='F'))
    return norm(binv - b_obs, ord=2)

def difffun(jfit,nEnergy=33,sx=109):
    '''used only for slsqp method'''
    # computes difference down columns (top to bottom)
    return 1e5-absolute(diff(jfit.reshape(nEnergy,sx,order='F'), n=1, axis=0)).max()
