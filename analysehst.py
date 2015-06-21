from __future__ import print_function, division
from plotsnew import getx0E0, plotB
from matplotlib.pyplot import figure
from matplotlib.ticker import MaxNLocator#,ScalarFormatter# ,LogFormatterMathtext, #for 1e4 -> 1 x 10^4, applied DIRECTLY in format=
from numpy import diff, empty
import h5py
from plotsnew import writeplots
from warnings import warn
#
from nans import nans

def analyseres(sim,cam,x,xp,Phifwd,Phifit,drn,dhat,vlim,makeplot,progms,verbose):
    if Phifwd is None or Phifit[0]['x'] is None:
        return
    '''
    we need to fill zeros in jfwd with machine epsilon, since to get avg we need
    to divide by jfwd and we'll get NaN if we divide by zero
    '''
    #jfwd[jfwd==0] = spacing(1) #a.k.a eps()

    sx = x.size
    nEnergy = Phifwd.shape[0]
    nit = len(Phifit)

#    with open('cord.csv','r') as e:
#        reader = csv.reader(e, delimiter=',', quoting = csv.QUOTE_NONE);
#        cord = [[r.strip() for r in row] for row in reader][0]
    afs=14; tfs=16

#%% brightness residual plot
    if 'optim' in makeplot:
        fg = figure()
        ax = fg.gca()
        ax.stem([f.fun for f in Phifit])
        ax.set_xlabel('instantiation',fontsize=afs)
        ax.set_ylabel('$||\hat{b} - b||_2$',fontsize=afs)
        ax.set_title('Residual $b$',fontsize=tfs)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        writeplots(fg,'error_boptim',9999,makeplot,progms)
#%% brightness plot -- plotting ALL at once to show evolution of dispersive event!
    try:
        if 'fwd' in makeplot and drn is not None:
            for i,b in enumerate(drn):
                plotB(b,sim,cam,vlim['b'],9999,makeplot,'$bfwdall',progms,verbose)
    # reconstructed brightness plot
        if 'optim' in makeplot and dhat is not None and len(dhat[0])>0:
            for i,b in enumerate(dhat):
                plotB(b['optim'],sim,cam,vlim['b'],9999,makeplot,'$bestall', progms,verbose)
    except Exception as e:
        warn('skipping plotting overall analysis plots of intensity.  {}'.format(e))
#%% energy flux plot amd calculations
    x0fwd = nans(nit); E0fwd = nans(nit); gx0fwd=nans(nit); gE0fwd=nans(nit)
    x0hat = nans(nit); E0hat = nans(nit); gx0hat=nans(nit); gE0hat=nans(nit)


    Eavgfwdx = nans((nit,sx))
    Eavghatx = nans((nit,sx))
    dE = empty(nEnergy)
    dE[0] = 9.952 #a priori for this pre-arranged EK
    dE[1:] = diff(Phifit[0]['EK']) #per Dahlgren matlab code line 276-280

#%% back to work
    for i,jf in enumerate(Phifit):
        #note even if array is F_CONTIGUOUS, argmax is C-order!!
        gx0fwd[i],gE0fwd[i], x0fwd[i],E0fwd[i] = getx0E0(Phifwd[...,i],jf['EK'],x,9999,progms,makeplot,verbose)
        gx0hat[i],gE0hat[i], x0hat[i],E0hat[i] = getx0E0(jf['x'],     jf['EK'],x,9999,progms,makeplot,verbose)

        print('t={} gaussian 2-D fits for (x,E). Fwd: {:.2f} {:.0f}'
              ' Optim: {:.2f} {:.0f}'.format(i, gx0fwd[i],gE0fwd[i],gx0hat[i],gE0hat[i]))

        trythis(Phifwd[...,i], jf['x'], jf['EK'],x,dE,makeplot,progms,verbose)
#%% average energy per x-location
    # formula is per JGR 2013 Dahlgren et al.
    #E_avg = sum(flux*E*dE) / sum(flux*dE)
        Eavgfwdx[i,:] = ((Phifwd[...,i] * jf['EK'][:,None] * dE[:,None]).sum(axis=0) /
                          (Phifwd[...,i] * dE[:,None]).sum(axis=0) )

        Eavghatx[i,:] =((jf['x'] * jf['EK'][:,None] * dE[:,None]).sum(axis=0) /
                         (jf['x'] * dE[:,None]).sum(axis=0)  )
    if 'fwd' in makeplot:
        fgf = figure()
        ax = fgf.gca()
        ax.semilogy(x,Eavgfwdx.T, marker='.')
        ax.set_xlabel('$B_\perp$ [km]')
        ax.set_ylabel('Expected Value $\overline{E}$ [eV]')
        ax.set_title('Fwd model: Average Energy $\overline{E}$ vs $B_\perp$')
        ax.legend(['{:0.0f} eV'.format(g) for g in E0fwd],loc='best',fontsize=9)
        writeplots(fgf,'Eavg_fwd',9999,makeplot,progms)

    if 'optim' in makeplot:
        fgo = figure()
        ax = fgo.gca()
        ax.semilogy(x,Eavghatx.T, marker='.')
        ax.set_xlabel('$B_\perp$ [km]')
        ax.set_ylabel('Expected Value $\overline{E}$ [eV]')
        ax.set_title('ESTIMATED Average Energy $\overline{\hat{E}}$ vs $B_\perp$')
        ax.legend(['{:0.0f} eV'.format(g) for g in E0fwd],loc='best',fontsize=9)
        writeplots(fgo,'Eavg_optim',9999,makeplot,progms)
#%% overall error
    gx0err = gx0hat-gx0fwd
    gE0err = gE0hat-gE0fwd

    print('B_\perp,0 Estimation-error =' + ' '.join(
                                          ['{:0.2f}'.format(h) for h in gx0err]))
    print('E_0 Estimation-error =' + ' '.join(
                                          ['{:0.1f}'.format(j) for j in gE0err]))

    if 'h5' in makeplot:
        fout = progms + '/fit_results.h5'
        with h5py.File(fout,'w',libver='latest') as f:
            f['/gx0/err']=gx0err
            f['/gx0/fwd']=gx0fwd
            f['/gx0/fit']=gx0hat
            
            f['/gE0/err']=gE0err
            f['/gE0/fwd']=gE0fwd
            f['/gE0/fit']=gE0hat

            f['/x0/err']=x0hat-x0fwd
            f['/x0/fit']=x0hat
            f['/x0/fwd']=x0fwd
            
            f['/E0/err']=E0hat-E0fwd
            f['/E0/fit']=E0hat
            f['/E0/fwd']=E0fwd

def trythis(jfwd,jfit,Ek,x,dE,makeplot,progms,verbose):
    #from numpy.testing import assert_allclose
    Eavgfwd = ((jfwd * Ek[:,None] * dE[:,None]).sum(axis=0) /
               (jfwd * dE[:,None]).sum(axis=0) )

    #xi = 43 #a priori from x=0.5km

    #jfwd1d = jfwd[:,xi]
    #Eavgfwd1d = (jfwd1d * Ek * dE).sum() / (jfwd1d * dE).sum()

    #assert_allclose(Eavgfwd1d,Eavgfwd[xi])

    #print('E_avg: {:0.1f}'.format(Eavgfwd1d))
    if verbose:
        print('E_avg: ',Eavgfwd[0]) #TODO how to pick
    if 'eavg' in makeplot:
        fg = figure()
        ax = fg.gca()
        #ax.loglog(Ek,jfwd1d,marker='.')
        ax.semilogy(x,Eavgfwd)
        ax.set_ylim(bottom=1)
        #ax.set_xlim(50,20e3)
        ax.set_title('2-D cut to 1-D ')
        #ax.set_xlabel('Energy [eV]')
        ax.set_xlabel('x [km]')
        ax.set_ylabel('diff. num. flux')

        writeplots(fg,'Eavg_fwd1d',9999,makeplot,progms)
