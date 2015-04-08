from __future__ import print_function, division
import csv
from plotsnew import getx0E0, plotB
from matplotlib.pyplot import figure,show
from matplotlib.ticker import MaxNLocator,ScalarFormatter# ,LogFormatterMathtext, #for 1e4 -> 1 x 10^4, applied DIRECTLY in format=
from numpy import diff, empty
import h5py
from plotsnew import writeplots

#
from nans import nans

def analyseres(sim,x,xp,usecam,cam,jfwd,jfit,drn,dhat,vlim,makeplot,progms):
    if all([j is None for j in jfwd]):
        return
    '''
    we need to fill zeros in jfwd with machine epsilon, since to get avg we need
    to divide by jfwd and we'll get NaN if we divide by zero
    '''
    #jfwd[jfwd==0] = spacing(1) #a.k.a eps()

    sx = x.size
    nEnergy = jfwd.shape[0]
    nit = len(jfit)

    with open('cord.csv','r') as e:
        reader = csv.reader(e, delimiter=',', quoting = csv.QUOTE_NONE);
        cord = [[r.strip() for r in row] for row in reader][0]
    afs=14; tfs=16

#%% brightness residual plot
    if 'optim' in makeplot:
        fg = figure()
        ax = fg.gca()
        ax.plot([f.fun for f in jfit],
                marker='.',linestyle='none')
        ax.set_xlabel('instantiation',fontsize=afs)
        ax.set_ylabel('$||\hat{b} - b||_2$',fontsize=afs)
        ax.set_title('Residual $b$',fontsize=tfs)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        writeplots(fg,'error_boptim',9999,makeplot,progms)
#%% brightness plot -- plotting ALL at once to show evolution of dispersive event!
    sfmt = ScalarFormatter(useMathText=True) #for 10^3 instead of 1e3
    sfmt.set_powerlimits((-2, 2))
    sfmt.set_scientific(True)
    sfmt.set_useOffset(False)
    try:
        if 'fwd' in makeplot:
            for i,b in enumerate(drn):
                plotB(b,usecam,sim.realdata,cam,vlim['b'],9999,makeplot,'$b_{fwd',
                                                          cord, #pass all cord or it will IndexError if using only some instatiations of phi0
                                                          [sfmt],8727,progms)
    # reconstructed brightness plot
        if 'optim' in makeplot and len(dhat[0])>0:
            for i,b in enumerate(dhat):
                plotB(b['optim'],usecam,sim.realdata,cam,vlim['b'],9999,makeplot,'$b_{optim',
                      cord[nit*i:nit*i+nit],[sfmt],8728,progms)
    except:
        print('** ERROR plotting overall analysis plots of intensity')
#%% energy flux plot amd calculations
    x0fwd = nans(nit); E0fwd = nans(nit)
    x0hat = nans(nit); E0hat = nans(nit)

    Eavgfwdx = nans((nit,sx))
    Eavghatx = nans((nit,sx))
    dE = empty(nEnergy)
    dE[0] = 9.952 #a priori for this pre-arranged EK
    dE[1:] = diff(jfit[0]['EK']) #per Dahlgren matlab code line 276-280

#%% back to work
    for ji,jf in enumerate(jfit):
        #note even if array is F_CONTIGUOUS, argmax is C-order!!
        gx0fwd,gE0fwd, x0fwd[ji],E0fwd[ji] = getx0E0(jfwd[...,ji],jf['EK'],x,9999,progms,makeplot)
        gx0hat,gE0hat, x0hat[ji],E0hat[ji] = getx0E0(jf['x'],     jf['EK'],x,9999,progms,makeplot)

        print('t={} gaussian 2-D fits for (x,E). Fwd: {:0.2f} {:0.1f} Optim: {:0.2f} {:0.1f}'.format(ji, gx0fwd,gE0fwd,gx0hat,gE0hat))

        trythis(jfwd[...,ji], jf['x'], jf['EK'],x,dE,makeplot,progms)
#%% average energy per x-location
    # formula is per JGR 2013 Dahlgren et al.
    #E_avg = sum(flux*E*dE) / sum(flux*dE)
        Eavgfwdx[ji,:] = ((jfwd[...,ji] * jf['EK'][:,None] * dE[:,None]).sum(axis=0) /
                          (jfwd[...,ji] * dE[:,None]).sum(axis=0) )

        Eavghatx[ji,:] =((jf['x'] * jf['EK'][:,None] * dE[:,None]).sum(axis=0) /
                         (jf['x'] * dE[:,None]).sum(axis=0)  )
    if 'fwd' in makeplot:
        fgf = figure()
        ax = fgf.gca()
        ax.semilogy(x,Eavgfwdx.T, marker='.')
        ax.set_xlabel('x [km]')
        ax.set_ylabel('Expected Value[Energy]')
        ax.set_title('Fwd model: Average Energy $\overline{E}$ at each x-location')
        ax.legend(['{:0.0f}'.format(g)+' eV' for g in E0fwd],loc='best',fontsize=9)
        writeplots(fgf,'Eavg_fwd',9999,makeplot,progms)

    if 'optim' in makeplot:
        fgo = figure()
        ax = fgo.gca()
        ax.semilogy(x,Eavghatx.T, marker='.')
        ax.set_xlabel('x [km]')
        ax.set_ylabel('Expected Value[Energy]')
        ax.set_title('ESTIMATED Average Energy $\overline{\hat{E}}$ at each x-location')
        ax.legend(['{:0.0f}'.format(g)+' eV' for g in E0fwd],loc='best',fontsize=9)
        writeplots(fgo,'Eavg_optim',9999,makeplot,progms)
#%% overall error
    x0err = x0hat-x0fwd
    E0err = E0hat-E0fwd

    print('Estimation-error in x-location=' + ' '.join(
                                          ['{:0.2f}'.format(h) for h in x0err]))
    print('Estimation-error in E0 estimate=' + ' '.join(
                                          ['{:0.1f}'.format(j) for j in E0err]))

    if 'h5' in makeplot:
        with h5py.File(progms + 'results.h5',libver='latest') as fid:
            fid.create_dataset('/x0/fit',data=x0hat)
            fid.create_dataset('/x0/fwd',data=x0fwd)
            fid.create_dataset('/E0/fit',data=E0hat)
            fid.create_dataset('/E0/fwd',data=E0fwd)

    if 'show' in makeplot:
        show()

def trythis(jfwd,jfit,Ek,x,dE,makeplot,progms):
    #from numpy.testing import assert_allclose
    Eavgfwd = ((jfwd * Ek[:,None] * dE[:,None]).sum(axis=0) /
               (jfwd * dE[:,None]).sum(axis=0) )

    #xi = 43 #a priori from x=0.5km

    #jfwd1d = jfwd[:,xi]
    #Eavgfwd1d = (jfwd1d * Ek * dE).sum() / (jfwd1d * dE).sum()

    #assert_allclose(Eavgfwd1d,Eavgfwd[xi])

    #print('E_avg: {:0.1f}'.format(Eavgfwd1d))

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
