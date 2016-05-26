import logging
from matplotlib.pyplot import figure,close,subplots
from matplotlib.ticker import MaxNLocator#,ScalarFormatter# ,LogFormatterMathtext, #for 1e4 -> 1 x 10^4, applied DIRECTLY in format=
from numpy import diff, empty,nan
import h5py
#
from .plotsnew import writeplots,getx0E0,plotB
from .nans import nans


def analyseres(sim,cam,x,xp,Phifwd,Phifit,drn,dhat,vlim,x0true=None,E0true=None,
               makeplot=[None], odir=None):
    #not Phifit tests for None and []
    if Phifwd is None or not Phifit or Phifit[0]['x'] is None or x0true is None or E0true is None:
        return
    """
    we need to fill zeros in jfwd with machine epsilon, since to get avg we need
    to divide by jfwd and we'll get NaN if we divide by zero

    gx0:  [:,0]: true x0 gaussian fit (true),   [:,1]: optim x0 gaussian fit (estimate)
    gE0:  [:,0]: true E0 gaussian fit (true),   [:,1]: optim E0 gaussian fit (estimate)
    """


    #jfwd[jfwd==0] = spacing(1) #a.k.a eps()
    sx = x.size
    nit = len(Phifit) if Phifit is not None else len(Phifwd)

#%% energy flux plot amd calculations
    gx0= nans((nit,2)); gE0 = nans((nit,2))
    Eavgfwdx = nans((nit,sx))
    Eavghatx = nans((nit,sx))

#%% back to work
    for i,jf in enumerate(Phifit):
        #note even if array is F_CONTIGUOUS, argmax is C-order!!
        gx0[i,:],gE0[i,:] = getx0E0(Phifwd[...,i], jf['x'], jf['EK'],x,9999,odir,makeplot)


        print('t={} gaussian 2-D fits for (x,E):\n'
              ' Fwdtrue: {:.2f} {:.0f}\n'
              ' Fwdfit: {:.2f} {:.0f}\n'
              ' Optim: {:.2f} {:.0f}\n'.format(i, x0true[i], E0true[i],
                                                gx0[i,0],  gE0[i,0],
                                                gx0[i,1],  gE0[i,1]))

        Eavgfwdx[i,:],Eavghatx[i,:] = avgcomp(Phifwd[...,i], jf['x'], jf['EK'],x,makeplot,odir)

#%% overall error
    gx0err = gx0[:,1] - x0true #-gx0[:,0]
    gE0err = gE0[:,1] - E0true #-gE0[:,0]
#%% plots
    #extplot(sim,cam,drn,dhat,vlim,makeplot,odir,verbose)


    doplot(x,Phifit,gE0,Eavgfwdx,Eavghatx, makeplot,odir)

    plotgauss(x0true,gx0,gE0,gx0err,gE0err,makeplot,odir)


def doplot(x,Phifit,gE0,Eavgfwdx,Eavghatx, makeplot,odir):
#    with open('cord.csv','r') as e:
#        reader = csv.reader(e, delimiter=',', quoting = csv.QUOTE_NONE);
#        cord = [[r.strip() for r in row] for row in reader][0]

    try:
        if 'optim' in makeplot:
            try:
                fg = figure()
                ax = fg.gca()
                ax.stem([f.fun for f in Phifit])
                ax.set_xlabel('instantiation')
                ax.set_ylabel('$||\hat{b} - b||_2$')
                ax.set_title('Residual $b$')
                ax.xaxis.set_major_locator(MaxNLocator(integer=True))
                writeplots(fg,'error_boptim',9999,makeplot,odir)
            except AttributeError:
                close(fg)
                pass

        if 'fwd' in makeplot and Eavgfwdx:
            fgf = figure()
            ax = fgf.gca()
            ax.semilogy(x,Eavgfwdx.T, marker='.')
            ax.set_xlabel('$B_\perp$ [km]')
            ax.set_ylabel('Expected Value $\overline{E}$ [eV]')
            ax.set_title('Fwd model: Average Energy $\overline{E}$ vs $B_\perp$')
            ax.legend(['{:.0f} eV'.format(g) for g in gE0[:,0]],loc='best',fontsize=9)
            writeplots(fgf,'Eavg_fwd',9999,makeplot,odir)

        if 'optim' in makeplot and Eavghatx:
            fgo = figure()
            ax = fgo.gca()
            ax.semilogy(x,Eavghatx.T, marker='.')
            ax.set_xlabel('$B_\perp$ [km]')
            ax.set_ylabel('Expected Value $\overline{E}$ [eV]')
            ax.set_title('ESTIMATED Average Energy $\overline{\hat{E}}$ vs $B_\perp$')
            ax.legend(['{:.0f} eV'.format(g) for g in gE0[:,0]],loc='best',fontsize=9)
            writeplots(fgo,'Eavg_optim',9999,makeplot,odir)
    except Exception as e:
        logging.info('skipping average energy plotting.   {}'.format(e))

def extplot(sim,cam,drn,dhat,vlim,makeplot,odir):
#%% brightness plot -- plotting ALL at once to show evolution of dispersive event!
    try:
        if 'fwd' in makeplot and drn:
            for i,b in enumerate(drn):
                plotB(b,cam,vlim['b'],9999,19999,makeplot,'$bfwdall',odir)
    # reconstructed brightness plot
        if 'optim' in makeplot and dhat is not None and len(dhat[0])>0:
            for i,b in enumerate(dhat):
                plotB(b,cam,vlim['b'],9999,29999,makeplot,'$bestall', odir)
    except Exception as e:
        logging.info('skipping plotting overall analysis plots of intensity.  {}'.format(e))

def plotgauss(x0true,gx0,gE0,gx0err,gE0err,makeplot,odir):
    fg,(axx,axE) = subplots(1,2,sharey=False)
    axx.stem(x0true,gx0err)
    axx.set_xlabel('$B_\perp$ [km]')
    axx.set_ylabel('$\hat{B}_{\perp,0}$ error [km]')
    axx.set_title('$\hat{B}_{\perp,0}$ error vs. time & position')
    axx.set_ylim(-0.5,0.5)
    axx.set_xlim(-7,7) #[km]

    axE.stem(x0true,gE0err)
    axE.set_xlabel('$B_\perp$ [km]')
    axE.set_ylabel('$\hat{E}_0$ error [km]',labelpad=-0.5)
    axE.set_title('$\hat{E}_0$ error vs. time & position')
    axE.set_ylim(-200,200)
    axE.set_xlim(-7,7) #[km]

    print('B_\perp,0 gaussfit-Estimation-error (fit-true) =' + ' '.join(
                                          ['{:.2f}'.format(h) for h in gx0err]))
    print('E_0 gaussfit-Estimation-error (fit-true) =' + ' '.join(
                                          ['{:.1f}'.format(j) for j in gE0err]))

    if 'h5' in makeplot and odir is not None:
        fout = odir/'fit_results.h5'
        with h5py.File(str(fout),'w',libver='latest') as f:
            f['/tind'] = gx0err.index.values

            f['/gx0/err']=gx0err.values.astype(float)
            f['/gx0/fwdfit']=gx0

            f['/gE0/err']=gE0err.values.astype(float)
            f['/gE0/fwdfit']=gE0

def avgcomp(Phifwd,Phifit,Ek,x,makeplot,odir):
    nEnergy = Phifwd.shape[0]

    dE = empty(nEnergy)
    dE[0] = 9.952 #a priori for this pre-arranged EK
    dE[1:] = diff(Ek) #per Dahlgren matlab code line 276-280
    #from numpy.testing import assert_allclose
    Eavgfwd = ((Phifwd * Ek[:,None] * dE[:,None]).sum(axis=0) /
               (Phifwd * dE[:,None]).sum(axis=0) )

    #xi = 43 #a priori from x=0.5km

    #jfwd1d = jfwd[:,xi]
    #Eavgfwd1d = (jfwd1d * Ek * dE).sum() / (jfwd1d * dE).sum()

    #assert_allclose(Eavgfwd1d,Eavgfwd[xi])

    #print('E_avg: {:0.1f}'.format(Eavgfwd1d))
    logging.info('E_avg: ',Eavgfwd[0]) #TODO how to pick
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

        writeplots(fg,'Eavg_fwd1d',9999,makeplot,odir)

    #%% average energy per x-location
    # formula is per JGR 2013 Dahlgren et al.
    #E_avg = sum(flux*E*dE) / sum(flux*dE)
        Eavgfwdx = ((Phifwd * Ek[:,None] * dE[:,None]).sum(axis=0) /
                          (Phifwd * dE[:,None]).sum(axis=0) )

        Eavghatx =((x * Ek[:,None] * dE[:,None]).sum(axis=0) /
                         (x * dE[:,None]).sum(axis=0)  )
        return Eavgfwdx,Eavghatx
    return nan,nan