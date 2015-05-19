from __future__ import print_function, division,absolute_import
from numpy import (in1d,s_,empty,empty_like,isnan,asfortranarray,linspace,outer,
                   sin,cos,pi,ones_like,array,nan,rint,unravel_index,meshgrid,
                   )
from matplotlib.pyplot import (figure,subplots, clf,text,draw)
#from matplotlib.cm import get_cmap
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext, MultipleLocator, ScalarFormatter #for 1e4 -> 1 x 10^4, applied DIRECTLY in format=
#from matplotlib.image import imsave
import h5py
from os.path import join
from scipy.interpolate import interp1d
#
#import seaborn as sns
#sns.set_style('whitegrid')
#sns.set_palette('husl')
#
try:
    import plotly.plotly as py
    from plotly.graph_objs import Data,Figure,XAxis,YAxis,Contour, Layout
except:
    pass
#
from gaussfitter import gaussfit,twodgaussian
from histutils.findnearest import find_nearest
#%% plot globals
afs = 20
tfs = 22
tkfs = 16
longtitle=False
#pcmcmap = get_cmap('jet')
#pcmcmap.set_under('white')
dymaj=100
dymin=20
format1d='png'
plotdpi=75
pstyle='contour'

#%%
def logfmt(makeplot,powlim=(-2,2)):
    """
    call this in EACH function--don't try to reuse or one plot will affect
    the other, particularly a problem with colorbar() (in matplotlib 1.4.2/1.4.3 at least)
    """
    cnorm = [None]
    sfmt = ScalarFormatter(useMathText=True) #for 10^3 instead of 1e3
    sfmt.set_powerlimits(powlim) #force scientific notation for numbers with 10^a where A<a<B
    sfmt.set_scientific(True)
    sfmt.set_useOffset(False)
    sfmt = [sfmt]
    if 'log' in makeplot: #it will plot log and linear!
        cnorm.append(LogNorm())
        sfmt.append(LogFormatterMathtext())
    return cnorm, sfmt

def placetxt(x,y,txt,ax):
    ax.text(x, y, txt,
            color='white',fontsize=45,
            va='bottom',ha='left',
            bbox=dict(boxstyle="round,pad=0.0",fc='black', alpha=0.25))

def goPlot(ParamFN,sim,arc,Fwd,cam,L,Tm,drn,dhat,ver,vfit,Phi0,
            fitp,rawdata,tInd,makeplot,progms,x1d,vlim,dbglvl):

    spfid = ''.join((progms,'dump_t{}.h5'.format(tInd)))
    cord = ['b','g','k','r','m','y']

    xKM = Fwd['x']
    zKM = Fwd['z']
    xp = Fwd['xPixCorn']
    zp = Fwd['zPixCorn']
    sx = Fwd['sx']
    sz = Fwd['sz']
    nCutPix = sim.nCutPix #FIXME assumes all cams same # of pixels

#%% get xind
    if not sim.realdata:
        if x1d is not None:
            Jxi = find_nearest(xKM,x1d)[0]
        else:
            Jxi = find_nearest(xKM,arc.x0)[0]
    else:
        Jxi = None

#%% eigenfunction
    if 'eig' in makeplot:
        ploteig(fitp['EKpcolor'],zKM,Tm,vlim['p'],sim,tInd,makeplot,'peig',progms)

    if 'eig1d' in makeplot:
        ploteig1d(fitp['EK'],zKM,Tm,vlim['p'],sim,tInd,makeplot,'peig1d',progms)

    if 'tphi0' in makeplot:
        plottphi0(Tm,Phi0,Jxi,fitp['EK'],zKM,vlim['p'],sim,tInd,makeplot,'tphi0',progms)
#%% optional show plots
    if 'realvid' in makeplot and sim.realdata:
        plotRealImg(sim,cam,rawdata,tInd,makeplot,'realdata','$I$',1830,
                               spfid,'gray',progms)

    if 'singleraw' in makeplot and sim.realdata:
        plotPlainImg(sim,cam,rawdata,tInd,makeplot,'realdata','$I$',1831,
                               spfid,'gray',progms)

#%% scatter plot of LOS
    if 'kml' in makeplot or 'kmlrays' in makeplot:
        planviewkml(cam,xKM,zKM,makeplot,5289,progms)
    if 'ray3' in makeplot:
        planview3(cam,xKM,zKM,makeplot,6289,progms)
    if 'sitemap' in makeplot:
        from sitemap import fancymap
        fancymap(cam)
#%%
    if 'ell' in makeplot:
        plotEachRay = False
        from EllLineLength import plotEll
        xzplot=None
        xFOVpixelEnds = empty((nCutPix, sim.nCamUsed),dtype=float)
        zFOVpixelEnds = empty_like(xFOVpixelEnds)
        xCam = empty(sim.nCamUsed,dtype=float)
        zCam = empty_like(xCam)
        for i,c in enumerate(cam):
            xFOVpixelEnds[:,i] = cam[c].xFOVpixelEnds
            zFOVpixelEnds[:,i] = cam[c].zFOVpixelEnds
            xCam[i] = cam[c].x_km
            zCam[i] = cam[c].z_km
        #we kept plotEll in EllLineLength.py for plotEachRay case :(
        plotEll(sim.nCamUsed,xFOVpixelEnds,zFOVpixelEnds,xCam,zCam,nCutPix,
                xp,zp,sz,sx, xzplot,sim.FwdLfn,plotEachRay, makeplot,vlim['p'])
#%% Picard plots
    if 'picardL' in makeplot:
        plotPicard(L,drn,'L')
    if 'picardLT' in makeplot:
        LxColInd = s_[Jxi*sz:(Jxi+1)*sz] #faster than fancy indexing
        Lx = L[:,LxColInd] #ellipses here doesn't work
        LT = Lx.dot(Tm)
        plotPicard(LT,drn,'LT')
#%% 1-D slice plots

    if 'bartrecon' in makeplot:
#        plotB(drn,cam,tInd,'$b$',cord,
#              axm,car,cac,'oneperplot')
#        plotB(dhat['art'],cam,tInd,'$\hat{b}$',cord[::-1],
#              axm,car,cac,['twinx','oneperplot'])
        plotBcompare(drn,dhat['artrecon'],cam,sim.nCamUsed,
                     nCutPix,'ART',cord,vlim['b'],tInd,progms)
#%% error plots
    if 'berror' in makeplot:
        plotB(drn - dhat['art'],sim.realdata,cam,nCutPix,vlim['b'],
                            tInd,makeplot,'$\Delta{b}$',cord,progms)
#%% characteristic energy determination for title labels
    #set_trace()
    gx0,gE0,x0,E0 = getx0E0(Phi0,fitp['EK'],xKM,tInd,progms,makeplot)
    if longtitle:
        fwdloctxt = ('\n$x_{{cam}}=${} $(x_0,E_0)=({:0.2f},{:0.0f})$ [km,eV]'.format(sim.allCamXkm,x0,E0))
    else:
        fwdloctxt=''

    bcomptxt = ('best')
               #'($x_0$,$E_0$)=(' +'{:0.2f}'.format(x0) + ',' + '{:0.0f}'.format(E0) + ') [km,eV]')
#%% ART Energy plots
    if 'jartadj' in makeplot:
        #fa,aa = subplots(nrows=1,ncols=2,sharex='col',num=992, figsize=(5,14))
        plotJ(sim,fitp['phiARTadjoint'],xKM,xp,fitp['EK'],fitp['EKpcolor'],vlim['j'],tInd,makeplot,'jartadj',
              '$\phi_{art,adj}[E,x]$ Estimated (ADJOINT) from ART VER',
              spfid,progms)
#%% Forward model plots
    if 'fwd' in makeplot or 'optim' in makeplot:
        plotnoise(cam,tInd,makeplot,'bnoise',progms)

    if 'fwd' in makeplot:
        plotfwd(sim,cam,drn,nCutPix,xKM,xp,zKM,zp,
                ver,Phi0,fitp,Jxi,vlim,tInd,makeplot,cord,fwdloctxt,spfid,progms)
#%% gaussian fit of optim
    if 'gaussian' in makeplot and 'fwd' in makeplot:
        plotBcompare(drn,dhat['gaussian'],cam,sim.nCamUsed,nCutPix,
                     bcomptxt, 'gaussfit',cord,vlim['b'],tInd, makeplot,progms)

        gx0hat,gE0hat,x0hat,E0hat = getx0E0(fitp['gaussian'],fitp['EK'],xKM,tInd,progms,makeplot)
#'Neval = {:d}'.format(fitp.nfev)
        plotJ(sim,fitp['gaussian'], xKM,xp, fitp['EK'],fitp['EKpcolor'], vlim['j'],tInd, makeplot,'jgaussian',
              ('$\widehat{\phi_{gaussian,optim}}[E,x]$ estimated diff. number flux' +
               fwdloctxt ),
              spfid,progms)

        print('Estimated $x_{{gauss,0}},E_{{gauss,0}}$={:0.2f}'.format(gx0hat) +
              ',' + '{:0.0f}'.format(gE0hat))

        plotVER(sim,vfit['gaussian'],xKM,xp,zKM,zp,vlim['p'],tInd,makeplot,'vgaussian',
              ('$\widehat{P_{gaussian,optim}}$ volume emission rate' + fwdloctxt),
              spfid,progms)
#%% optimize search plots
    if 'optim' in makeplot:
        plotoptim(sim,fitp,xKM,xp,zKM,zp,vlim,tInd,makeplot,spfid,progms)
#%% Adjoint TJ plots
    if 'adj' in makeplot:
        plotadj(sim,cam,drn,dhat,nCutPix,xKM,xp,zKM,zp,vfit,fitp,vlim,
                                 bcomptxt,tInd,makeplot,cord,spfid,progms)
#%% maximum entropy
    if 'phimaxent' in makeplot:
        plotJ(sim,fitp['maxent'],xKM,xp,fitp['EK'],fitp['EKpcolor'],vlim['j'],tInd,makeplot,'jme',
                              '$\Phi_{maxent}$ estimated diff. number flux',
                            spfid,progms)

    if 'phimaxent1d' in makeplot and Jxi is not None:
        plotJ1D(sim,fitp['maxent'][:,Jxi],fitp['EK'],vlim['j'],tInd,makeplot,'jme_1D',
                ('Differential Number flux at $B_\perp$={:0.2f} [km]'.format(xKM[Jxi])),
                spfid,progms)
    if 'bmaxent' in makeplot:
        plotBcompare(drn,dhat['fit_maxent'],cam,sim.nCamUsed, nCutPix,
                     bcomptxt, 'maxent',cord,vlim['b'],tInd,makeplot,progms)
#        plotB(drn,cam,tInd,'$b$',cord,
#              axm,car,cac,'oneperplot')
#        plotB(dhat['fit_maxent'],cam,tInd,'$b_{fit_maxent}$',cord[::-1],
#              axm,car,cac,['twinx','oneperplot'])
    if 'pmaxent' in makeplot:
        plotVER(sim,vfit['maxent'],xKM,xp,zKM,zp,vlim['p'],tInd,makeplot,'maxent',
              '$\hat{v}_{maxent}$ from maximum entropy regularization',
              spfid,progms)
    if 'pmaxent1d' in makeplot and Jxi is not None:
        plotVER1D(sim,vfit['maxent'][:,Jxi],zKM,vlim['p'][2:],tInd,makeplot,
                  'vermaxent_1D',
                  ('$p_{optim,maxent}$  $B_\perp$={:0.2f}'.format(xKM[Jxi]) +' [km] ' + fwdloctxt),
                  13994,spfid,progms)
        if 'pfwd1d' in makeplot:
            try:
                cax = figure().gca()
                cax.plot(vfit['maxent'][:,Jxi], zKM,
                         label='pmaxent',color='red')
                cax.plot(ver[:,Jxi],zKM,
                         label='pfwd', color='blue',linestyle='--')
                cax.legend(loc='best')
                cax.set_xlabel('VER',labelpad=0)
                cax.set_ylabel('z [km]',labelpad=0)
                cax.set_title('VER forward model vs. VER maximum entropy reconstruction')
                cax.grid(True)
                cax.yaxis.set_major_locator(MultipleLocator(dymaj))
                cax.yaxis.set_minor_locator(MultipleLocator(dymin))
                print('max diff vmaxent-vfwd=' + str((vfit['maxent'][:,Jxi]-ver[:,Jxi]).max()))
            except Exception as e:
                print('* could not plot vfwd vmaxent comparison.  {}'.format(e))
#%% diff number flux from ART
    if 'jart' in makeplot:
        plotJ(sim,fitp['art'],xKM,xp,fitp['EK'],fitp['EKpcolor'],vlim['j'],tInd,makeplot,'jart',
                        '$J_{art}[E,x]$ Estimated J from Kaczmarz ART on LT and b',
                         spfid,progms)
    if 'vart' in makeplot:
        assert isnan(vfit['art']).any() == False
        plotVER(sim,vfit['art'],xKM,xp,zKM,zp,vlim['p'],tInd,makeplot,'vart',
              '$\hat{P}_{art}$ from ART estimation of $J$',
              spfid,progms)
    if 'bart' in makeplot:
        plotBcompare(drn,dhat['fit_art'],cam,sim.nCamUsed,nCutPix,
                     bcomptxt,'art',cord,vlim['b'],tInd, makeplot,progms)

#%%
def plotfwd(sim,cam,drn,nCutPix,xKM,xp,zKM,zp,
            ver,Phi0,fitp,Jxi,vlim,tInd,makeplot,cord,fwdloctxt,spfid,progms):

    plotB(drn,sim.realdata,cam,nCutPix,vlim['b'],tInd,makeplot,'$br',cord,progms)

    if not sim.realdata:
        # Forward model VER
        plotVER(sim,ver,xKM,xp,zKM,zp,vlim['p'],tInd,makeplot,'pfwd',
            ('$P_{{fwd}}$ volume emission rate' + fwdloctxt),
                    spfid,progms)

#       print('max |diff(phifwd)| = ' + str(np.abs(np.diff(phiInit, n=1, axis=0)).max()))
        plotJ(sim,Phi0,xKM,xp,fitp['EK'],fitp['EKpcolor'],vlim['j'],tInd,makeplot,'phifwd',
              ('$\Phi_{{top}}$ input diff. number flux'+ fwdloctxt),
                        spfid,progms)

        if not 'optim' in makeplot:
            plotJ1D(sim,Phi0[:,Jxi],None,fitp['EK'],vlim['j'],tInd,makeplot,'phifwd1d',
                 ('Differential Number flux at $B_\perp$={:0.2f} [km]'.format(xKM[Jxi])),
                       spfid,progms)


        if not 'optim' in makeplot:
            plotVER1D(sim,ver[:,Jxi],None,zKM,vlim['p'][2:],tInd,makeplot,'pfwd1d',
              ('$P_{{fwd}}$ at $B_\perp$={:0.2f}'.format(xKM[Jxi]) +
               ' [km]' + fwdloctxt),
              119900,spfid,progms)
#%%
def plotoptim(sim,cam,drn,dhat,nCutPix,bcomptxt,fwdloctxt,ver,Phi0,Jxi,
              vfit,fitp,xKM,xp,zKM,zp,vlim,tInd,cord,makeplot,spfid,progms):
    print(fitp.message)

    plotBcompare(drn,dhat['optim'],cam,sim.nCamUsed,nCutPix,
                 bcomptxt, 'est',cord,vlim['b'],tInd,makeplot,progms)

    plotVER(sim,vfit['optim'],xKM,xp,zKM,zp,vlim['p'],tInd,makeplot,'pest',
          ('$\widehat{P}$ volume emission rate' + fwdloctxt),
          spfid,progms)

    #print('max |diff(phi)| = ' + str(np.abs(np.diff(fitp.x, n=1, axis=0)).max()))
    gx0hat,gE0hat,x0hat,E0hat = getx0E0(fitp.x,fitp['EK'],xKM,tInd,progms,makeplot)
    print('Estimated $x_{{gauss,0}},E_{{gauss,0}}$={:0.2f}'.format(gx0hat) +
          ',' + '{:0.0f}'.format(gE0hat))

    plotJ(sim,fitp.x,xKM,xp,fitp['EK'],fitp['EKpcolor'],vlim['j'],tInd,makeplot,'phiest',
          ('$\widehat{\phi}[E,x]$ estimated diff. number flux' + fwdloctxt ),
          spfid,progms)
          #'Neval = {:d}'.format(fitp.nfev)

    if Jxi is not None:
        plotVER1D(sim,ver[:,Jxi],vfit['optim'][:,Jxi],zKM,vlim['p'][2:],tInd,makeplot,'pest1d',
                  ('$\widehat{{P}}$ at $B_\perp$={:0.2f} [km]'.format(xKM[Jxi])
                   + fwdloctxt),
                  119900,spfid,progms)

        plotJ1D(sim,Phi0[:,Jxi],fitp.x[:,Jxi],fitp['EK'],vlim['j'],tInd,makeplot,'phiest1d',
                   ('Differential Number flux at $B_\perp$={:0.2f} [km]'.format(xKM[Jxi])),
                           spfid,progms)
#%%
def plotadj(sim,cam,drn,dhat,nCutPix,xKM,xp,zKM,zp,vfit,fitp,vlim,
                                 bcomptxt,tInd,makeplot,cord,spfid,progms):
    if 'phiadj' in makeplot:
        plotJ(sim,fitp['phiTJadjoint'],xKM,xp,fitp['EK'],fitp['EKpcolor'],vlim['j'],tInd,makeplot,'jadj',
                    '$\Phi_{adj}[E,x]$ Estimated diff. number flux from $LT\Phi=B$',
                     spfid,progms)
    if 'padj' in makeplot:
        plotVER(sim,vfit['adj'],xKM,xp,zKM,zp,vlim['p'],tInd,makeplot,'vadj',
              '$\hat{P}_{adj}$ from unfiltered backprojection $A^+b$',
              spfid,progms)
    if 'badj' in makeplot:
#        plotB(drn,cam,tInd,'$b$',cord,
#              'oneperplot')
#        plotB(dhat['fit_adj'],cam,tInd,'$b_{fit_adj}$',cord[::-1],
#              axm,car,cac,['twinx','oneperplot'])
        plotBcompare(drn,dhat['fit_adj'],cam,sim.nCamUsed,nCutPix,
                     bcomptxt,'adj',cord,vlim['b'],tInd,makeplot,progms)
#    if 'vbackproj' in makeplot:
#        plotVER(sim,Phat['vBackProj'],xKM,xp,zKM,zp,vlim['p'],tInd,makeplot,'vadj',
#              '$\hat{v}_{back_projection}$ from unfiltered backprojection $A^+b$',
#              spfid,progms)
#%% ############################################################################
def plotnoise(cam,tInd,makeplot,prefix,progms):
    fg = figure()
    ax = fg.add_subplot(311)

    for c in cam:
        ax.plot(cam[c].dnoise,label=cam[c].name)
        ax.set_ylabel('amplitude')
        #ax.set_xlabel('pixel number')
        ax.set_title('Noise injected into raw intensity data')
        ax.grid(True)
    ax.legend(loc='best')

    ax2 = fg.add_subplot(312)
    for c in cam:
        ax2.plot(cam[c].noisy,label=cam[c].name)
        ax2.set_ylabel('amplitude')
        ax2.set_xlabel('pixel number')
        ax2.set_title('Noisy data')
        ax2.grid(True)
    ax2.legend(loc='best')

    ax3 = fg.add_subplot(313)
    for c in cam:
        ax3.plot(cam[c].nonneg,label=cam[c].name)
        ax3.set_ylabel('amplitude')
        ax3.set_xlabel('pixel number')
        ax3.set_title('Noisy, zero-truncated data')
        ax3.grid(True)
    ax3.legend(loc='best')


    writeplots(fg,prefix,tInd,makeplot,progms,format1d)

def plottphi0(Tm,Phi0,Jxi,Ek,zKM,vlim,sim,tInd,makeplot,prefix,progms):
    fg = figure()
    ax = fg.gca()

    for i,e in enumerate(Ek):
        #try:
            ax.plot(Tm[:,i]*Phi0[i,Jxi], zKM, label='{:0.0f}'.format(e) + 'eV')
        #except ValueError:
        #    print('skipped plotting {:0.0f}'.format(e) +' eV due to zero diff. num flux at that energy')
    ax.set_xscale('log')
    ax.set_title('individual eigenprofile contributions to $P_{fwd}$',fontsize=tfs)
    ax.set_xlabel('volume emission rate [photons cm$^{-3}$ s$^{-1}$]',fontsize=afs)
    ax.set_ylabel('$B_\parallel$ [km]',fontsize=afs)
    ax.grid(True)

    ax.legend(loc='upper left',framealpha=0.5)

    ax.yaxis.set_major_locator(MultipleLocator(dymaj))
    ax.yaxis.set_minor_locator(MultipleLocator(dymin))
    ax.tick_params(axis='both', which='both', direction='in', labelsize=tkfs)
    ax.set_xlim(vlim[4:])
    ax.set_ylim(vlim[2:-2])

    writeplots(fg,prefix,tInd,makeplot,progms)

def ploteig(EKpcolor,zKM,Tm,vlim,sim,tInd,makeplot,prefix,progms):

    fg = figure(); ax = fg.gca()
    pcm = ax.pcolormesh(EKpcolor, zKM, Tm,
                        edgecolors='none',#cmap=pcmcmap,
                        norm=LogNorm(),
                        vmin=vlim[4], vmax=vlim[5])
    ax.set_xlabel('Energy [eV]',fontsize=afs)
    ax.set_ylabel('$B_\parallel$ [km]',fontsize=afs)
    ax.autoscale(True,tight=True)
    ax.set_xscale('log')
    ax.yaxis.set_major_locator(MultipleLocator(dymaj))
    ax.yaxis.set_minor_locator(MultipleLocator(dymin))
    mptitle = '$P_{eig}$, filter: ' + sim.opticalfilter
    mptitle += str(sim.reacreq)
#    if sim.loadver:
#        mptitle+= '\n' + sim.loadverfn
#    else:
#        mptitle += '\n' + sim.transcarpath
    ax.set_title(mptitle,fontsize=tfs)
    cbar = fg.colorbar(pcm,ax=ax)
    cbar.set_label('[photons cm$^{-3}$s$^{-1}$]',labelpad=0,fontsize=afs)
    cbar.ax.tick_params(labelsize=afs)
    cbar.ax.yaxis.get_offset_text().set_size(afs)

    ax.tick_params(axis='both', which='both', direction='out', labelsize=tkfs)
    ax.set_ylim(vlim[2:4])
    writeplots(fg,prefix,tInd,makeplot,progms)

def ploteig1d(Ek,zKM,Tm,vlim,sim,tInd,makeplot,prefix,progms):
    firstbeamind=0

    beamsel = s_[firstbeamind:-firstbeamind-1]
    #print(len(EK))
    fg = figure()
    ax = fg.gca()

    for i,e in enumerate(Ek[beamsel],firstbeamind):
        ax.semilogx(Tm[:,i], zKM,label='{:0.0f}'.format(e)+'eV')#,marker='.')
    ax.yaxis.set_major_locator(MultipleLocator(dymaj))
    ax.yaxis.set_minor_locator(MultipleLocator(dymin))
    ax.set_ylabel('$B_\parallel$ [km]',fontsize=afs)
    ax.set_xlabel('volume emission rate [photons cm$^{-3}$s$^{-1}$]',fontsize=afs)
    ax.tick_params(axis='both', which='both', direction='in', labelsize=tkfs)
    ax.set_xlim(vlim[4:])
    ax.set_ylim(vlim[2:-2])
    titletxt = '$P_{eig}$ '
    titletxt += str(sim.reacreq)

    titletxt += ' Opt. Filter: ' + sim.opticalfilter
    ax.set_title(titletxt, fontsize=tfs)
    ax.legend(loc='upper left',framealpha=0.5)
    #ax.autoscale(True,tight=True)
    ax.grid(True)
    writeplots(fg,prefix,tInd,makeplot,progms)


def plotPlainImg(sim,cam,rawdata,t,makeplot,prefix,titletxt,figh,spfid,progms):
    """http://stackoverflow.com/questions/22408237/named-colors-in-matplotlib"""
    for c in cam:
        figure(figh).clf()
        fg = figure(figh)
        ax = fg.gca() #Axes(fg,[0,0,1,1])
        ax.set_axis_off()
        ax.imshow(rawdata[c][t,:,:],
                  origin='lower',
                  vmin=cam[c].plotminmax[0], vmax=cam[c].plotminmax[1],
                  cmap='gray')
        ax.text(0.05, 0.075, cam[c].tKeo[t].strftime('%Y-%m-%dT%H:%M:%S.%f')[:23],
                     ha='left',
                     va='top',
                     transform=ax.transAxes,
                     color='limegreen',
                     #weight='bold',
                     size=24
                    )
        if in1d(('rawpng','save'),makeplot).any():
            writeplots(fg,'cam{}rawFrame'.format(c),t,'png',progms)
        else:
            draw()

#%%
def plotRealImg(sim,cam,rawdata,t,makeplot,prefix,titletxt,figh,spfid,progms):
    showcb = True
    #alltReq = sim.alltReq
    figure(figh).clf()
    fg,axm = subplots(nrows=1,ncols=2,sharey='row',num=figh, figsize=(11.5,5))
    for c,ax in zip(cam,axm):
        #fixme this would need help if one of the cameras isn't plotted (this will probably never happen)

        #plotting raw uint16 data
        hi = ax.imshow(rawdata[c][t,:,:],
                         origin='lower',
                         vmin=cam[c].plotminmax[0], vmax=cam[c].plotminmax[1],
                         cmap='gray')
        if showcb:
            hc = fg.colorbar(hi, ax=ax) #not cax!
            hc.set_label(str(rawdata[c].dtype) + ' data numbers')
        ax.set_title('Cam{} time: {}'.format(c,cam[c].tKeo[t]) ,fontsize=9)
        #ax.set_xlabel('x-pixel')
        if c=='0':
            ax.set_ylabel('y-pixel') #save horizontal space
    #%% plotting 1D cut line
        ax.plot(cam[c].cutcol,cam[c].cutrow,
                 marker='.',linestyle='none',markersize=1)
        #plot magnetic zenith
        ax.plot(cam[c].cutcol[cam[c].angleMagzenind],
                cam[c].cutrow[cam[c].angleMagzenind],
                marker='*',linestyle='none',color='red',markersize=10)
    #%% plot cleanup
        ax.autoscale(True,tight=True) #fills existing axes


    if in1d(('rawpng','save'),makeplot).any():
        writeplots(fg,'rawFrame',t,'png',progms)
    else:
        draw() #This draw must be here for multiplot animations or this window won't update!


def plotPicard(A,b,cvar=None):
    from picard import picard
    from scipy.sparse import issparse
    from numpy.linalg import svd
    print('computing SVD for Picard plot: ' + cvar + '  ...',end='')
    if issparse(A):
        [U,s,V] = svd(A.todense(order='F'),full_matrices=False,compute_uv=1) #V is not used!
    else:
        [U,s,V] = svd(A,full_matrices=False,compute_uv=1) #V is not used!
    picard(asfortranarray(U), s, b)
    print('DONE')

def plotJ1D(sim,PhiFwd,PhiInv,Ek,vlim,tInd,makeplot,prefix,titletxt,spfid,progms):
    anno = False
    #fontsizes

    fg = figure()
    ax = fg.gca()

    for Phi,l in zip((PhiFwd,PhiInv),('$\Phi$','$\widehat{{\Phi}}$')):
        if Phi is not None:
            iE0est = Phi.argmax()
            E0est = Ek[iE0est]

            try:
                ax.loglog(Ek,Phi,marker='.',
                          label= l+ ', {:0.0f} eV'.format(E0est))
                ax.axvline(E0est,linestyle='--',color='red')
                if anno:
                    ax.annotate(('$\hatE_0$={:0.0f}'.format(E0est) + ' eV'),
                             xy=(E0est,Phi[iE0est]),
                             xytext=(E0est*0.75, Phi[iE0est]*0.2),
                             ha='right', va='top',
                             arrowprops={'facecolor':'red', 'shrink':0.2, 'headwidth':8, 'width':2},)
            except ValueError as e:
                print('** could not plot Jfwd1D due to non-positive value(s) in Jflux, t= '+str(tInd) + ' ' + titletxt)
                print('* did you pick the correct --x1d ?   {}'.format(e))

    ax.grid(True)
    ax.autoscale(True,tight=False)
    ax.set_ylim(1, vlim[1])
    ax.set_xlim([Ek[0]*0.98, Ek[-1]*1.05])

    ax.set_xlabel('particle energy [eV]', fontsize=afs,labelpad=0)
    ax.set_ylabel('Differential Number Flux  [cm$^{-2}$s$^{-1}$eV$^{-1}$]',fontsize=afs,labelpad=0)
    ax.set_title(titletxt,fontsize=tfs)

    ax.tick_params(axis='both', which='both',labelsize=tkfs)

    ax.legend(loc='lower left')

    writeplots(fg,prefix,tInd,makeplot,progms,format1d)

    if 'h5' in makeplot:
        dumph5(sim,spfid,prefix,tInd,PhiFwd=PhiFwd,PhiInv=PhiInv,Ek=Ek)
#%%
def plotJ(sim,Jflux,x,xp,Ek,EKpcolor,vlim,tInd,makeplot,prefix,titletxt,spfid,progms):

    cnorm,sfmt = logfmt(makeplot)
    plt3 = 'octave'         #'octave' #'mpl'#'mayavi'
    #fontsizes

#%% 2-D
    pre = ('lin','log')
    vmin = (0.,1.)
    for c,s,v,p in zip(cnorm,sfmt,vmin,pre):
        #determine lowest level to plot
        vmin = getmin(v,vlim[0])

        if 'plotly' in makeplot:
            dpy = Data([Contour(x=x,y=Ek,z=Jflux,
                                xaxis='$B_\perp$ [km]',
                                yaxis='$B_\parallel$ [km]')])
            dlay= Layout(xaxis=XAxis(range=[vmin,vlim[1]]),
                         yaxis=YAxis(range=vlim[2:4]))

            dfg = Figure(data=dpy,layout=dlay)

            py.plot(dfg, filename='{}_{}_{}'.format(progms,prefix,tInd),auto_open=False)
            #print(plot_url)
        else:
            fg = figure()
            ax = fg.gca()
            if pstyle == 'pcolor':
                hc = ax.pcolormesh(xp,EKpcolor,Jflux,edgecolors='none',
                       norm=c, rasterized=False, #cmap=pcmcmap,
                       vmin=vmin, vmax=vlim[1]) #vmin can't be 0 when using LogNorm!
                       #my recollection is that rasterized=True didn't really help savefig speed
            elif pstyle=='contour':
                hc = ax.contour(x,Ek,Jflux,linewidths=2,
                                norm=c,vmin=v,vmax=vlim[1],)
            cbar = fg.colorbar(hc,ax=ax,format=s)
            cbar.set_label('[cm$^{-2}$s$^{-1}$eV$^{-1}$]', labelpad=0, fontsize=afs)
            cbar.ax.tick_params(labelsize=afs)
            cbar.ax.yaxis.get_offset_text().set_size(afs)
            cbar.ax.yaxis.get_offset_text().set_position((-1, 0))
            ax.set_yscale('log')

            fg.subplots_adjust(top=0.85)
            doJlbl(ax,titletxt)
            writeplots(fg,prefix+p,tInd,makeplot,progms)

#%% 3-D
    if '3d' in makeplot:
        try:
            ax3 = plotJ3(xp,EKpcolor,Jflux,plt3)
        except Exception as e:
            print('* falling back to matplotlib.  {}'.format(e))
            ax3 = plotJ3(xp,EKpcolor,Jflux,'mpl')
        doJlbl(ax3,titletxt)

    if 'h5' in makeplot:
        dumph5(sim,spfid,prefix,tInd,jfwd=Jflux,xp=xp,Ek=EKpcolor)

def doJlbl(ax,titletxt):
    ax.set_ylabel('Energy [eV]',fontsize=afs,labelpad=0)
    ax.set_xlabel('$B_\perp$ [km]',fontsize=afs,labelpad=0)
    ax.set_title(titletxt, fontsize=tfs, y = 1.02) # + '  $t_i=' + str(tInd) + '$'

    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.2))

    ax.tick_params(axis='both', which='both', direction='out',labelsize=afs)
    ax.autoscale(True,tight=True) #need this to fill axes (does not resize axes)
#%%
def plotJ3(xp,EKpcolor,Jflux,plt3):
    if plt3 == 'mpl':
        x,e = meshgrid(xp[:-1],EKpcolor) #TODO don't use corners here
        ax3 = figure().gca(projection='3d')
        ax3.plot_wireframe(x,e,Jflux)
        #ax3.yaxis.set_scale('log')
    elif plt3=='mayavi':
        ax3 = None
        from mayavi import mlab
        mlab.figure()
        # this .axes command not working
        #axm = mlab.axes(xlabel='x [km]',ylabel='Beam Energy [eV]',x_axis_visibility=True, y_axis_visibility=True)

        #mobj = mlab.mesh(e,x,Jflux, representation='wireframe')
        #mobj.actor.actor.scale = [1,0.0001,1]
        # http://docs.enthought.com/mayavi/mayavi/auto/mlab_helper_functions.html#mayavi.mlab.mesh
        # surf is preferred over mesh in this case

        mlab.surf(xp[:-1],EKpcolor,Jflux, warp_scale='auto',
                  representation='wireframe')
    elif plt3=='octave':
        from oct2py import octave
        ax3=None
        x,e = meshgrid(xp[:-1],EKpcolor) #TODO don't use corners here
        octave.Jmesh(x,e,Jflux)
    return ax3

def getmin(vs,vu):
    if vu is not None:
        return max(vs,vu)
    else:
        return vs

def plotVER1D(sim,pfwd,pinv,zKM,vlim,tInd,makeplot,prefix,titletxt,figh,spfid,progms):
    figure(figh).clf()
    fg = figure(figh)
    ax = fg.gca()

    for p,l in zip((pfwd,pinv),('fwd','est')):
        if p is not None:
            ax.semilogx(p,zKM,label=l)

    ax.legend(loc='upper right')
    ax.grid(True)
    ax.set_xlabel('Volume emission rate [photons cm$^{-3}$s$^{-1}$]',fontsize=afs,labelpad=0)
    ax.set_ylabel('$B_\parallel$ [km]',fontsize=afs,labelpad=-1)

    if not sim.loadver:
        titletxt += '\nReactions: {}'.format(sim.reacreq)

    ax.set_title(titletxt, fontsize=tfs, y=1.02)

    ax.yaxis.set_major_locator(MultipleLocator(dymaj))
    ax.yaxis.set_minor_locator(MultipleLocator(dymin))
    ax.tick_params(axis='both', which='major', direction='in',labelsize=tkfs)

    ax.set_ylim(vlim[:2])
    ax.set_xlim(vlim[2:])

    writeplots(fg,prefix,tInd,makeplot,progms,format1d)

    if 'h5' in makeplot: #a separate stanza
        dumph5(sim,spfid,prefix,tInd,pfwd=pfwd,pinv=pinv,z=zKM)
#%%
def plotVER(sim,ver,x,xp,z,zp,vlim,tInd,makeplot,prefix,titletxt,spfid,progms):

    cnorm,sfmt = logfmt(makeplot)
    '''
    # http://stackoverflow.com/questions/13943217/how-to-add-colorbars-to-scatterplots-created-like-this
    # http://stackoverflow.com/questions/18856069/how-to-shrink-a-subplot-colorbar
    '''



    vmin=(0.,1.)
    pre = ('lin','log')

    if ver is not None:
        for c,s,v,p in zip(cnorm,sfmt,vmin,pre):
          #determine lowest level to plot
          vmin = getmin(v,vlim[0])

          if 'plotly' in makeplot:
            dpy = Data([Contour(x=x,y=z,z=ver)])
            dlay= Layout(xaxis=XAxis(range=[vmin,vlim[1]], title='$B_\perp$ [km]'),
                         yaxis=YAxis(range=vlim[2:4],title='$B_\parallel$ [km]'))

            dfg = Figure(data=dpy,layout=dlay)

            py.plot(dfg, filename='{}_{}_{}'.format(progms,prefix,tInd), auto_open=False)
            #print(plot_url)
          else:
            fg = figure()
            ax = fg.gca()
            if pstyle=='pcolor':
                hc = ax.pcolormesh(xp,zp,ver,edgecolors='none',
                                    norm=c, #cmap=pcmcmap,n
                                    vmin=v,vmax=vlim[5],
                                    rasterized=True)
                ax.autoscale(True,tight=True) #need this to fill axes (does not resize axes)

            elif pstyle=='contour':
                print(p)
                hc = ax.contour(x,z,ver,linewidths=2,norm=c,vmin=vmin,vmax=vlim[5],)
                               # cmap=sns.dark_palette("palegreen", as_cmap=True))
                ax.clabel(hc,fontsize=6,inline=True)
                #hc.ax.get_children().set_linewidths(5.0)

            cbar = fg.colorbar(hc,ax=ax,format=s)
            cbar.set_label('[photons cm$^{-3}$s$^{-1}$]',labelpad=0,fontsize=afs)
            cbar.ax.tick_params(labelsize=afs)
            cbar.ax.yaxis.get_offset_text().set_size(afs)
            cbar.ax.yaxis.get_offset_text().set_position((-1, 0))
            if pstyle=='contour':
                cbar.add_lines(hc)

            ax.yaxis.set_major_locator(MultipleLocator(dymaj))
            ax.yaxis.set_minor_locator(MultipleLocator(dymin))
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.1))
            ax.tick_params(axis='both', which='both', direction='out', labelsize=tfs) #not needed with seaborn

            ax.set_xlabel('$B_\perp$ [km]',fontsize=afs)
            ax.set_ylabel('$B_\parallel$ [km]',fontsize=afs)

            ax.set_xlim(vlim[:2])
            ax.set_ylim(vlim[2:4])
            ax.set_title(titletxt,fontsize=tfs)
            writeplots(fg,prefix+p,tInd,makeplot,progms)
    else:
        text(0,0,'Using Actual Data (ver=None)')

    if 'h5' in makeplot: #a separate stanza
        dumph5(sim,spfid,prefix,tInd,ver=ver,x=x,xp=xp,z=z,zp=zp)
#%%
def plotBcompare(braw,bfit,cam,nCam,nCutPix,prefix,fittxt,cord,
                                                vlim,tInd,makeplot,progms):
    dosubtract = False

    fg = figure()
    ax1 = fg.gca()

#%% plot raw
    cnorm,sfmt = logfmt(makeplot,(-3,3))
    ax1.get_yaxis().set_major_formatter(sfmt[0]) #only need lin

    ax1.yaxis.get_offset_text().set_size(afs)
    ax1.tick_params(axis='both', which='both', direction='out',labelsize=afs)

    icm = 0
    for c in cam:
        cInd = cam[c].ind
        ax1.plot(cam[c].angle_deg,braw[cInd],
                 label=('$\mathbf{{B}}_{{camP{:d}}}$'.format(c)),
                 color=cord[icm])#, marker='.')
        icm+=1
#%% plot fit
    # do we need twinax? Let's find out if they're within factor of 10
    maxfit = bfit.max(); maxraw = braw.max()
    if 10*maxraw > maxfit > 0.1*maxraw:
        singax = True
        ax2 = ax1
        ax1.set_ylabel('$\mathbf{B}$ [photons sr$^{-1}$ s$^{-1}$]',labelpad=0,fontsize=afs)
    else:
        singax = False
        ax2 = ax1.twinx()
        ax2.get_yaxis().set_major_formatter(sfmt)
        ax1.set_ylabel('$\mathbf{B}$ [photons sr$^{-1}$ s$^{-1}$]',labelpad=0,fontsize=afs)
        ax2.set_ylabel('$\mathbf{\widehat{B}}$ [photons sr$^{-1}$ s$^{-1}$]',labelpad=0,fontsize=afs)
#%% now plot each camera
    for c in cam:
        cInd = cam[c].ind
        ax2.plot(cam[c].angle_deg,bfit[cInd],
                 label=('$\mathbf{{\widehat{{B}}}}_{{cam{}}}$'.format(c)),
                 color=cord[icm])#, marker='.')
        icm+=1
    if singax:
        ax1.legend(loc='upper left', fontsize=afs)
    else:
        '''
        legend for twinx http://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend
        '''
        h1, l1 = ax1.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax1.legend(h1+h2, l1+l2, loc='upper right')

    ax1.set_title('$\mathbf{B}$ ground-observed intensity', fontsize=tfs,y=1.04)

    ax1.set_xlabel('local view angle [deg.]',fontsize=afs)
    ax1.xaxis.set_major_locator(MultipleLocator(1))

    ax1.grid(True)
    ax1.autoscale(True,tight=True)
    print(vlim)
    ax1.set_ylim(vlim) #must come AFTER autoscale!
#%% do more detailed comparison
    if dosubtract:
        bias=[]
        ax3 =figure().gca()
        for c in cam:
            cInd = cam[c].ind
            bias.append(bfit[cInd].max() - braw[cInd].max())
            ax3.plot(cam[c].angle_deg,braw[cInd] - (bfit[cInd]))#-bias[iCam]))

        ax3.set_title('error $\mathbf{B}_{fwdf} - ' + fittxt + ', bias={}'.format(bias),fontsize=12)
        #  $t_i=' + str(tInd) + '$'
        ax3.set_xlabel('local view angle [deg.]')
        ax3.xaxis.set_major_locator(MultipleLocator(1))
        ax3.set_ylabel('error',labelpad=0)
        ax3.get_yaxis().set_major_formatter(sfmt)


    writeplots(fg,prefix,tInd,makeplot,progms,format1d)

#%%
def plotB(bpix,isrealdata,cam,nCutPix,vlim,tInd,makeplot,labeltxt,cord,progms):

    cnorm,sfmt = logfmt(makeplot)

    fgb = figure()
    ax1 = fgb.gca()

    std = []
#%% do we need twinax? Let's find out if they're within factor of 10
#    thiscammax = empty(len(cam))
#    for c in cam:
#        thiscammax[c] = bpix[cam[c].ind].max()
#    maxraw = thiscammax.max(); minmaxraw = thiscammax.min()
#    if 10*maxraw > minmaxraw > 0.1*maxraw:
#        singax = True
#        ax2 = ax1
#        ax1.set_ylabel('$b$ Data Numbers',labelpad=0,fontsize=afs)
#    else:
#        singax = False
#        ax2 = ax1.twinx()
#        ax2.get_yaxis().set_major_formatter(sfmt)
#        ax1.set_ylabel('$b_{raw}$ data numbers',labelpad=0,fontsize=afs)
#        ax2.set_ylabel( fittxt + ' Data Numbers',labelpad=0,fontsize=afs)
#%% make plot
    for c in cam:
        std.append('{:0.1e}'.format(cam[c].noiselam))
        cInd = cam[c].ind
        ax1.plot(cam[c].angle_deg,
                bpix[cInd],
                 label=(labeltxt + ',cam{}$'.format(c)),
                 #marker='.',
                 color=cord[c])
    doBlbl(ax1,isrealdata,sfmt[0],vlim,labeltxt,std) #b is never log
    writeplots(fgb,'b'+labeltxt[4:-2],tInd,makeplot,progms,format1d)

def doBlbl(axb,isrealdata,sfmt,vlim,labeltxt,noiselam):
    axb.legend(loc='upper left',fontsize=afs)
    axb.set_title('Observer Intensity',fontsize=tfs)

    axb.yaxis.get_offset_text().set_size(afs)

#    if not isrealdata and noiseStd:
#        ax.text(0.05,0.95,'noise std. dev.=' +noiseStd,fontsize=afs,
#                transform=ax.transAxes,va='top',ha='left')
    axb.set_xlabel('local view angle [deg.]',fontsize=afs)

    axb.xaxis.set_major_locator(MultipleLocator(1))
    axb.xaxis.set_minor_locator(MultipleLocator(0.2))
    axb.tick_params(axis='both', which='both', direction='out',labelsize=afs)

    #sfmt.set_powerlimits((-6, 6))
    axb.yaxis.set_major_formatter(sfmt)

    axb.set_ylabel(labeltxt + '}$ [photons sr$^{-1}$ s$^{-1}$]',labelpad=-1,fontsize=afs) # yes the }$ is appropriate

    axb.grid(True,which='major')
    axb.autoscale(True,tight=True)
    axb.set_ylim(vlim) #must come AFTER autoscale!

def allLinePlot(simfilelist):
#%% priming
    nsim = len(simfilelist)
    nt = len(simfilelist[0])
    with h5py.File(simfilelist[0][0],'r',libver='latest') as fid:
        J0 = fid['/Jn/j0'].value
        nE = J0.shape[0]
        Ek = fid['/EK'].value

    Jme = empty((nsim,nE,nt),dtype=float,order='F')
#%% working
    for ti in range(nt):
        for si,fnlist in enumerate(simfilelist):
            with h5py.File(fnlist[ti],'r',libver='latest') as fid:
                Jme[si,:,ti] = fid['/Jn/jme_1D'].value

        ax = figure(ti).gca()
        ax.plot(Ek,Jme[:,:,ti])
#%%
def planview3(cam,xKM,zKM,makeplot,figh,progms):
    from coordconv3d import aer2ecef, geodetic2ecef
    fg = figure(figh)#; fg.clf()
    ax = fg.gca(projection='3d')
    az = 106.444916022 #TODO arbitrary should be from real camera pointing
    srange = 500e3 #[m]
    decimaterayfactor = 16
    clr=['r','b','g','m']

    for c in cam:
        el =  cam[c].angle_deg[::decimaterayfactor] #yes double colon
        Np = el.size
        x0,y0,z0 = geodetic2ecef(cam[c].lat, cam[c].lon, cam[c].alt_m)
        # get LLA of pixel rays at 100km altitude
        xray,yray,zray = aer2ecef(az,el,srange,
                                  cam[c].lat, cam[c].lon, cam[c].alt_m)
        #camera rays
        for cri in range(Np):
            ax.plot((x0,xray[cri]),(y0,yray[cri]),(z0,zray[cri]),color=clr[c])

    #plot sphere
    earthrad = 6371e3 #[m]
    Nps = 50
    u = linspace(0, 2 * pi, Nps)
    v = linspace(0, pi, Nps)

    x = earthrad * outer(cos(u), sin(v))
    y = earthrad * outer(sin(u), sin(v))
    z = earthrad * outer(ones_like(u), cos(v))
    ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='g')
#%%
def planviewkml(cam,xKM,zKM,makeplot,figh,progms):
    """
    https://developers.google.com/kml/documentation/cameras
    """
    from coordconv3d import aer2geodetic#, aer2ecef
    import simplekml as skml
    decimaterayfactor = 16

    az = 106.444916022 #TODO arbitrary should be from real camera pointing
    srange = 500e3 #[m]

    clr=['r','b','g','m']
    kclr = ['ff5c5ccd','ffff0000']
    ax = figure(figh).gca(); clf()

    cnorm,sfmt = logfmt(makeplot,(-6,6))


    kml1d = skml.Kml()
    #setup camera (I preferred LookAt)
#        camview = skml.Camera(latitude=67.2,
#                              longitude=-147.2,
#                              altitude=100e3, #meters, int
#                              heading = 180., tilt = 10., roll = 0.,
#                              altitudemode = skml.AltitudeMode.relativetoground)
    #setup LookAt
    lkat = skml.LookAt(latitude=65.111,
                       longitude=-147.465,
                       altitude=0,
                       heading=180,
                       range=4e3,
                       tilt=45)

    lla = []
    for c in cam:
        #az is the temporary scalar defined above FIXME
        el = cam[c].angle_deg[::decimaterayfactor] #double colon
        Np = el.size
        lat0 = cam[c].lat
        lon0 = cam[c].lon
        alt0 = cam[c].alt_m #[meters]
        lla.append((lon0,lat0))
        # get LLA of pixel rays at 100km altitude
        latre,lonre,altre = aer2geodetic(az,el,srange,lat0,lon0,alt0)
        # get ECEF of center ray at 90km (bottom of model space)
        #centazray = az #TODO
        #centelray = cam[ci].angle_deg[Np//2]
        #xrc,yrc,zrc = aer2ecef(centazray,centelray,zbottom,lat0,lon0,alt0)

        #camera location points
        bpnt = kml1d.newpoint(name='HiST{}'.format(c),
                              description='camera {} location'.format(c),
                              coords=[(lon0,lat0)],
                              altitudemode=skml.AltitudeMode.clamptoground)
        bpnt.style.iconstyle.icon.href='http://maps.google.com/mapfiles/kml/shapes/arrow.png' # 'http://maps.google.com/mapfiles/kml/paddle/pink-blank.png'
        bpnt.style.iconstyle.scale = 2.0
        bpnt.style.labelstyle.size= 2.5
        #bpnt.camera = camview
        bpnt.lookat = lkat

        #camera rays
        if 'kmlrays' in makeplot:
            for cri in range(Np):
                linestr = kml1d.newlinestring(name='')
                linestr.coords = [(lon0,       lat0,       alt0),
                                  (lonre[cri], latre[cri], altre[cri])]
                linestr.altitudemode = skml.AltitudeMode.relativetoground
                linestr.style.linestyle.color = kclr[c]


        #plot
        ax.plot(lonre,latre,'x',color=clr[c],markersize=6)
        ax.plot(lon0,lat0,'o',color=clr[c],markersize=12,label='cam{}'.format(c))


    ax.set_ylabel('WGS84 latitude [deg.]')
    ax.set_xlabel('WGS84 longitude [deg.]')
    ax.set_title('pixel locations at 100km altitude')
    ax.legend(fontsize=afs)

        #setup line on ground connecting sites
    """
    https://developers.google.com/kml/faq#linestyle
    https://simplekml.readthedocs.org/en/latest/geometries.html#simplekml.LineString
    """
    ls = kml1d.newlinestring(name='3 km', coords=lla)
    ls.style.linestyle.width=5
    ls.style.linestyle.color=skml.Color.yellow
    ls.style.labelstyle.scale=2.5
    ls.style.labelstyle.color = skml.Color.white
    ls.style.labelstyle.gxlabelvisibility=1
    ls.visiblity=1

    #if 'kml' in makeplot or 'kmlrays' in makeplot: #write KML
    kmlfn = join(progms,'cam'+'.kml')
    print('saving ' + kmlfn)
    kml1d.save(kmlfn)
#%%
def dumph5(sim,fn,prefix,tInd,**writevar): #used in other .py too
    print('dumping to '+ fn)
    with h5py.File(fn,'a',libver='latest') as f:
        for k,v in writevar.items():
            try:
                f['/'+k+'/'+prefix]=v
            except Exception as e:
                print('failed to write {}.  {}'.format(k,e))
#%%
#def writenoax(img,plotprefix,tInd,method,progms,minmax):
#    tmpl = ('eps','jpg','png','pdf')
#    used = in1d(tmpl,method)
#    if used.any():
#        fmt = array(tmpl)[used][0]
#        cn = join(progms,(plotprefix + '_t{:03d}'.format(tInd) + '.' + fmt))
#        print('saving ' + cn + '...')
#        imsave(cn,img,vmin=minmax[0],vmax=minmax[1],format=fmt,cmap='gray')

def writeplots(fg,plotprefix,tInd,method,progms,overridefmt=None):
    #TIF was not faster and was 100 times the file size!
    #PGF is slow and big file,
    #RAW crashes
    #JPG no faster than PNG

    tmpl = ('eps','jpg','png','pdf')
    used = in1d(tmpl,method)
    if used.any():
        if overridefmt is not None:
            fmt = overridefmt
        else:
            fmt = array(tmpl)[used][0]
        cn = join(progms,(plotprefix + '_t{:03d}.{:s}'.format(tInd,fmt)))
        print('save {}...'.format(cn),end='')
        fg.savefig(cn,bbox_inches='tight',dpi=plotdpi,format=fmt)  # this is slow and async
#%%
def getx0E0(jf,Ek,x,tInd,progms,makeplot):
    if jf is None or isnan(jf).any() or not 'optim' in makeplot:
        return nan, nan,nan,nan
#%% interp log to linear
    Eklin = linspace(Ek[0],Ek[-1],200)
    ff = interp1d(Ek,jf,kind='linear',axis=0)
    jflin = ff(Eklin)
#%% first guess of 2-D peak, take region of pixels near the peak to fit
    #note that unravel_index for jf must be order='C'
    pkrow,pkcol = unravel_index(jflin.argmax(axis=None),  jflin.shape,   order='C')
    nh = (5,20) #MxN region to extract
    rcpluck = pkrow-nh[0], pkcol-nh[1]
    gpix = jflin[rcpluck[0]:pkrow+nh[0],rcpluck[1]:pkcol+nh[1]]
    #set_trace()
    try:
        gparam = gaussfit(gpix,returnfitimage=False)
    except ValueError:
        print(gpix.shape)
        print(pkrow,pkcol)
        print('gaussian fit doesnt work when peak is against edge of model space')
        return nan, nan,nan,nan

    Ghcol = rint(gparam[2] + rcpluck[1])
    Ghrow = rint(gparam[3] + rcpluck[0]) + 1# altitude of peak

    gparamshift = gparam.copy()
    gparamshift[2] = gparam[2] + rcpluck[1]
    gparamshift[3] = gparam[3] + rcpluck[0]
    gpix = twodgaussian(gparamshift,shape=jflin.shape)
    try: #gx0, gE0 and the rest need to be inside this try:
        gx0 = x[Ghcol]
        gE0 = Eklin[Ghrow]

        if 'gfit' in makeplot:
            fgb = figure()
            axb = fgb.gca()
            hb = axb.pcolormesh(x,Eklin,jflin)
            fgb.colorbar(hb)
            axb.set_title('base')
            axb.set_xlabel('x [km]')
            axb.set_ylabel('Beam Energy [eV]')
            axb.set_yscale('log')
            axb.autoscale(True,tight=True)

            fg = figure()
            axg = fg.gca()
            hp = axg.pcolormesh(x,Eklin,gpix)
            fg.colorbar(hp)
            axg.set_title('gaussian fit: $x_{{gauss,0}},E_{{gauss0}}$={:0.2f},{:0.0f})'.format(gx0,gE0))
            axg.set_xlabel('x [km]')
            axg.set_ylabel('Beam Energy [eV]')
            axg.set_yscale('log')
            axg.autoscale(True,tight=True)

            writeplots(fgb,'fluxlininterp',tInd,makeplot,progms)
            writeplots(fg,'gaussfitlin',tInd,makeplot,progms)

        #Ghrow,Ghcol = np.unravel_index(gfit.argmax(axis=None),gfit.shape, order='C')

        return gx0, gE0, x[pkcol], Eklin[pkrow]
    except IndexError:
        print('gaussian fit was outside model space')
        return nan, nan, nan, nan

