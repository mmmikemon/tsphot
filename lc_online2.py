#!/usr/bin/env python 
import matplotlib
#matplotlib.use('TkAgg')
from matplotlib.pyplot import *
from matplotlib.transforms import offset_copy
import numpy as np
import scipy 
from scipy.signal import lombscargle
import scipy.optimize as sco
import matplotlib.gridspec as gridspec
from matplotlib.widgets import CheckButtons
from matplotlib.ticker import MaxNLocator

global iline, icts

def padxlim2(ax1,xdata):
    """
    Set limits on current axis object in a SM-like manner
    """
    pad = 0.05
    xmax = max(xdata)
    xmin = min(xdata)
    deltax = xmax - xmin
    dx = pad * deltax
    x0 = xmin - dx
    x1 = xmax + dx
    print x0, x1
    ax1.set_xlim(x0,x1)

def fwhm_fit2(aplist,targs):
    nlc = len(targs[:,0])
    fwhm_vec = []
    for j in np.arange(nlc):
        apvec = targs[j,:]
        f = apvec
        ndim = len(f)
        i = np.arange(ndim)
        ip = i + 1
        aplist2 = np.concatenate([ [0],aplist])
        apvec2 = np.concatenate([ [0],apvec])
        dapvec2 = (apvec2[ip]-apvec2[i])/(aplist2[ip]**2-aplist2[i]**2)
        aplist_shifted = (2./3.)*(aplist2[ip]**3-aplist2[i]**3)/(aplist2[ip]**2-aplist2[i]**2)
        s0 = 10.
        s1 = dapvec2[0]
        w = 3.
        pinitial = np.array([ s0, s1, w ])
        popt, pcov = sco.curve_fit(psf, aplist_shifted, dapvec2, p0=pinitial)
        w = popt[2]
        fwhm = 2. * np.sqrt(2.*np.log(2.)) * w
        fwhm_vec.append(fwhm)

    fwhm_vec = np.array(fwhm_vec)
    apfine = np.arange(0,max(aplist),0.1)
    psf_fit = psf((apfine),*popt)
    print 'FWHM = ',fwhm
    fig=figure(3,figsize=(7,4.7))
    ax1 = fig.add_subplot()
    plot(aplist_shifted,dapvec2,'o')
    plot(apfine,psf_fit,'-')
    xlabel('radius (pixels)')
    ylabel('profile')
    tstring = 'FWHM is {0:.3f} pixels\n'.format(fwhm)
    totstring = tstring 
    x1,x2=xlim()
    y1,y2=ylim()
    xs=0.5
    ys=0.8
    xpos = x1 + xs*(x2-x1)
    ypos = y1 + ys*(y2-y1)
    text(xpos,ypos,totstring, horizontalalignment='center', verticalalignment='center')

    tight_layout()
    #show()

    #psffile='psf_fit.pdf'
    #savefig(psffile,transparent=True,bbox_inches='tight')
    #close()
    #print 'PSF fit stored in',psffile,'\n'
    return fwhm_vec

def fwhm_fit(aplist,apvec,is_first_iter):
    f = apvec
    ndim = len(f)
    i = np.arange(ndim)
    ip = i + 1
    aplist2 = np.concatenate([ [0],aplist])
    apvec2 = np.concatenate([ [0],apvec])
    dapvec2 = (apvec2[ip]-apvec2[i])/(aplist2[ip]**2-aplist2[i]**2)
    aplist_shifted = 0.5*(aplist2[ip]+aplist2[i])
    s0 = 10.
    s1 = dapvec2[0]
    w = 3.
    pinitial = np.array([ s0, s1, w ])
    popt, pcov = sco.curve_fit(psf, aplist_shifted, dapvec2, p0=pinitial)
    w = popt[2]
    fwhm = 2. * np.sqrt(2.*np.log(2.)) * w
    apfine = np.arange(0,max(aplist),0.1)
    psf_fit = psf((apfine),*popt)
    print 'FWHM = ',fwhm
    fig=figure(3)
    ax1 = fig.add_subplot()
    plot(aplist_shifted,dapvec2,'o')
    plot(apfine,psf_fit,'-')
    xlabel('radius (pixels)')
    ylabel('profile')
    tstring = 'FWHM is {0:.3f} pixels\n'.format(fwhm)
    totstring = tstring 
    x1,x2=xlim()
    y1,y2=ylim()
    xs=0.5
    ys=0.8
    xpos = x1 + xs*(x2-x1)
    ypos = y1 + ys*(y2-y1)
    text(xpos,ypos,totstring, horizontalalignment='center', verticalalignment='center')

    if is_first_iter:
        show(block=False)
        pause(1.0)
    else:
        draw()
        pause(1.0)
    tight_layout()

    #psffile='psf_fit.pdf'
    #savefig(psffile,transparent=True,bbox_inches='tight')
    #close()
    #print 'PSF fit stored in',psffile,'\n'

# Gaussian functional form assumed for PSF fits
def psf((rr),s0,s1,w):
    elow = -50.
    #arg = - (rr/w)**2
    arg = - rr**2/(2.*w**2)
    arg[arg <= elow] = elow
    intensity = s0 + s1 * np.exp( arg )
    fwhm = 2. * np.sqrt(2.*np.log(2.)) * w
    return intensity

# Total integrated flux and fwhm for the assumed Gaussian PSF
def psf_flux(s0,s1,x0,y0,w):
    flux = np.pi * s1/np.abs(w)
    fwhm = 2. * np.sqrt(np.log(2.)/np.abs(w))
    return fwhm, flux

def comb_string(combs):
    ic = 0
    for i in combs:
        if ic == 0:
            cstring = str(i)
        else:
            cstring = cstring + ' + ' + str(i)
        ic = ic + 1
    return cstring

def list_powerset2(lst):
    return reduce(lambda result, x: result + [subset + [x] for subset in result], lst, [[]])

# Make a plot of light curve scatter versus aperture size
#def applot(aplist,sigvec,apmin,cstring):
def applot(apfile, aplist,sigvec,apmin,cstring):
    #apfile = 'aperture.pdf'
    #fig=figure(2)
    fig=figure(2,figsize=(7,4.7))
    ax1 = fig.add_subplot()
    mask = np.where(aplist > 0.)
    ap_filt = aplist[mask]
    sig_filt = sigvec[mask]
    iarg = np.argsort(ap_filt)
    ap_filt_s = ap_filt[iarg]
    sig_filt_s = sig_filt[iarg]
    #plot(ap_filt,sig_filt,'-o')
    plot(ap_filt_s,sig_filt_s,'-o')
    xlabel('Aperture size')
    ylabel('Scatter')
    tstring = 'optimal aperture is {0} pixels\n'.format(apmin)
    cstring = 'optimal comparison star =\n' + cstring
    totstring = tstring + '\n' + cstring
    x1,x2=xlim()
    y1,y2=ylim()
    xs=0.5
    ys=0.8
    xpos = x1 + xs*(x2-x1)
    ypos = y1 + ys*(y2-y1)
    text(xpos,ypos,totstring, horizontalalignment='center', verticalalignment='center')

    tight_layout()
    #show()

    #savefig(apfile,transparent=True,bbox_inches='tight')
    #close()
    #print 'Aperture optimization stored in',apfile
    
# Make a plot of the optimal light curve file
def lcplot(flc_pdf, time,target,comp,sky,fwhm_vec,cstring,is_first_iter):
#def lcplot(time,target,comp,sky,fwhm_vec,cstring):
    global iline, icts
    ratio = target/comp
    ratio_norm = ratio/np.mean(ratio) - 1.
    #scipy.convolve(y, ones(N)/N)
    nfilt=5
    ratio_norm_smooth = scipy.convolve(ratio_norm, np.ones(nfilt)/nfilt)
    time_smooth = scipy.convolve(time, np.ones(nfilt)/nfilt)
    target_norm = target/np.mean(target) - 1.
    comp_norm = comp/np.mean(comp) - 1.
    #sky_norm = sky/np.mean(sky) - 1.
    sky_norm = sky/np.amin(sky)
    ndata = len(time)
    dt = time[1]-time[0]
    fmax = 1./(2.*dt)
    #fmax = 0.01
    df = 1./(6.*(time[-1]-time[0]))
    #df = 1./(2.*(time[-1]-time[0]))
    fmin = df
    fmin = 0.01*df
    #df = (fmax-fmin)/1000.
    farray = np.arange(fmin,fmax,df)
    omarray = 2.* np.pi * farray
    pow = lombscargle(time,ratio_norm,omarray)
    pow = pow * 4/ndata
    amp = np.sqrt(pow)

    fig=figure(1,figsize=(11, 11))
    gs1 = gridspec.GridSpec(6, 1,hspace=0.0)
    ax1 = fig.add_subplot(gs1[0])
    l1, = ax1.plot(time,ratio_norm,'o')
    #ax1.set_xlim([-500,15000])
    ax1.yaxis.set_major_locator( MaxNLocator(nbins = 6) )
    ylabel(r'$\delta I/I$',size='x-large')
    leg1=ax1.legend(['target'],'best',fancybox=True,shadow=False,markerscale=0.0,frameon=False)
    setp(ax1.get_xticklabels(), visible=False)
    #padxlim2(ax1,time)

    ax2 = fig.add_subplot(gs1[1],sharex=ax1)
    l2, = plot(time_smooth[0:-nfilt],ratio_norm_smooth[:-nfilt],'o')
    ax2.yaxis.set_major_locator( MaxNLocator(nbins = 6, prune = 'upper') )
    ylabel(r'$\delta I/I$',size='x-large')
    leg=ax2.legend(['target smoothed'],'best',fancybox=True,shadow=False,markerscale=0.0,frameon=False)
    setp(ax2.get_xticklabels(), visible=False)

    ax3 = fig.add_subplot(gs1[2],sharex=ax1)
    l3, = plot(time,comp_norm,'o')
    ax3.yaxis.set_major_locator( MaxNLocator(nbins = 6, prune = 'upper') )
    ylabel(r'$\delta I/I$',size='x-large')
    comstring = 'comparison\n' + '= ' + cstring
    leg=ax3.legend([comstring],'best',fancybox=True,shadow=False,markerscale=0.0,frameon=False)
    setp(ax3.get_xticklabels(), visible=False)

    ax4 = fig.add_subplot(gs1[3],sharex=ax1)
    pl1 =ax4.plot(time,sky_norm,'o',label = 'Sky')
    l4, = pl1
    ax4.yaxis.set_major_locator( MaxNLocator(nbins = 5, prune = 'upper') )
    ylabel('Sky',size='large')
    xlabel('time (sec)',size='large')
    fwhm_string = 'FWHM={0:.3f}'.format(fwhm_vec[-1])
    ax6 = ax4.twinx()
    pl2 = ax6.plot(time, fwhm_vec, 'ro',label = fwhm_string)
    l5, = pl2
    pls = pl1 + pl2
    ax6.set_ylabel('FWHM', color='r')
    for tl in ax6.get_yticklabels():
        tl.set_color('r')

    labs = [l.get_label() for l in pls]
    leg=ax4.legend(pls, labs, 'best', frameon = False, ncol=2, handletextpad=0.0, columnspacing=0.0)

    gs2 = gridspec.GridSpec(4, 1,hspace=0.1)
    ax5 = fig.add_subplot(gs2[3])

    freqs = farray * 1.e+6
    plot(freqs,amp)
    xmax=max(freqs)
    xlim(0.,xmax)
    xlabel(r'Frequency ($\mu$Hz)',size='large')
    ylabel('Amplitude',size='large')
    leg=ax5.legend(['FT'],'best',fancybox=True,shadow=False,handlelength=0.0)
    leg.draw_frame(False)

    rax = axes([0.46, 0.905, 0.11, 0.07])
    lin = ' lines'
    cts = ' counts'
    check = CheckButtons(rax, (lin,' counts'), (False, False))
    iline = 1
    icts = 1
    lset = [l1, l2, l3, l4, l5]
    def stylefunc(label):
        global iline, icts
        if label == lin: 
            if iline == 1:
                for ll in lset:
                    ll.set_linestyle('-')
                    ll.set_marker('None')
            else:
                for ll in lset:
                    ll.set_linestyle('o')
                    ll.set_marker('o')
            iline = -1 * iline
        if label == cts: 
            if icts == 1:
                l3.set_ydata(comp)
                ax3.set_ylabel('Comp',size='large')
                l4.set_ydata(sky)
            else:
                l3.set_ydata(comp_norm)
                ax3.set_ylabel(r'$\delta I/I$',size='x-large')
                l4.set_ydata(sky_norm)
            ax3.relim()
            ax3.autoscale_view()
            ax4.relim()
            ax4.autoscale_view()
            icts = -1 * icts
        draw()

    check.on_clicked(stylefunc)

    #tight_layout()
    if is_first_iter:
        show(block=False)
    else:
        draw()

    #filebase = 'lc.pdf'
    #savefig(filebase,transparent=True,bbox_inches='tight')
    #close()
    #print 'Optimal light curve plot stored in',filebase,'\n'

# ysig is normalized so that it represents the point-to-point 
# scatter sigma_i, assuming uncorrelated, random noise
def scatter(lcvec):
    ndata = len(lcvec)
    ivec = np.arange(0,ndata-1)
    ivecp = ivec + 1
    dy    = lcvec[ivecp] - lcvec[ivec]
    ysig  = np.sqrt(np.dot(dy,dy)/(2.*(ndata-1.)))
    return ysig

#if __name__ == '__main__':
def main(args,is_first_iter):

    #lcfile = 'lightcurve.app'
    #lcfile = 'lightcurve_test.app'
    lcfile = 'tmp2.app'

    # Get number of stars
    f = open(args.flc,'r')
    line = f.readline()
    s = line.split()
    nstars = int(s[-1])
    f.close()

    print '\nThe Number of stars is',nstars

    cols=range(0,3*nstars+1)
    var = np.loadtxt(args.flc,usecols=cols)

    jdarr  = var[:,0]
    aparr  = var[:,1]
    fluxes = var[:,2:2+nstars]
    sky    = var[:,2+nstars:2+2*nstars]
    pos    = var[:,2+2*nstars:2+4*nstars]
    fwhm   = var[:,2+4*nstars:2+5*nstars]

    # Get the list of unique aperture sizes
    apset = set(aparr)
    aplist = list(apset)
    aplist = np.array(aplist)

    print '\nUsing the following apertures:'
    print aplist

    # Cyle through all possible combinations of comparison stars and choose the one 
    # that minimizes the point-to-point scatter ysig

    starlist =  np.arange(1,nstars)
    pset = list_powerset2(starlist)
    del pset[0]

    # Now that we've got the list of comparison stars, let's find the optimal aperture

    # create arrays to store lightcurves for different apertures
    # Use median aperture for this
    apmed = float(int(np.median(aplist)))
    mask = np.where(aparr == apmed)
    stars = fluxes[mask,:][0]
    target = stars[:,0]
    ntarg = len(target)
    nap   = len(aplist)
    ncom  = len(pset)
    nlcs = nap*ncom

    targs   = np.zeros((ntarg,nap))
    comps   = np.zeros((ntarg,nap,ncom))
    skyvals = np.zeros((ntarg,nap))
    ysigarr = np.ones((nap,ncom))
    ysigarr = 1000000.*ysigarr
    # counter for apertures
    iap = 0
    for app in aplist:

        mask = np.where(aparr == app)

        jd =  jdarr[mask]
        ap =  aparr[mask]
        stars = fluxes[mask,:][0]
        skys  = sky[mask,:][0]
        xpos  = pos[mask,:][0]
        fw =   fwhm[mask]

        ta = stars[:,0]

        # store target and sky lightcurves
        targs[:,iap]=ta
        skyvals[:,iap]=skys[:,0]

        # Loop over all possible comparison stars
        # counter for combination of comparison stars
        nc = 0
        for term in pset:
            compstar = 0.*ta
            for ii in term:
                compstar = compstar + stars[:,ii]
            
            # divide by comparison, normalize, and shift to zero
            compstar[compstar < 10.] = 10.
            ratio = ta/compstar
            ratio = ratio/np.mean(ratio) - 1.0
            ysig = scatter(ratio)

            #sigvec.append(ysig)
            ysigarr[iap,nc] = ysig
            # print iap, nc, term, ysig

            # store comparison lightcurve
            comps[:,iap,nc]=compstar
            nc = nc + 1

        iap = iap + 1

    # Find optimal aperture and its index
    sigmin = np.amin(ysigarr)
    isigmin = np.argmin(ysigarr)
    iarg = np.unravel_index(isigmin,(nap,ncom))
    iapmin, ncmin = iarg
    combs = pset[ncmin]

    apmin = aplist[iapmin]
    report  = '\nThe optimal aperture is {0} pixels, the point-to-point scatter is {1:0.5f},'.format(apmin,sigmin)
    print report
    print 'and the optimal comparison star combination is',pset[ncmin],'\n'

    # Store "combination string" for the optimal combination of comparison stars
    cstring = comb_string(combs)

    # Make plot of scatter versus aperture size
    sigvec = ysigarr[:,ncmin]
    applot(args.fap_pdf, aplist,sigvec,apmin,cstring)

    time = 86400.*(jd-jd[0])
    target = targs[:,iapmin]
    compstar = comps[:,iapmin,ncmin]
    sky0 = skyvals[:,iapmin]

    # Calculate scatter for composite comparison star
    ncompstar = compstar/np.mean(compstar) - 1.0
    comp_ysig = scatter(ncompstar)
    print 'Scatter: target=',sigmin,' Comparison=',comp_ysig

    # Do FWHM fits for star at all points of light curve.
    # Choose comparison or target star based on lower level of scatter.
    if comp_ysig <= sigmin:
        fwhm_vec = fwhm_fit2(aplist,comps[:,:,ncmin])   # FWHM of composite comparison star
    else:
        fwhm_vec = fwhm_fit2(aplist,targs)   # FWHM of target star

    # Make online plot of lightcurves, sky, and the FT
    lcplot(args.flc_pdf, time,target,compstar,sky0,fwhm_vec,cstring,is_first_iter)

    return None

if __name__ == '__main__':
    arg_default_map = {}
    arg_default_map['flc'] = "lightcurve.txt"
    arg_default_map['flc_pdf'] = "lc.pdf"
    arg_default_map['fap_pdf'] = "aperture.pdf"
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=("Read lightcurve file and create plots."))
    parser.add_argument("--flc",
                        default=arg_default_map['flc'],
                        help=(("Input fixed-width-format text file"
                               +" with columns of star intensities by aperture size.\n"
                               +"Default: {fname}").format(fname=arg_default_map['flc'])))
    parser.add_argument("--flc_pdf",
                        default=arg_default_map['flc_pdf'],
                        help=(("Output .pdf file with plots of the lightcurve.\n"
                              +"Default: {fname}").format(fname=arg_default_map['flc_pdf'])))
    parser.add_argument("--fap_pdf",
                        default=arg_default_map['fap_pdf'],
                        help=(("Output .pdf file with plot of scatter vs aperture size.\n"
                               +"Default: {fname}").format(fname=arg_default_map['fap_pdf'])))
    parser.add_argument("--verbose", "-v",
                        action='store_true',
                        help=("Print 'INFO:' messages to stdout."))
    args = parser.parse_args()
    if args.verbose:
        print "INFO: Arguments:"
        for arg in args.__dict__:
            print ' ', arg, args.__dict__[arg]
    if not os.path.isfile(args.flc):
        raise IOError(("File does not exist: {fname}").format(fname=args.flc))
    main(args=args)
