########################################################################
# Code to produce a variety of model lensing statistics.               #
########################################################################

import sys
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import quad
from scipy.interpolate import splev, splrep
from scipy.special import jn,jn_zeros

def main():

# Parameters
  doproj = True  # Model projected statistics (or angular statistics)
  if (doproj):
    opt = 1      # 1) Delta Sigma 2) w_p 3) xi
    rmin = 0.5   # Minimum projected separation of model [Mpc/h]
    rmax = 50.   # Maximum projected separation of model [Mpc/h]
    nrbin = 10   # Number of bins of model
    rlinlog = 2  # Space bins as 1) linear or 2) logarithmic
    pimax = 60.  # Limit of integral over Pi [Mpc/h]
  else:
    opt = 1      # 1) xi_+ 2) xi_- 3) gamma_t 4) w
    thmin = 0.03 # Minimum angular separation of model [degrees]
    thmax = 3.   # Maximum angular separation of model [degrees]
    nthbin = 20  # Number of bins of model
    thlinlog = 2 # Space bins as 1) linear or 2) logarithmic
    lmin = 1.    # Minimum multipole to pre-compute
    lmax = 1.e+5 # Maximum multipole to pre-compute
    nl = 500     # Number of multipole bins to pre-compute
# Cosmological model
  om = 0.286     # Omega_m
  h = 0.7        # Hubble parameter
  fb = 0.164335664336 # Baryon fraction
  sig8 = 0.82    # sigma_8
  b = 1.         # Galaxy bias
# Files containing the source and lens redshift distributions
  pzsourcefile = 'pz_sources_buzzard_zp0pt80_1pt00.dat'
  pzlensfile = 'pz_lenses_buzzard_lrg_zs0pt40_0pt60.dat'
# File containing the model power spectrum
  pkcambfile,npcamb,nkcamb = 'pkcambhalofit_zrange_buzzard.dat',21,922
  klinlog,kminmod,kmaxmod,nkmod = 2,1.e-4,0.999e+4,800
# Pre-compute model as a function of r
  rminmod,rmaxmod,nrmod = 1.e-4,1.e+3,70
# Limits of integral over k
  kminint,kmaxint = 1.e-4,0.999e+4
# astropy cosmology
  cosmo = FlatLambdaCDM(H0=100.,Om0=om)
# Read in redshift distributions
  zsourcearr,pzsourcearr,zminsource,zmaxsource = readnz(pzsourcefile)
  zlensarr,pzlensarr,zminlens,zmaxlens = readnz(pzlensfile)
  tcksource = splrep(zsourcearr,pzsourcearr)
  tcklens = splrep(zlensarr,pzlensarr)
# Read in power spectra
  kmod,zmod,pkzmod = readpkcamb(pkcambfile,npcamb,nkcamb,om,h,fb,sig8,klinlog,kminmod,kmaxmod,nkmod)
# Determine radial functions
  zsourcearr,fsourcearr = getfsource(0.,zmaxsource,zminsource,zmaxsource,tcksource,cosmo)
  zlensarr,flensarr = getflens(zminlens,zmaxlens,tcklens,cosmo)
# Compute projected lensing correlation functions
  if (doproj):
    z = 0.5*(zminlens+zmaxlens)
    pkmod = intpkzmod(z,zmod,pkzmod)
# Generate lensing function
    rbin,lensprojmod = getlensprojmod(opt,rmin,rmax,nrbin,rlinlog,pimax,rminmod,rmaxmod,nrmod,kminint,kmaxint,kmod,pkmod,om)
# Compute angular lensing correlation functions
  else:
# Determine angular power spectra
    larr,clshapearr = getclshapearr(opt,lmin,lmax,nl,zminlens,zmaxlens,zlensarr,flensarr,zminsource,zmaxsource,zsourcearr,fsourcearr,kmod,zmod,pkzmod,om,b,cosmo)
# Generate lensing function
    thdeg,lensshapemod = getlensangmod(opt,thmin,thmax,nthbin,thlinlog,lmin,lmax,larr,clshapearr)
  return

########################################################################
# Read in redshift distribution file.                                  #
########################################################################

def readnz(nzfile):
  print nzfile
  f = open(nzfile,'r')
  fields = f.readline().split()
  zmin,zmax,nz = float(fields[0]),float(fields[1]),int(fields[2])
  zarr,pzarr = np.empty(nz),np.empty(nz)
  for i in range(nz):
    fields = f.readline().split()
    zarr[i],pzarr[i] = float(fields[0]),float(fields[1])
  f.close()
  return zarr,pzarr,zmin,zmax

########################################################################
# Read in and interpolate power spectrum file.                         #
########################################################################

def readpkcamb(cambfile,npcamb,nkcamb,om,h,fb,sig8,klinlog,kminmod,kmaxmod,nkmod):
  print '\nReading in power spectra...'
  if (klinlog == 1):
    kmod = np.linspace(kminmod,kmaxmod,nkmod)
  else:
    kmod = np.logspace(np.log10(kminmod),np.log10(kmaxmod),nkmod)
  omp,hp,fbp,sig8p,zmod = np.empty(npcamb),np.empty(npcamb),np.empty(npcamb),np.empty(npcamb),np.empty(npcamb)
  nom,nh,nfb,nsig8,nz = 0,0,0,0,0
  karr,pkparr = np.empty(nkcamb),np.empty(shape=(npcamb,nkcamb))
  print cambfile
  f = open(cambfile,'r')
  for ip in range(npcamb):
    fields = f.readline().split()
    iom,ih,ifb,isig8,iz = int(fields[0]),int(fields[1]),int(fields[2]),int(fields[3]),int(fields[4])
    nom,nh,nfb,nsig8,nz = max(iom,nom),max(ih,nh),max(ifb,nfb),max(isig8,nsig8),max(iz,nz)
    omp[iom-1],hp[ih-1],fbp[ifb-1],sig8p[isig8-1],zmod[nz-iz] = float(fields[5]),float(fields[6]),float(fields[7]),float(fields[8]),float(fields[9])
    for ik in range(nkcamb):
      fields = f.readline().split()
      if (ip == 0):
        karr[ik] = float(fields[0])
      pkparr[ip,ik] = float(fields[1])
  f.close()
  if ((om < omp[0]) or (om > omp[nom-1]) or (fb < fbp[0]) or (fb > fbp[nfb-1]) or (h < hp[0]) or (h > hp[nh-1]) or (sig8 < sig8p[0]) or (sig8 > sig8p[nsig8-1])):
    print 'Parameters outside grid!!'
    print om,omp[0],omp[nom-1]
    print fb,fbp[0],fbp[nfb-1]
    print h,hp[0],hp[nh-1]
    print sig8,sig8p[0],sig8p[nsig8-1]
    sys.exit()
  iom,ih,ifb,isig8 = 0,0,0,0
  omf,hf,fbf,sig8f = 0.,0.,0.,0.
  if (nom > 1):
    dom = omp[1]-omp[0]
    iom = int(np.floor((1.0001*(om-omp[0]))/dom))
    omf = (om-omp[iom])/dom
  if (nh > 1):
    dh = hp[1]-hp[0]
    ih = int(np.floor((1.0001*(h-hp[0]))/dh))
    hf = (h-hp[ih])/dh
  if (nfb > 1):
    dfb = fbp[1]-fbp[0]
    ifb = int(np.floor((1.0001*(fb-fbp[0]))/dfb))
    fbf = (fb-fbp[ifb])/dfb
  if (nsig8 > 1):
    dsig8 = sig8p[1]-sig8p[0]
    isig8 = int(np.floor((1.0001*(sig8-sig8p[0])/dsig8)))
    sig8f = (sig8-sig8p[isig8])/dsig8
  pk1,pk2,pk3,pk4 = np.empty(2),np.empty(shape=(2,2)),np.empty(shape=(2,2,2)),np.empty(shape=(2,2,2,2))
  pkzmod,pkarr = np.empty(shape=(nz,nkmod)),np.empty(nkcamb)
  for iz in range(nz):
    for ik in range(nkcamb):
      for iom1 in range(2):
        for ih1 in range(2):
          for ifb1 in range(2):
            for isig81 in range(2):
              if (nom == 1):
                iom1 = 0
              if (nh == 1):
                ih1 = 0
              if (nfb == 1):
                ifb1 = 0
              if (nsig8 == 1):
                isig81 = 0
              ip = (iom+iom1)*nh*nfb*nsig8*nz + (ih+ih1)*nfb*nsig8*nz + (ifb+ifb1)*nsig8*nz + (isig8+isig81)*nz + iz
              pk4[iom1,ih1,ifb1,isig81] = pkparr[ip,ik]
      pk3[:,:,:] = pk4[:,:,:,0] + sig8f*(pk4[:,:,:,1]-pk4[:,:,:,0])
      pk2[:,:] = pk3[:,:,0] + fbf*(pk3[:,:,1]-pk3[:,:,0])
      pk1[:] = pk2[:,0] + hf*(pk2[:,1]-pk2[:,0])
      pkarr[ik] = pk1[0] + omf*(pk1[1]-pk1[0])
    pkzmod[iz,:] = np.interp(kmod,karr,pkarr)
  return kmod,zmod,pkzmod

########################################################################
# Interpolate P(k,z) file to P(k) at a given z.                        #
########################################################################

def intpkzmod(z,zmod,pkzmod):
  dz = zmod[1]-zmod[0]
  iz = int(np.floor((z-zmod[0])/dz))
  zf = (z-zmod[iz])/dz
  pkmod = np.empty(pkzmod.shape[1])
  pkmod[:] = pkzmod[iz,:] + zf*(pkzmod[iz+1,:]-pkzmod[iz,:])
  return pkmod

########################################################################
# Interpolate P(k,z) file to a value at a given (k,z).                 #
########################################################################

def intpkzmod2(k,z,kmod,dk,zmod,dz,pkzmod):
  iz = int(np.floor((z-zmod[0])/dz))
  zf = (z-zmod[iz])/dz
  ik = int(np.floor((np.log10(k)-np.log10(kmod[0]))/dk))
  kf = (k-kmod[ik])/(kmod[ik+1]-kmod[ik])
  pk1 = pkzmod[iz,ik] + zf*(pkzmod[iz+1,ik]-pkzmod[iz,ik])
  pk2 = pkzmod[iz,ik+1] + zf*(pkzmod[iz+1,ik+1]-pkzmod[iz,ik+1])
  pkmod = pk1 + kf*(pk2-pk1)
  return pkmod

########################################################################
# Construct radial selection function for sources.                     #
########################################################################

def getfsource(zmin,zmax,zminint,zmaxint,tcksource,cosmo):
  print 'Determining source radial function...'
  print 'zmin =',zmin,'zmax =',zmax
  dzarr = 0.01
  nzarr = int(np.rint((zmax-zmin)/dzarr)) + 1
  zsourcearr,fsourcearr = np.linspace(zmin,zmax,nzarr),np.empty(nzarr)
  zarr = np.linspace(min(zmin,zminint),max(zmax,zmaxint),1000)
  rarr = cosmo.comoving_distance(zarr).value
  for i in range(nzarr):
    z = zsourcearr[i]
    x = cosmo.comoving_distance(z).value
    z1,z2 = max(z,zminint),zmaxint
    ans,err = quad(fsourceint,z1,z2,args=(tcksource,zarr,rarr,x))
    fsourcearr[i] = ans
  return zsourcearr,fsourcearr

def fsourceint(zp,tcksource,zarr,rarr,x):
  xp = np.interp(zp,zarr,rarr)
  pz = splev(zp,tcksource)
  return pz*((xp-x)/xp)

########################################################################
# Construct radial selection function for lenses.                      #
########################################################################

def getflens(zmin,zmax,tcklens,cosmo):
  print 'Determining lens radial function...'
  print 'zmin =',zmin,'zmax =',zmax
  dzarr = 0.01
  nzarr = int(np.rint((zmax-zmin)/dzarr)) + 1
  zlensarr = np.linspace(zmin,zmax,nzarr)
  rarr = cosmo.comoving_distance(zlensarr).value
  drdzarr = 2997.9*cosmo.inv_efunc(zlensarr)
  pz = splev(zlensarr,tcklens)
  flensarr = pz/(rarr*drdzarr)
  return zlensarr,flensarr

########################################################################
# Compute projected lensing correlation functions.                     #
# lensopt: 1) Delta Sigma 2) w_p 3) xi                                 #
########################################################################

def getlensprojmod(lensopt,rmin,rmax,nbin,linlog,pimax,rminmod,rmaxmod,nrmod,kminint,kmaxint,kmod,pkmod,om):
  if (lensopt == 1):
    print '\nGenerating model for Delta Sigma...'
  elif (lensopt == 2):
    print '\nGenerating model for w_p...'
    print 'pimax =',pimax
  elif (lensopt == 3):
    print '\nGenerating model for xi...'
    rminmod,rmaxmod,nrmod = 0.9*rmin,1.1*rmax,100
  if (linlog == 1):
    dr = (rmax-rmin)/float(nbin)
    rbin = np.linspace(rmin+0.5*dr,rmax-0.5*dr,nbin)
  else:
    dr = (np.log10(rmax)-np.log10(rmin))/float(nbin)
    rbin = np.logspace(np.log10(rmin)+0.5*dr,np.log10(rmax)-0.5*dr,nbin)
  rmod,ximod = getxiarr(rminmod,rmaxmod,nrmod,kminint,kmaxint,kmod,pkmod)
  print '\nGenerating projected clustering...'
  lrminmod,lrmod = np.log(rminmod),np.log(rmod)
  dlrmod = lrmod[1]-lrmod[0]
  nlrint = 10
  rint = np.logspace(np.log10(rminmod),np.log10(rmaxmod),nlrint+1)
  lensmod = np.zeros(nbin)
  for i in range(nbin):
    r = rbin[i]
    if (lensopt == 1):
      term1,term2 = 0.,0.
      for j in range(nlrint):
        lr1,lr2 = np.log(rint[j]),np.log(rint[j+1])
        ans,err = quad(term1int,lr1,lr2,args=(lrmod,ximod,lrminmod,dlrmod))
        term1 += ans
        ans,err = quad(term2int,lr1,lr2,args=(r,lrmod,ximod,lrminmod,dlrmod))
        term2 += ans
      term1 *= 4./(r**2)
      term2 *= 4./(r**2)
    if (lensopt != 3):
      ans,err = quad(term3int,0.,pimax,args=(r,lrmod,ximod,lrminmod,dlrmod))
      term3 = 2.*ans
    if (lensopt == 1):
      lensmod[i] = 0.277518*om*(term1 - term2 - term3)
    elif (lensopt == 2):
      lensmod[i] = term3
    elif (lensopt == 3):
      lr = np.log(r)
      lensmod[i] = intxiarr(lr,lrmod,ximod,lrminmod,dlrmod)
    print r,lensmod[i]
  return rbin,lensmod

def term1int(lx,lrmod,ximod,lrminmod,dlrmod):
  x = np.exp(lx)
  xi = intxiarr(lx,lrmod,ximod,lrminmod,dlrmod)
  return (x**3)*xi

def term2int(lx,r,lrmod,ximod,lrminmod,dlrmod):
  x = np.exp(lx)
  if (x > r):
    xi = intxiarr(lx,lrmod,ximod,lrminmod,dlrmod)
    temp = (x**2)*xi*np.sqrt(x**2-r**2)
  else:
    temp = 0.
  return temp

def term3int(pi,rp,lrmod,ximod,lrminmod,dlrmod):
  lr = 0.5*np.log(pi**2 + rp**2)
  xi = intxiarr(lr,lrmod,ximod,lrminmod,dlrmod)
  return xi

########################################################################
# Generate spatial correlation function from power spectrum.           #
########################################################################

def getxiarr(rminmod,rmaxmod,nrmod,kminint,kmaxint,kmod,pkmod):
  rmod = np.logspace(np.log10(rminmod),np.log10(rmaxmod),nrmod)
  rdamp = 0.
  ximod = np.empty(nrmod)
  print '\nGenerating correlation function...'
  print 'rminmod =',rminmod,'rmaxmod =',rmaxmod,'nrmod =',nrmod
  print 'kminint =',kminint,'kmaxint =',kmaxint
  for i in range(nrmod):
    r = rmod[i]
    if (r < 50.):
      fracmin = 1.e-3
    elif (r < 100.):
      fracmin = 1.e-2
    elif (r < 500.):
      fracmin = 3.e-2
    elif (r < 1000.):
      fracmin = 5.e-2
    else:
      fracmin = 7.e-2
    ximod[i] = pktoxi(r,kminint,kmaxint,kmod,pkmod,rdamp,fracmin)
    print r,ximod[i]
  return rmod,ximod

def intxiarr(lr,lrarr,xiarr,lrmin,dlr):
  i = int((lr-lrmin)/dlr)
  fact = (lr-lrarr[i])/dlr
  xi = xiarr[i] + fact*(xiarr[i+1]-xiarr[i])
  return xi

########################################################################
# Convert power spectrum to correlation function.                      #
########################################################################

def pktoxi(s,kmin,kmax,kmod,pkmod,sdamp,fracmin):
  if (fracmin > 0.):
    epsabs,epsrel = 0.1*fracmin,0.1*fracmin
  else:
    epsabs,epsrel = 1.49e-8,1.49e-8
  lkmod = np.log(kmod)
  lkmin,dlk = lkmod[0],lkmod[1]-lkmod[0]
  lkstep = 1.
  xi = 0.
  istart = int(np.floor(kmin*s/np.pi)) + 1
  iend = int(np.floor(kmax*s/np.pi))
  kstart = istart*np.pi/s
  kend = iend*np.pi/s
  frac = 1.
  if (kstart > kmax):
    xi,err = quad(xiint,np.log(kmin),np.log(kmax),args=(s,lkmod,pkmod,lkmin,dlk,sdamp),epsabs=epsabs,epsrel=epsrel)
  else:
    nstep = int(round((np.log(kstart)-np.log(kmin))/lkstep))
    lkstep = (np.log(kstart)-np.log(kmin))/nstep
    for i in range(1,nstep+1):
      lk1 = np.log(kmin) + lkstep*(i-1)
      lk2 = np.log(kmin) + lkstep*i
      ans,err = quad(xiint,lk1,lk2,args=(s,lkmod,pkmod,lkmin,dlk,sdamp),epsabs=epsabs,epsrel=epsrel)
      xi += ans
      frac = ans/xi
      if (abs(frac) < fracmin):
        return xi/(2.*np.pi*np.pi*s)
    for i in range(istart,iend):
      k1 = i*np.pi/s
      k2 = (i+1)*np.pi/s
      ans,err = quad(xiint,np.log(k1),np.log(k2),args=(s,lkmod,pkmod,lkmin,dlk,sdamp),epsabs=epsabs,epsrel=epsrel)
      xi += ans
      frac = ans/xi
      if (abs(frac) < fracmin):
        return xi/(2.*np.pi*np.pi*s)
    if ((iend > 0) & (kend < kmax)):
      ans,err = quad(xiint,np.log(kend),np.log(kmax),args=(s,lkmod,pkmod,lkmin,dlk,sdamp),epsabs=epsabs,epsrel=epsrel)
      xi += ans
      frac = ans/xi
  return xi/(2.*np.pi*np.pi*s)

def xiint(lk,s,lkmod,pkmod,lkmin,dlk,sdamp):
  k = np.exp(lk)
  pk = np.interp(lk,lkmod,pkmod)
  return k*k*np.sin(k*s)*pk*np.exp(-((k*sdamp)**2))

########################################################################
# Pre-compute integrands for correlation functions.                    #
# lensopt: 1) xi_+ 2) xi_- 3) gamma_t 4) w(theta)                      #
########################################################################

def getclshapearr(lensopt,lmin,lmax,nl,zminlens,zmaxlens,zlensarr,flensarr,zminsource,zmaxsource,zsourcearr,fsourcearr,kmod,zmod,pkzmod,om,b,cosmo):
  if (lensopt == 1):
    print '\nComputing integrand for xi_+...'
  elif (lensopt == 2):
    print '\nComputing integrand for xi_-...'
  elif (lensopt == 3):
    print '\nComputing integrand for gamma_t...'
  elif (lensopt == 4):
    print '\nComputing integrand for w...'
  ch0 = 2997.9
  normg = b
  normk = (3.*om)/(2.*ch0*ch0)
  larr,clshapearr = np.logspace(np.log10(lmin),np.log10(lmax),nl),np.empty(nl)
  zmin0 = 0.005
  if ((lensopt == 1) or (lensopt == 2)):
    zmin,zmax = zmin0,zmaxsource
    norm = normk**2
  elif (lensopt == 3):
    zmin,zmax = max(zminlens,zmin0),min(zmaxlens,zmaxsource)
    norm = normk*normg
  elif (lensopt == 4):
    zmin,zmax = max(zminlens,zmin0),zmaxlens
    norm = normg**2
  print 'zmin =',zmin,'zmax =',zmax
  kmin,kmax = lmin/cosmo.comoving_distance(zmax).value,lmax/cosmo.comoving_distance(zmin).value
  if (kmin < kmod[0]):
    print 'Decrease value of kminmod!!'
    print kmin,kmod[0]
    sys.exit()
  if (kmax > kmod[-1]):
    print 'Increase value of kmaxmod!!'
    print kmax,kmod[-1]
    sys.exit()
  zarr = np.linspace(zmin,zmax,1000)
  rarr = cosmo.comoving_distance(zarr).value
  drdzarr = ch0*cosmo.inv_efunc(zarr)
  dz,dk = zmod[1]-zmod[0],np.log10(kmod[1])-np.log10(kmod[0])
  for i in range(nl):
    if ((lensopt == 1) or (lensopt == 2)):
      ans,err = quad(clkkint,zmin,zmax,args=(larr[i],zarr,rarr,drdzarr,zsourcearr,fsourcearr,kmod,dk,zmod,dz,pkzmod))
    elif (lensopt == 3):
      ans,err = quad(clgkint,zmin,zmax,args=(larr[i],zarr,rarr,drdzarr,zlensarr,flensarr,zsourcearr,fsourcearr,kmod,dk,zmod,dz,pkzmod))
    elif (lensopt == 4):
      ans,err = quad(clggint,zmin,zmax,args=(larr[i],zarr,rarr,drdzarr,zlensarr,flensarr,kmod,dk,zmod,dz,pkzmod))
    clshapearr[i] = norm*ans
    print larr[i],clshapearr[i]
  return larr,clshapearr

def clkkint(z,l,zarr,rarr,drdzarr,zsourcearr,fsourcearr,kmod,dk,zmod,dz,pkzmod):
  x = np.interp(z,zarr,rarr)
  dxdz = np.interp(z,zarr,drdzarr)
  fsource = np.interp(z,zsourcearr,fsourcearr)
  k = l/x
  pk = intpkzmod2(k,z,kmod,dk,zmod,dz,pkzmod)
  clkkint = dxdz*pk*(fsource**2)*((1.+z)**2)
  return clkkint

def clgkint(z,l,zarr,rarr,drdzarr,zlensarr,flensarr,zsourcearr,fsourcearr,kmod,dk,zmod,dz,pkzmod):
  x = np.interp(z,zarr,rarr)
  dxdz = np.interp(z,zarr,drdzarr)
  flens = np.interp(z,zlensarr,flensarr)
  fsource = np.interp(z,zsourcearr,fsourcearr)
  k = l/x
  pk = intpkzmod2(k,z,kmod,dk,zmod,dz,pkzmod)
  clgkint = dxdz*pk*flens*fsource*(1.+z)
  return clgkint

def clggint(z,l,zarr,rarr,drdzarr,zlensarr,flensarr,kmod,dk,zmod,dz,pkzmod):
  x = np.interp(z,zarr,rarr)
  dxdz = np.interp(z,zarr,drdzarr)
  flens = np.interp(z,zlensarr,flensarr)
  k = l/x
  pk = intpkzmod2(k,z,kmod,dk,zmod,dz,pkzmod)
  clggint = dxdz*pk*(flens**2)
  return clggint

########################################################################
# Compute angular lensing correlation functions.                       #
# lensopt: 1) xi_+ 2) xi_- 3) gamma_t 4) w(theta)                      #
########################################################################

def getlensangmod(lensopt,thmin,thmax,nbin,linlog,lmin,lmax,larr,clshapearr):
  if (lensopt == 1):
    print '\nGenerating model for xi_+...'
    lbess = 0
  elif (lensopt == 2):
    print '\nGenerating model for xi_-...'
    lbess = 4
  elif (lensopt == 3):
    print '\nGenerating model for gamma_t...'
    lbess = 2
  elif (lensopt == 4):
    print '\nGenerating model for w...'
    lbess = 0
  print 'lmin =',lmin,'lmax =',lmax
  if (linlog == 1):
    dth = (thmax-thmin)/float(nbin)
    thdeg = np.linspace(thmin+0.5*dth,thmax-0.5*dth,nbin)
  else:
    dth = (np.log10(thmax)-np.log10(thmin))/float(nbin)
    thdeg = np.logspace(np.log10(thmin)+0.5*dth,np.log10(thmax)-0.5*dth,nbin)
  thrad = np.radians(thdeg)
  nzero = 20000
  lensmod = np.zeros(nbin)
  for i in range(nbin):
    lzero = jn_zeros(lbess,nzero)/thrad[i]
    if (lzero[-1] < lmax):
      print lzero
      print lmax
      print 'Generate more zeros!!'
      sys.exit()
    lint = np.concatenate((np.array([lmin]),lzero[lzero < lmax],np.array([lmax])))
    for j in range(len(lint)-1):
      l1,l2 = lint[j],lint[j+1]
      dlensmod,err = quad(lensmodint,l1,l2,args=(thrad[i],lbess,larr,clshapearr))
      lensmod[i] += dlensmod/(2.*np.pi)
    print thdeg[i],lensmod[i]
  return thdeg,lensmod

def lensmodint(l,thrad,lbess,larr,clshapearr):
  cl = np.interp(l,larr,clshapearr)
  return l*cl*jn(lbess,l*thrad)

if __name__ == '__main__':
  main()
