************************************************************************
* Code to measure the galaxy-galaxy lensing statistic Delta Sigma(R).  *
* using a source shape catalogue behind a lens catalogue.              *
*                                                                      *
* Compile using:                                                       *
*                                                                      *
* gfortran -o delsigest delsigest.f                                    *
************************************************************************

      program delsigest
      implicit none

      integer maxsnumshape,maxsnumdens,maxnbin,maxnjack
      parameter (maxsnumshape=1000000)
      parameter (maxsnumdens=100000)
      parameter (maxnbin=30)
      parameter (maxnjack=200)
      integer snumshapej(maxnjack),snumdensj(maxnjack)
      integer numj(maxnjack)
      real srasdens(maxsnumdens),sdecdens(maxsnumdens)
      real sweidens(maxsnumdens),sjackdens(maxsnumdens)
      real sreddens(maxsnumdens)
      real srasshape(maxsnumshape),sdecshape(maxsnumshape)
      real sweishape(maxsnumshape),se1shape(maxsnumshape)
      real se2shape(maxsnumshape),sredshape(maxsnumshape)
      real sjackshape(maxsnumshape)
      real rminj(maxnjack),rmaxj(maxnjack)
      real dminj(maxnjack),dmaxj(maxnjack)
      real perparr(maxnbin),numbins(maxnbin),denbins(maxnbin)
      real corr(maxnbin),jcorr(maxnbin),jerror(maxnbin)
      real tjcorr(maxnbin),tjerror(maxnbin)
      real xjcorr(maxnbin),xjerror(maxnbin)
      real corrsave(maxnjack,maxnbin)
      double precision tnumjack(maxnbin,maxnjack,maxnjack)
      double precision xnumjack(maxnbin,maxnjack,maxnjack)
      double precision denjack(maxnbin,maxnjack,maxnjack)
      integer nbin,linlog,nxjack,nyjack,ii,jj,kk,njack,snumshape
      integer snumdens,ijack,icorr
      real minperp,maxperp,dperp,omfid,rmin,rmax,dmin,dmax
      real delz,mu,sig,sumrweidens,zminl,zmaxl,zmins,zmaxs
      real rminshape,rmaxshape,dminshape,dmaxshape,rmindens,rmaxdens
      real dmindens,dmaxdens
      double precision dsum,dsumsq,dtemp,numsum,densum
      character*100 densfile,shapefile
      common /cosmo/ omfid

C----------------------------------------------------------------------C
C Input data files for sources and lenses.                             C
C----------------------------------------------------------------------C

      shapefile = 'sources_buzzard_pix128_zp0pt80_1pt00.dat'
      densfile = 'lenses_buzzard_lrg_pix128_zs0pt40_0pt60.dat'

C----------------------------------------------------------------------C
C Parameters.                                                          C
C----------------------------------------------------------------------C

      zmins = 0.8     ! Minimum source redshift
      zmaxs = 1.0     ! Maximum source redshift
      zminl = 0.4     ! Minimum lens redshift
      zmaxl = 0.6     ! Maximum lens redshift
      delz = 0.       ! Make measurement for z_s > z_l + delz
      minperp = 0.5   ! Minimum transverse separation in Mpc/h
      maxperp = 50.   ! Maximum transverse separation in Mpc/h
      nbin = 10       ! Number of transverse bins
      linlog = 2      ! (1) linear (2) logarithmic binning
      nxjack = 7      ! Number of jack-knife regions in x-direction
      nyjack = 7      ! Number of jack-knife regions in y-direction
      omfid = 0.3     ! Omega_m for fiducial cosmology

C----------------------------------------------------------------------C
C Initializations.                                                     C
C----------------------------------------------------------------------C

      njack = nxjack*nyjack
      dperp = (maxperp-minperp)/real(nbin)
      do ii=1,nbin
        if (linlog.eq.1) then
          perparr(ii) = minperp + dperp*(real(ii)-0.5)
        else
          perparr(ii) = minperp*
     >     exp(((real(ii)-0.5)/real(nbin))*log(maxperp/minperp))
        endif
      enddo
      if (nbin.gt.maxnbin) then
        write(*,'("Increase value of maxnbin!!")')
        stop
      endif
      if (njack.gt.maxnjack) then
        write(*,'("Decrease value of njack!!")')
        stop
      endif
      write(*,'(" ")')
      write(*,'("Measuring galaxy-galaxy lensing...")')
      write(*,'("zminl = ",F4.2)')zminl
      write(*,'("zmaxl = ",F4.2)')zmaxl
      write(*,'("zmins = ",F4.2)')zmins
      write(*,'("zmaxs = ",F4.2)')zmaxs

C----------------------------------------------------------------------C
C Read in shape sample.                                                C
C----------------------------------------------------------------------C

      call readshapesim(shapefile,zmins,zmaxs,maxsnumshape,srasshape,
     >                  sdecshape,sredshape,sweishape,se1shape,se2shape,
     >                  snumshape,rminshape,rmaxshape,dminshape,
     >                  dmaxshape)

C----------------------------------------------------------------------C
C Read in density sample.                                              C
C----------------------------------------------------------------------C

      call readdenssim(densfile,zminl,zmaxl,maxsnumdens,srasdens,
     >                 sdecdens,sreddens,sweidens,snumdens,rmindens,
     >                 rmaxdens,dmindens,dmaxdens)

C----------------------------------------------------------------------C
C Determine jack-knife partitions for data samples.                    C
C----------------------------------------------------------------------C

      rmin = min(rminshape,rmindens)
      rmax = max(rmaxshape,rmaxdens)
      dmin = min(dminshape,dmindens)
      dmax = max(dmaxshape,dmaxdens)
      call sort6col(sdecshape,srasshape,sredshape,se1shape,se2shape,
     >              sweishape,snumshape)
      call getjackreg(srasshape,sdecshape,snumshape,nxjack,nyjack,
     >                 rmin,rmax,dmin,dmax,rminj,rmaxj,dminj,dmaxj)
      call ptstojack(srasshape,sdecshape,snumshape,nxjack,nyjack,
     >               rminj,rmaxj,dminj,dmaxj,sjackshape,snumshapej)
      call sort4col(sdecdens,srasdens,sreddens,sweidens,snumdens)
      call ptstojack(srasdens,sdecdens,snumdens,nxjack,nyjack,rminj,
     >               rmaxj,dminj,dmaxj,sjackdens,snumdensj)

C----------------------------------------------------------------------C
C Message to screen.                                                   C
C----------------------------------------------------------------------C

      write(*,'(" ")')
      write(*,'("ndens   = ",I8)')snumdens
      write(*,'("nshape  = ",I8)')snumshape
      write(*,'("Sample boundaries...")')
      write(*,'(F7.3," < R.A. < ",F7.3)')rmin,rmax
      write(*,'(F7.3," < Dec. < ",F7.3)')dmin,dmax
      write(*,'(" ")')
      write(*,'("Measuring projected correlation function...")')
      write(*,'("minperp = ",F5.1)')minperp
      write(*,'("maxperp = ",F5.1)')maxperp
      write(*,'("nbin    = ",I2)')nbin
      write(*,'("linlog  = ",I1)')linlog
      write(*,'("delz    = ",F4.1)')delz

C----------------------------------------------------------------------C
C Measure the cross-correlation function of the tangential and cross   C
C component of shapes, with the position of the lenses, divided into   C
C jack-knife regions.                                                  C
C----------------------------------------------------------------------C

      write(*,'(" ")')
      write(*,'("Generating gtD and gxD...")')
      call binshape(srasdens,sdecdens,sreddens,sweidens,sjackdens,
     >              snumdens,srasshape,sdecshape,sredshape,sweishape,
     >              se1shape,se2shape,sjackshape,snumshape,minperp,
     >              maxperp,linlog,nbin,tnumjack,xnumjack,denjack,delz)

C----------------------------------------------------------------------C
C Construct the correlation functions and errors from the measurements C
C in jack-knife regions.  icorr=1 uses the tangential component of the C
C shape, and icorr=2 uses the cross component of the shape.            C
C----------------------------------------------------------------------C

      do icorr=1,2

C----------------------------------------------------------------------C
C The total correlation function in each bin.                          C
C----------------------------------------------------------------------C

        do ii=1,nbin
          numsum = 0D0
          densum = 0D0
          do jj=1,njack
            do kk=1,njack
              if (icorr.eq.1) then
                numsum = numsum + tnumjack(ii,jj,kk)
                densum = densum + denjack(ii,jj,kk)
              else if (icorr.eq.2) then
                numsum = numsum + xnumjack(ii,jj,kk)
                densum = densum + denjack(ii,jj,kk)
              endif
            enddo
          enddo
          numbins(ii) = real(numsum)
          denbins(ii) = real(densum)
        enddo
        do ii=1,nbin
          corr(ii) = -1.
          if (denbins(ii).gt.0.) then
            corr(ii) = numbins(ii)/denbins(ii)
          endif
        enddo

C----------------------------------------------------------------------C
C The correlation function deleting each jack-knife region.            C
C----------------------------------------------------------------------C

        do ijack=1,njack
          do ii=1,nbin
            numsum = 0D0
            densum = 0D0
            do jj=1,njack
              do kk=1,njack
                if ((jj.ne.ijack).and.(kk.ne.ijack)) then
                  if (icorr.eq.1) then
                    numsum = numsum + tnumjack(ii,jj,kk)
                    densum = densum + denjack(ii,jj,kk)
                  else if (icorr.eq.2) then
                    numsum = numsum + xnumjack(ii,jj,kk)
                    densum = densum + denjack(ii,jj,kk)
                  endif
                endif
              enddo
            enddo
            numbins(ii) = real(numsum)
            denbins(ii) = real(densum)
          enddo
          do ii=1,nbin
            jcorr(ii) = -1.
            if (denbins(ii).gt.0.) then
              jcorr(ii) = numbins(ii)/denbins(ii)
            endif
          enddo
          do ii=1,nbin
            corrsave(ijack,ii) = jcorr(ii)
          enddo
        enddo

C----------------------------------------------------------------------C
C Use these correlation functions to determine the jack-knife error.   C
C----------------------------------------------------------------------C
       
        do ii=1,nbin
          dsum = 0D0
          dsumsq = 0D0
          do ijack=1,njack
            dtemp = dble(corrsave(ijack,ii))
            dsum = dsum + dtemp
            dsumsq = dsumsq + dtemp**2
          enddo
          mu = real(dsum)/real(njack)
          sig = sqrt((real(dsumsq)/real(njack))-(mu**2))
          jcorr(ii) = mu
          jerror(ii) = sig*sqrt(real(njack-1))
        enddo
        do ii=1,nbin
          if (icorr.eq.1) then
            tjcorr(ii) = corr(ii)
            tjerror(ii) = jerror(ii)
          else if (icorr.eq.2) then
            xjcorr(ii) = corr(ii)
            xjerror(ii) = jerror(ii)
          endif
        enddo
      enddo

C----------------------------------------------------------------------C
C Display correlation function measurements.                           C
C----------------------------------------------------------------------C

      write(*,'(" ")')
      write(*,'("Delta Sigma measurements:")')
      write(*,'("----------------------------------------")')
      write(*,'("       R    DS_t     err    DS_x     err")')
      write(*,'("----------------------------------------")')
      do ii=1,nbin
        write(*,'(5F8.4)')perparr(ii),tjcorr(ii),tjerror(ii),xjcorr(ii),
     >   xjerror(ii)
      enddo
      write(*,'("----------------------------------------")')

C----------------------------------------------------------------------C
C Output data to file.                                                 C
C----------------------------------------------------------------------C

      write(*,'(" ")')
      write(*,'("Outputting results to file...")')
      open(1,file='delsigest.dat',status='new')
      do ii=1,nbin
        write(1,*)perparr(ii),tjcorr(ii),tjerror(ii),xjcorr(ii),
     >   xjerror(ii)
      enddo
      close(1)

      end

************************************************************************
* Read in shape catalogue.                                             *
************************************************************************

      subroutine readshapesim(infile,zmins,zmaxs,maxsnumshape,srasshape,
     >                        sdecshape,sredshape,sweishape,se1shape,
     >                        se2shape,snumshape,rmin,rmax,dmin,dmax)
      implicit none

      real srasshape(*),sdecshape(*),sredshape(*),sweishape(*)
      real se1shape(*),se2shape(*)
      integer maxsnumshape,snumshape,ii,jj,ioerr
      real zmins,zmaxs,rmin,rmax,dmin,dmax,rr,dd,zz,e1,e2
      character*100 infile,line

      write(*,'(" ")')
      write(*,'("Reading in simulated shape catalogue...")')
      write(*,*)infile
      open(1,file=infile,status='old')
      do ii=1,2
        read(1,'(a)')line
      enddo
      read(1,*)snumshape
      write(*,'("nshape = ",I8)')snumshape
      ioerr = 0
      ii = 0
      jj = 0
      do while (ioerr.eq.0)
        read(1,'(a)',iostat=ioerr)line
        if (ioerr.eq.0) then
          read(line,*)rr,dd,zz,e1,e2
          jj = jj + 1
          if ((zz.gt.zmins).and.(zz.lt.zmaxs)) then
            ii = ii + 1
            srasshape(ii) = rr
            sdecshape(ii) = dd
            sredshape(ii) = zz
            se1shape(ii) = e1
            se2shape(ii) = e2
            sweishape(ii) = 1.
          endif
        endif
      enddo
      close(1)
      snumshape = ii
      write(*,'("snumshape = ",I8," of ",I8," for ",F5.2," < z < ",
     > F4.2)')snumshape,jj,zmins,zmaxs
      if (snumshape.gt.maxsnumshape) then
        write(*,'("readshapesim - increase value of maxsnumshape!!")')
        stop
      endif
      rmin = 999999.
      dmin = 999999.
      rmax = -999999.
      dmax = -999999.
      do ii=1,snumshape
        rr = srasshape(ii)
        dd = sdecshape(ii)
        if (rr.lt.rmin) rmin = rr
        if (rr.gt.rmax) rmax = rr
        if (dd.lt.dmin) dmin = dd
        if (dd.gt.dmax) dmax = dd
      enddo
      write(*,'("Boundaries of shape sample:")')
      write(*,'(F7.3," < R.A. < ",F7.3)')rmin,rmax
      write(*,'(F7.3," < Dec. < ",F7.3)')dmin,dmax

      end

************************************************************************
* Read in lens catalogue.                                              *
************************************************************************

      subroutine readdenssim(densfile,zminl,zmaxl,maxsnumdens,srasdens,
     >                       sdecdens,sreddens,sweidens,snumdens,rmin,
     >                       rmax,dmin,dmax)
      implicit none

      real srasdens(*),sdecdens(*),sreddens(*),sweidens(*)
      integer seed,maxsnumdens,snumdens,ioerr,ii,jj
      real zminl,zmaxl,rmin,rmax,dmin,dmax,rr,dd,zz
      character*100 densfile,line

      write(*,'(" ")')
      write(*,'("Reading in simulated density data catalogue...")')
      write(*,'("densfile = ",a60)')densfile
      open(1,file=densfile,status='old')
      do ii=1,2
        read(1,'(a)')line
      enddo
      read(1,*)snumdens
      ioerr = 0
      ii = 0
      jj = 0
      do while (ioerr.eq.0)
        read(1,'(a)',iostat=ioerr)line
        if (ioerr.eq.0) then
          read(line,*)rr,dd,zz
          jj = jj + 1
          if ((zz.gt.zminl).and.(zz.lt.zmaxl)) then
            ii = ii + 1
            srasdens(ii) = rr
            sdecdens(ii) = dd
            sreddens(ii) = zz
            sweidens(ii) = 1.
          endif
        endif
      enddo
      close(1)
      snumdens = ii
      write(*,'("snumdens = ",I8," of ",I8," for ",F5.2," < z < ",
     > F4.2)')snumdens,jj,zminl,zmaxl
      if (snumdens.gt.maxsnumdens) then
        write(*,'("readdenssim - increase value of maxsnumdens!!")')
        stop
      endif
      rmin = 999999.
      dmin = 999999.
      rmax = -999999.
      dmax = -999999.
      do ii=1,snumdens
        rr = srasdens(ii)
        dd = sdecdens(ii)
        if (rr.lt.rmin) rmin = rr
        if (rr.gt.rmax) rmax = rr
        if (dd.lt.dmin) dmin = dd
        if (dd.gt.dmax) dmax = dd
      enddo
      write(*,'(" ")')
      write(*,'("Boundaries of density sample:")')
      write(*,'(F7.3," < R.A. < ",F7.3)')rmin,rmax
      write(*,'(F7.3," < Dec. < ",F7.3)')dmin,dmax

      end

************************************************************************
* Define jack-knife regions as (ra,dec) boundaries containing equal    *
* number of sources.                                                   *
************************************************************************

      subroutine getjackreg(sras,sdec,snum,nxjack,nyjack,rmin,rmax,
     >                      dmin,dmax,rminj,rmaxj,dminj,dmaxj)
      implicit none

      integer maxsnum1
      parameter (maxsnum1=1000000)
      real sras1(maxsnum1),sdec1(maxsnum1)
      real sras(*),sdec(*),rminj(*),rmaxj(*),dminj(*),dmaxj(*)
      integer nxjack,nyjack,snum,ijack,idyjack,iyjack,ii,snum1,idxjack
      integer ixjack,snum2
      real rmin,rmax,dmin,dmax,r1,r2,d1,d2,xmin,xmax,ymin,ymax

      write(*,'(" ")')
      write(*,'("Defining equal-number 2D jack-knife regions...")')
      do ii=2,snum
        if (sdec(ii).lt.sdec(ii-1)) then
          write(*,*)ii,sdec(ii-1),sdec(ii)
          write(*,
     > '("getjackreg - sort catalogue in ascending declination!!")')
          stop
        endif
      enddo
      ijack = 0
      idyjack = nint(real(snum)/real(nyjack))
      do iyjack=1,nyjack
        if (nyjack.eq.1) then
          d1 = dmin
          d2 = dmax
        else if (iyjack.eq.1) then
          ii = idyjack
          d1 = dmin
          d2 = 0.5*(sdec(ii)+sdec(ii+1))
        else if (iyjack.eq.nyjack) then
          ii = idyjack*(nyjack-1)
          d1 = 0.5*(sdec(ii)+sdec(ii+1))
          d2 = dmax
        else
          ii = idyjack*(iyjack-1)
          d1 = 0.5*(sdec(ii)+sdec(ii+1))
          ii = idyjack*iyjack
          d2 = 0.5*(sdec(ii)+sdec(ii+1))
        endif
        snum1 = 0
        do ii=1,snum
          if ((sdec(ii).gt.d1).and.(sdec(ii).lt.d2)) then
            snum1 = snum1 + 1
            sras1(snum1) = sras(ii)
            sdec1(snum1) = sdec(ii)
          endif
        enddo
        if (snum1.gt.maxsnum1) then
          write(*,'("snum1 = ",I7)')snum1
          write(*,'("getjackreg - increase value of maxsnum1!!")')
          stop
        endif
        call sort2col(sras1,sdec1,snum1)
        idxjack = nint(real(snum1)/real(nxjack))
        do ixjack=1,nxjack
          if (nxjack.eq.1) then
            r1 = rmin
            r2 = rmax
          else if (ixjack.eq.1) then
            ii = idxjack
            r1 = rmin
            r2 = 0.5*(sras1(ii)+sras1(ii+1))
          else if (ixjack.eq.nxjack) then
            ii = idxjack*(nxjack-1)
            r1 = 0.5*(sras1(ii)+sras1(ii+1))
            r2 = rmax
          else
            ii = idxjack*(ixjack-1)
            r1 = 0.5*(sras1(ii)+sras1(ii+1))
            ii = idxjack*ixjack
            r2 = 0.5*(sras1(ii)+sras1(ii+1))
          endif
          snum2 = 0
          do ii=1,snum
            if ((sras(ii).gt.r1).and.(sras(ii).lt.r2).and.
     >          (sdec(ii).gt.d1).and.(sdec(ii).lt.d2)) then
              snum2 = snum2 + 1
            endif
          enddo
          ijack = ijack + 1
          rminj(ijack) = r1
          rmaxj(ijack) = r2
          dminj(ijack) = d1
          dmaxj(ijack) = d2
        enddo
      enddo

      end

************************************************************************
* Assign jack-knife regions to the data points using equal-number      *
* jack-knife binning.                                                  *
************************************************************************

      subroutine ptstojack(sras,sdec,snum,nxjack,nyjack,rminj,rmaxj,
     >                     dminj,dmaxj,sjack,numj)
      implicit none

      integer numj(*)
      real sras(*),sdec(*),sjack(*),rminj(*),rmaxj(*),dminj(*),dmaxj(*)
      integer nxjack,nyjack,snum,njack,ii,jj,kk
      real d1

      do ii=2,snum
        if (sdec(ii).lt.sdec(ii-1)) then
          write(*,'("ptstojack - sort in ascending declination!!")')
          write(*,*)ii,sdec(ii-1),sdec(ii)
          stop
        endif
      enddo
      njack = nxjack*nyjack
      do ii=1,njack
        numj(ii) = 0
      enddo
      d1 = dmaxj(1)
      jj = 1
      do ii=1,snum
        if (sdec(ii).gt.d1) then
          jj = jj + 1
          kk = (nxjack*(jj-1)) + 1
          d1 = dmaxj(kk)
        endif
        kk = (nxjack*(jj-1)) + 1
        do while ((sras(ii).lt.rminj(kk)).or.(sras(ii).gt.rmaxj(kk)))
          kk = kk + 1
        enddo
        sjack(ii) = real(kk)
        numj(kk) = numj(kk) + 1
      enddo

      end

************************************************************************
* Perform the bin counts in jack-knife regions required to estimate    *
* the galaxy-galaxy lensing statistic Delta Sigma(R) for the           *
* tangential and cross component of the shape around each lens.        *
* Sources are only binned in their redshift z_s > z_l + delz.          *
************************************************************************

      subroutine binshape(rasdens,decdens,reddens,weidens,jackdens,
     >                    ndens,rasshape,decshape,redshape,weishape,
     >                    e1shape,e2shape,jackshape,nshape,minperp,
     >                    maxperp,linlog,nperp,tnumjack,xnumjack,
     >                    denjack,delz)
      implicit none

      integer maxnn,maxnbin,maxnjack
      parameter (maxnn=1000000)
      parameter (maxnbin=30)
      parameter (maxnjack=200)
      double precision tnumjack(maxnbin,maxnjack,maxnjack)
      double precision xnumjack(maxnbin,maxnjack,maxnjack)
      double precision denjack(maxnbin,maxnjack,maxnjack)
      real distdens(maxnn),distshape(maxnn)
      real rasdens(*),decdens(*),reddens(*),weidens(*),jackdens(*)
      real rasshape(*),decshape(*),redshape(*),weishape(*),e1shape(*)
      real e2shape(*),jackshape(*)
      integer ndens,nshape,linlog,nperp,ii,jj,kk,j1,j2,ibin,incr
      real minperp,maxperp,rl,dl,xl,wl,rs,ds,xs,ws,e1,e2,zz,minsq,maxsq
      real dperp,lminperp,lmaxperp,pi,fact,rfact,dfact,cosd,perpsq,xfact
      real lperp,phi,et,ec,maxsep,perp,sigc,wls,ch0,sigcnorm,zl,zs,wtot
      real temp,delz
      real getxx
      external getxx

      incr = 1000
      if (max(ndens,nshape).gt.maxnn) then
        write(*,'("binshape - increase value of maxnn!!")')
        stop
      endif
      if (nperp.gt.maxnbin) then
        write(*,'("binshape - increase value of maxnbin!!")')
        stop
      endif
      do ii=2,ndens
        if (decdens(ii).lt.decdens(ii-1)) then
          write(*,'("Bad declination ordering!")')
          stop
        endif
      enddo
      do ii=2,nshape
        if (decshape(ii).lt.decshape(ii-1)) then
          write(*,'("Bad declination ordering!")')
          stop
        endif
      enddo
      minsq = minperp*minperp
      maxsq = maxperp*maxperp
      if (linlog.eq.1) then
        dperp = (maxperp-minperp)/real(nperp)
      else
        lminperp = log(minperp)
        lmaxperp = log(maxperp)
        dperp = (lmaxperp-lminperp)/real(nperp)
      endif
      do ii=1,nperp
        do jj=1,maxnjack
          do kk=1,maxnjack
            tnumjack(ii,jj,kk) = 0D0
            xnumjack(ii,jj,kk) = 0D0
            denjack(ii,jj,kk) = 0D0
          enddo
        enddo
      enddo
      pi = 4.*atan(1.)
      fact = pi/180.
      ch0 = 2997.9
      sigcnorm = 554.65 ! units h M_sol pc^-2
      do ii=1,ndens
        zz = reddens(ii)
        distdens(ii) = getxx(zz)
      enddo
      do ii=1,nshape
        zz = redshape(ii)
        distshape(ii) = getxx(zz)
      enddo
      do ii=1,ndens
        if (mod(ii,incr).eq.0) then
          write(*,'("Counted ",I6,"/",I6," objects")')ii,ndens
        endif
        rl = rasdens(ii)
        dl = decdens(ii)
        xl = distdens(ii)
        zl = reddens(ii)
        wl = weidens(ii)
        j1 = nint(jackdens(ii))
        xfact = xl*fact
        maxsep = maxperp/xfact
        do jj=1,nshape
          ds = decshape(jj)
          dfact = ds-dl
          if ((-dfact).le.maxsep) then
            if (dfact.gt.maxsep) goto 200
            xs = distshape(jj)
            zs = redshape(jj)
            if ((zs.gt.(zl+delz)).and.(xs.gt.xl)) then
              dfact = abs(dfact)
              rs = rasshape(jj)
              rfact = abs(rl-rs)
              if (rfact.gt.180.) rfact = 360.-rfact
              cosd = cos(0.5*fact*(dl+ds))
              rfact = cosd*rfact
              perpsq = (xfact**2)*(rfact**2 + dfact**2)
              if ((perpsq.gt.minsq).and.(perpsq.lt.maxsq)) then
                ws = weishape(jj)
                j2 = nint(jackshape(jj))
C Sigma_c = (2/3) D_H rho_c sigc = (554.65 h M_sol pc^-2) sigc
                sigc = (ch0*xs)/((1.+zl)*xl*(xs-xl))
                wls = (1./sigc)**2
                if (linlog.eq.1) then
                  perp = sqrt(perpsq)
                  ibin = int((perp-minperp)/dperp) + 1
                else
                  lperp = 0.5*log(perpsq)
                  ibin = int((lperp-lminperp)/dperp) + 1                
                endif
                wtot = wl*ws*wls
                temp = sigcnorm*sigc
                e1 = e1shape(jj)
                e2 = e2shape(jj)
                call septoangle(rfact,dfact,rl,rs,dl,ds,phi)
                call angletoellip(phi,e1,e2,et,ec)
                tnumjack(ibin,j1,j2) = tnumjack(ibin,j1,j2) +
     >           dble(wtot*et*temp)
                xnumjack(ibin,j1,j2) = xnumjack(ibin,j1,j2) +
     >           dble(wtot*ec*temp)
                denjack(ibin,j1,j2) = denjack(ibin,j1,j2) + dble(wtot)
              endif
            endif
          endif
        enddo
200   enddo

      end

************************************************************************
* Rotate galaxy ellipticity.                                           *
* phi = angle of separation line to R.A.=0 [-90 deg < phi < 90 deg]    *
************************************************************************

      subroutine angletoellip(phi,e1,e2,et,ec)
      implicit none

      real phi,e1,e2,et,ec,twophi

      twophi = 2.*phi
      et = - e1*cos(twophi) - e2*sin(twophi)
      ec = e1*sin(twophi) - e2*cos(twophi)

      end

************************************************************************
* Angle of separation line between (r1,d1) and (r2,d2) to positive R.A.*
* phi is constrained to lie in -90 deg < phi < 90 deg.                 *
************************************************************************

      subroutine septoangle(rfact,dfact,r1,r2,d1,d2,phi)
      implicit none

      real rfact,dfact,r1,r2,d1,d2,phi,temp

C rfact and dfact are always positive, returns value in range 0 to pi/2
      temp = atan(dfact/rfact)
      if (d2.gt.d1) then
        if (r2.gt.r1) then
          phi = temp
        else
          phi = -temp
        endif
      else
        if (r2.gt.r1) then
          phi = -temp
        else
          phi = temp
        endif
      endif

      end

************************************************************************
* Co-moving distance for a given redshift.                             *
************************************************************************

      function getxx(zz)
      implicit none

      real getxx,zz,getdxdz,ans
      external getdxdz

      call qsimp(getdxdz,0.,zz,ans)
      getxx = ans

      end

      function getdxdz(zz)
      implicit none

      real getdxdz,zz,ch0,om
      common /cosmo/ om

      ch0 = 2997.9
      getdxdz = ch0/sqrt((om*((1.+zz)**3)) + (1.-om))

      end

**********************************************************************
* Sort a survey of 2 columns in order of increasing value of the     *
* first column.                                                      *
**********************************************************************

      subroutine sort2col(col1,col2,num)
      implicit none

      integer maxnum
      parameter (maxnum=1000000)
      real col1(*),col2(*)
      integer indexord(maxnum)
      real temparray(maxnum)
      integer num,ii,jj,kk,ll,mm
      real temp

      if (num.gt.maxnum) then
        write(*,'("sort2col - Too many points!!")')
        write(*,*)num
        stop
      endif
      do ii=1,num
        indexord(ii)=ii
      enddo
      jj=num/2+1
      kk=num
300   continue
      if (jj.gt.1) then
        jj=jj-1
        ll=indexord(jj)
        temp=col1(ll)
      else
        ll=indexord(kk)
        temp=col1(ll)
        indexord(kk)=indexord(1)
        kk=kk-1
        if (kk.eq.1) then
          indexord(1)=ll
          go to 500
        endif
      endif
      mm=jj
      ii=jj+jj
400   if (ii.le.kk) then
        if (ii.lt.kk) then
          if (col1(indexord(ii)).lt.col1(indexord(ii+1))) then
            ii=ii+1
          endif
        endif  
        if (temp.lt.col1(indexord(ii))) then
          indexord(mm)=indexord(ii)
          mm=ii
          ii=ii+ii
        else
          ii=kk+1
        endif
        go to 400
      endif
      indexord(mm)=ll
      go to 300
500   do ii=1,num
        temparray(ii) = col1(ii)
      enddo
      do ii=1,num
        jj = indexord(ii)
        col1(ii) = temparray(jj)
      enddo
      do ii=1,num
        temparray(ii) = col2(ii)
      enddo
      do ii=1,num
        jj = indexord(ii)
        col2(ii) = temparray(jj)
      enddo

      end

**********************************************************************
* Sort a survey of 4 columns in order of increasing value of the     *
* first column.                                                      *
**********************************************************************

      subroutine sort4col(col1,col2,col3,col4,num)
      implicit none

      integer maxnum
      parameter (maxnum=1000000)
      integer indexord(maxnum)
      real realarray(maxnum),col1(*),col2(*),col3(*),col4(*)
      integer num,ii,jj,kk,ll,mm
      real temp

      if (num.gt.maxnum) then
        write(*,'("sort4col - Too many points!!")')
        write(*,*)num
        stop
      endif
      do ii=1,num
        indexord(ii)=ii
      enddo
      jj=num/2+1
      kk=num
300   continue
      if (jj.gt.1) then
        jj=jj-1
        ll=indexord(jj)
        temp=col1(ll)
      else
        ll=indexord(kk)
        temp=col1(ll)
        indexord(kk)=indexord(1)
        kk=kk-1
        if (kk.eq.1) then
          indexord(1)=ll
          go to 500
        endif
      endif
      mm=jj
      ii=jj+jj
400   if (ii.le.kk) then
        if (ii.lt.kk) then
          if (col1(indexord(ii)).lt.col1(indexord(ii+1))) then
            ii=ii+1
          endif
        endif  
        if (temp.lt.col1(indexord(ii))) then
          indexord(mm)=indexord(ii)
          mm=ii
          ii=ii+ii
        else
          ii=kk+1
        endif
        go to 400
      endif
      indexord(mm)=ll
      go to 300
500   do ii=1,num
        realarray(ii) = col1(ii)
      enddo
      do ii=1,num
        jj = indexord(ii)
        col1(ii) = realarray(jj)
      enddo
      do ii=1,num
        realarray(ii) = col2(ii)
      enddo
      do ii=1,num
        jj = indexord(ii)
        col2(ii) = realarray(jj)
      enddo
      do ii=1,num
        realarray(ii) = col3(ii)
      enddo
      do ii=1,num
        jj = indexord(ii)
        col3(ii) = realarray(jj)
      enddo
      do ii=1,num
        realarray(ii) = col4(ii)
      enddo
      do ii=1,num
        jj = indexord(ii)
        col4(ii) = realarray(jj)
      enddo

      end

**********************************************************************
* Sort a survey of 6 columns in order of increasing value of the     *
* first column.                                                      *
**********************************************************************

      subroutine sort6col(col1,col2,col3,col4,col5,col6,num)
      implicit none

      integer maxnum
      parameter (maxnum=1000000)
      integer indexord(maxnum)
      real realarray(maxnum),col1(*),col2(*),col3(*),col4(*),col5(*)
      real col6(*)
      integer num,ii,jj,kk,ll,mm
      real temp

      if (num.gt.maxnum) then
        write(*,'("sort6col - Too many points!!")')
        write(*,*)num
        stop
      endif
      do ii=1,num
        indexord(ii)=ii
      enddo
      jj=num/2+1
      kk=num
300   continue
      if (jj.gt.1) then
        jj=jj-1
        ll=indexord(jj)
        temp=col1(ll)
      else
        ll=indexord(kk)
        temp=col1(ll)
        indexord(kk)=indexord(1)
        kk=kk-1
        if (kk.eq.1) then
          indexord(1)=ll
          go to 500
        endif
      endif
      mm=jj
      ii=jj+jj
400   if (ii.le.kk) then
        if (ii.lt.kk) then
          if (col1(indexord(ii)).lt.col1(indexord(ii+1))) then
            ii=ii+1
          endif
        endif  
        if (temp.lt.col1(indexord(ii))) then
          indexord(mm)=indexord(ii)
          mm=ii
          ii=ii+ii
        else
          ii=kk+1
        endif
        go to 400
      endif
      indexord(mm)=ll
      go to 300
500   do ii=1,num
        realarray(ii) = col1(ii)
      enddo
      do ii=1,num
        jj = indexord(ii)
        col1(ii) = realarray(jj)
      enddo
      do ii=1,num
        realarray(ii) = col2(ii)
      enddo
      do ii=1,num
        jj = indexord(ii)
        col2(ii) = realarray(jj)
      enddo
      do ii=1,num
        realarray(ii) = col3(ii)
      enddo
      do ii=1,num
        jj = indexord(ii)
        col3(ii) = realarray(jj)
      enddo
      do ii=1,num
        realarray(ii) = col4(ii)
      enddo
      do ii=1,num
        jj = indexord(ii)
        col4(ii) = realarray(jj)
      enddo
      do ii=1,num
        realarray(ii) = col5(ii)
      enddo
      do ii=1,num
        jj = indexord(ii)
        col5(ii) = realarray(jj)
      enddo
      do ii=1,num
        realarray(ii) = col6(ii)
      enddo
      do ii=1,num
        jj = indexord(ii)
        col6(ii) = realarray(jj)
      enddo

      end

**********************************************************************
* Numerical Recipes integration routines.                            *
**********************************************************************

      SUBROUTINE qsimp(func,a,b,s)
      INTEGER JMAX
      REAL a,b,func,s,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20)
CU    USES trapzd
      INTEGER j
      REAL os,ost,st
      ost=-1.e30
      os= -1.e30
      do 11 j=1,JMAX
        call trapzd(func,a,b,st,j)
        s=(4.*st-ost)/3.
        if (j.gt.5) then
          if (abs(s-os).lt.EPS*abs(os).or.(s.eq.0..and.os.eq.0.)) return
        endif
        os=s
        ost=st
11    continue
      write(*,*)a,b
      write(*,'("too many steps in qsimp")')
      stop
      END

      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
