! Data block with parameter structure
      block data text1r_data_block
        implicit none
        include 'text1r.fh'
      end block data

! Initialize parameter structure
      subroutine text1r_init(ttc,p,nu0,r, itype)
        implicit none
        real*8 ttc,p,nu0,r
        integer itype, i, zz,hh
        include 'text1r.fh'
        include 'he3.f90h'

        text_a = he3_text_a(ttc,p)
        text_lg2 = he3_text_lg2(ttc,p)
        text_lg1 = text_lg2 * (2 + he3_text_delta(ttc,p))
        text_lhv = he3_text_lhv(ttc,p)
        text_lsg = 3D0
        text_ld  = he3_ld(ttc,p)
        text_lo  = 0D0
        text_d = he3_text_d(ttc,p)
        text_r = r
        text_h = const_2pi*nu0/he3_gyro

        ! initial conditions for alpha_n, beta_n
        do i=1,MAXN
          hh = itype/2
          zz = (itype - hh*2)*2-1 ! parity of itype -1 or 1
          text_an(i) = -zz*acos(0.5D0)
          text_bn(i) = 2*acos(0D0)*hh + acos(-zz*1D0/sqrt(5D0)) * &
                       dble(i-1)/dble(MAXN-1)
          text_rr(i)=text_r*dble(i-1)/dble(MAXN-1)
          text_apsi(i)=0D0
          text_vr(i)=0D0
          text_vz(i)=0D0
          text_vf(i)=0D0
          text_lr(i)=0D0
          text_lz(i)=0D0
          text_lf(i)=0D0
          text_w(i)=0D0
        enddo
      end

! Vortex cluster profile
      subroutine text1r_set_vortex_cluster(omega, omega_v)
        implicit none
        include 'text1r.fh'
        integer i
        real*8 omega,omega_v,rr
        do i=1,MAXN
          rr = dble(i-1)/dble(MAXN-1)*text_r
          ! flow velocity
          text_vr(i)=0D0
          text_vf(i)=0D0
          text_vz(i)=0D0
          ! vortex direction
          text_lr(i)=0D0
          text_lf(i)=0D0
          text_lz(i)=1D0
          text_w(i)=2D0*omega
          ! fortex-free part
          if (rr > text_r*SQRT(omega_v/omega)) then
            text_vf(i) = rr*(omega - omega_v*(text_r/rr)**2)
            text_w(i)=0D0
          endif
        enddo
      end

! Uniform vortex profile
      subroutine text1r_set_vortex_uniform(omega, omega_v)
        implicit none
        include 'text1r.fh'
        integer i
        real*8 omega,omega_v,rr
        do i=1,MAXN
          rr = dble(i-1)/dble(MAXN-1)*text_r
          ! flow velocity
          text_vr(i)=0D0
          text_vf(i)=(omega-omega_v)*rr
          text_vz(i)=0D0
          ! vortex direction
          text_lr(i)=0D0
          text_lf(i)=0D0
          text_lz(i)=1D0
          text_w(i)=2D0*omega_v
        enddo
      end

! Twisted-state velocity profile
      subroutine text1r_set_vortex_twisted(omega, kr)
        implicit none
        include 'text1r.fh'
        integer i
        real*8 omega,kr,rr,r
        do i=1,MAXN
          rr = dble(i-1)/dble(MAXN-1)*text_r
          r  = text_r
          text_vr(i)=0D0
          text_vf(i)=omega*rr-(kr**2/LOG(1+kr**2))*omega*rr/(1+(kr*rr/r)**2)
          text_vz(i)=-omega*r*(kr/(LOG(1+kr**2)*(1+(kr*rr/r)**2))-1/kr)
          text_lr(i)=0D0
          text_lf(i)=(kr*rr/r)/SQRT(1+(kr*rr/r)**2)
          text_lz(i)=1D0/SQRT(1+(kr*rr/r)**2)
          text_w(i)=2D0*omega*(kr**2/LOG(1+kr**2))*(1+(kr*rr/r)**2)**(-1.5)
        enddo
      end

! Print texture
      subroutine text1r_print(fname)
        implicit none
        integer fd,i
        real*8 r2d
        character fname*(*)
        include 'text1r.fh'
  101   format (A,E12.4,A)
  102   format (A,I4)
        fd=100
        write(*,*) '<<<', fname(1:10), '>>>', LEN(fname)
        open(fd, file='result.dat')
        write (fd,101) '# Texture parameters: '
        write (fd,102) '#  Number of points n = ', MAXN
        write (fd,101) '#  Cell radius r = ', text_r,   ' cm'
        write (fd,101) '#  Mag.field   H = ', text_H,   ' G'
        write (fd,101) '#          (nu_0 = ', text_H*20.0378/2/acos(0D0), ' kHz)'
        write (fd,101) '#  lambda_D      = ', text_ld,  ' erg/cm3'
        write (fd,101) '#  Field par. a  = ', text_a,   ' erg/cm^3 1/G^2'
        write (fd,101) '#  lambda_G1     = ', text_lg1, ' erg/cm'
        write (fd,101) '#  lambda_G2     = ', text_lg2, ' erg/cm'
        write (fd,101) '#        (delta  = ', text_lg1/text_lg2-2D0, ''
        write (fd,101) '#  lambda_HV     = ', text_lhv, ' erg/cm3 1/G2 1/(cm/s)^2'
        write (fd,101) '#  lambda_SG     = ', text_lsg, ' ??'
        write (fd,101) '#  lambda/Omega  = ', text_lo,  ' s/rad'
        write (fd,101) '#  surface par d = ', text_d,   ' erg/cm^2 1/G^2'
        write (fd,'(A)') '#'
        r2d=90D0/acos(0D0)
  103   format ('# ' A6 '  ' A9 ' ' A9 '  ' A9                           &
      &              '  ' A9 ' ' A9 ' ' A9                               &
      &              '  ' A9 ' ' A9 ' ' A9                               &
      &              '  ' A9)
  104   format (F7.5 '  ' F10.5 ' ' F9.5 '  ' F9.5                       &
      &              '  ' F9.5 ' ' F9.5 ' ' F9.5                         &
      &              '  ' F9.5 ' ' F9.5 ' ' F9.5                         &
      &              '  ' F9.5)
        write(fd, 103) 'r,cm', 'a_n,deg', 'b_n,deg', 'apsi',             &
      &                'vr,cm/s', 'vz,cm/s', 'vf,cm/s',                  &
      &                'lr', 'lz', 'lf','w'
        do i=1,MAXN
          write (fd,104) text_rr(i),                                     &
      &                  text_an(i)*r2d, text_bn(i)*r2d,                 &
      &                  text_apsi(i),                                   &
      &                  text_vr(i),text_vz(i),text_vf(i),               &
      &                  text_lr(i),text_lz(i),text_lf(i),               &
      &                  text_w(i)
        enddo
        close(fd)
      end

! Vary alpha_n and beta_n to get equilibrium texture
      subroutine text1r_minimize(msglev)
        implicit none
        include 'text1r.fh'
        integer maxnpar, lw, ierror, msglev
        parameter (maxnpar = 2*MAXN-2)
        parameter (lw = 14*maxnpar)
        real*8 x(maxnpar), g(maxnpar), f
        real*8 w(lw)
        external text1r_mfunc
        call text1r_text2x(MAXN, text_an, text_bn, x)
        call tn(ierror,maxnpar,x,f,g,w,lw,text1r_mfunc,msglev)
        call text1r_x2text(MAXN, text_an, text_bn, x)
      end

!      subroutine text1r_minimize_btn()
!        implicit none
!        include 'text1r.fh'
!        integer iflag
!        integer maxnpar=2*MAXN-1
!        integer nprocs=8
!        integer lw= 3*maxnpt + 3*nprocs &
!                  + 4*nprocs*nprocs + 7 *(maxnpt*nprocs)
!        real*8 x(maxnpar), g(maxnpar)
!        real*8 w(lw)
!        call ab2x(MAXN, text_an, text_bn, x)
!        call btnez(n,x,f,g, w, lw, text1r_mfunc, iflag)
!        call x2ab(MAXN, text_an, text_bn, x)
!      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Convert alpha and beta vectors to x vector
! x(1:2n-1) = [alpha(1:n) beta(2:n)]
! beta(1) is thrown away (it is 0)
      subroutine text1r_text2x(n,a,b, x)
        implicit none
        integer i,n
        real*8 a(n),b(n),x(2*n-2)
        do i=2,n
          x(i-1)       = a(i)
          x(n-1 + i-1) = b(i)
        enddo
      end

! Convert x to alpha, beta vectors
! beta(0) is set to 0
      subroutine text1r_x2text(n,a,b, x)
        implicit none
        integer i,n
        real*8 a(n),b(n),x(2*n-2)
        do i=2,n
          a(i) = x(i-1)
          b(i) = x(i-1+n-1)
        enddo
        a(1)=a(2)
        b(1)=0D0
      end

! Wrapper for egrad function for using in the TN.
! Calculate f and g from x values.
! x as array of both alpha and beta values
! g is array of both ga, gb
      subroutine text1r_mfunc(nx,x,f,g)
        !! Wrapper for egrad function for using in the TN.
        implicit none
        integer i,nx,n
        real*8 f, x(nx),g(nx)
        real*8 a(nx/2+1), b(nx/2+1)
        real*8 ga(nx/2+1), gb(nx/2+1)
        n=nx/2+1
        call text1r_x2text(n, a,b, x)
        call text1r_egrad(n,a,b,f,ga,gb)
        call text1r_text2x(n, ga,gb, g)
      end

!!! Calculate free energy E and derivatives dE/da(i), dE/db(i)
!
! E = Int e(r, a(r),b(r), da/dr, db/dr,...) r dr
! Gaussian quadrature is used for calculating integral. For each
! mesh point i = 1:N-1 energies em,ep is calculated in points
! i+sm, i+sp (sm/sp = (3 +/- 3/sqrt(3))/6), and
! (em*rm + ep*rp)*dr/2 is added to the sum
!
! Derivatives dE/da(i), dE/db(i) are also needed.
! Change DA of a(i) (or b(i)) affects 4 terms in E intergral:
!   at i-sp, i-sm (for i!=0), i+sp, i+sm (for i!=n)
! Changes of a in these points are DA*sm, DA*sp, DA*sp, DA*sm
! Changes of a' in these points are DA/dr, DA/dr, -DA/dr, -DA/dr
! We need to calculate dE/DA = Sum(dE/da * da + dE/da' * da')/DA
! We also need r*dr/2 factor as in energy calculation
!
! Strightforward approach is to calculate sum for these 4 points
!   (Ea*sm*dr - Eda)*r/2 for i+sp, i!=n
!   (Ea*sp*dr - Eda)*r/2 for i+sm, i!=n
!   (Ea*sp*dr + Eda)*r/2 for i-sm, i!=0
!   (Ea*sm*dr + Eda)*r/2 for i-sp, i!=0
! but we can calculate E* only in two points instead of 4
! and add some terms to both dE/da(i) and dE/da(i+1)
      subroutine text1r_egrad(n,a,b,e,ga,gb)

        implicit none
        include 'text1r.fh'
        integer i,n
        real*8 a(n), b(n), ga(n), gb(n)
        real*8 rp,rm,bp,bm,ap,am,e
        real*8 apsip, vzp,vrp,vfp,lzp,lrp,lfp,wp
        real*8 apsim, vzm,vrm,vfm,lzm,lrm,lfm,wm
        real*8 da,db
        real*8 E0,Ea,Eb,Eda,Edb
        real*8 dx,sp,sm

        ! points for Gaussian quadrature
        sp = (3D0 + sqrt(3D0))/6D0
        sm = (3D0 - sqrt(3D0))/6D0
        ! grid step
        dx = 1D0/(n-1)

        do i=1,n
           ga(i)=0D0
           gb(i)=0D0
        enddo
        e=0D0

        do i=1,n-1
          rp = (i-1+sp)*dx
          rm = (i-1+sm)*dx
          ap = sp*a(i+1)+sm*a(i)
          am = sm*a(i+1)+sp*a(i)
          bp = sp*b(i+1)+sm*b(i)
          bm = sm*b(i+1)+sp*b(i)

          da = (a(i+1)-a(i))/dx
          db = (b(i+1)-b(i))/dx

          apsip=sp*text_apsi(i+1)+sm*text_apsi(i)
          apsim=sm*text_apsi(i+1)+sp*text_apsi(i)

          vzp=sp*text_vz(i+1)+sm*text_vz(i)
          vzm=sm*text_vz(i+1)+sp*text_vz(i)
          vrp=sp*text_vr(i+1)+sm*text_vr(i)
          vrm=sm*text_vr(i+1)+sp*text_vr(i)
          vfp=sp*text_vf(i+1)+sm*text_vf(i)
          vfm=sm*text_vf(i+1)+sp*text_vf(i)
          lzp=sp*text_lz(i+1)+sm*text_lz(i)
          lzm=sm*text_lz(i+1)+sp*text_lz(i)
          lrp=sp*text_lr(i+1)+sm*text_lr(i)
          lrm=sm*text_lr(i+1)+sp*text_lr(i)
          lfp=sp*text_lf(i+1)+sm*text_lf(i)
          lfm=sm*text_lf(i+1)+sp*text_lf(i)
          wp=sp*text_w(i+1)+sm*text_w(i)
          wm=sm*text_w(i+1)+sp*text_w(i)

          call text1r_ebulk(rp, ap,bp,da,db, apsip, vzp,vrp,vfp, &
                   lzp,lrp,lfp, wp, E0,Ea,Eb,Eda,Edb)
          ga(i) = ga(i) + (Ea*sm*dx - Eda)*rp/2.0
          gb(i) = gb(i) + (Eb*sm*dx - Edb)*rp/2.0
          ga(i+1) = ga(i+1) + (Ea*sp*dx + Eda)*rp/2.0
          gb(i+1) = gb(i+1) + (Eb*sp*dx + Edb)*rp/2.0
          e = e + rp*e0*0.5*dx

          call text1r_ebulk(rm, am,bm,da,db, apsim, vzm,vrm,vfm, &
                   lzm,lrm,lfm, wm, E0,Ea,Eb,Eda,Edb)
          ga(i) = ga(i) + (Ea*sp*dx - Eda)*rm/2.0
          gb(i) = gb(i) + (Eb*sp*dx - Edb)*rm/2.0
          ga(i+1) = ga(i+1) + (Ea*sm*dx + Eda)*rm/2.0
          gb(i+1) = gb(i+1) + (Eb*sm*dx + Edb)*rm/2.0
          e = e + rm*E0*0.5*dx

        enddo
        ! surface terms
        call text1r_esurf(a(n),b(n),E0,Ea,Eb)
        e = e + E0
        ga(n) = ga(n) + Ea
        gb(n) = gb(n) + Eb
      end


!! Calculate He3 bulk energy E and derivatives dE/da, dE/db, dE/da', dE/db'
!! in 1d radial coordinated as a function of r, n-vector angles a and b,
!! and derivatives a'=da/dr, b'=db/dr.
!! Energy is divided by (a H^2).
!! parameters used:
!!   chia*(nub/nu0)^2, for non-zero apsi
!!   lo for non-zero rotation
!!   vd for non-zero flow
!!   de and xi
      subroutine text1r_ebulk(r,a,b,da,db, &
                         apsi, vz,vr,vf, lz,lr,lf, w, &
                          E,Ea,Eb,Eda,Edb)
        implicit none
        include 'text1r.fh'
        real*8 r,a,b,da,db,E,Ea,Eb,Eda,Edb
        real*8 apsi, vz,vr,vf, lz,lr,lf, w
        real*8 nz,nr,nf, rzz,rzr,rzf
        real*8 sin_a, sin_b, cos_a, cos_b, cos2b, sin2b
        real*8 con1, con2, help, c,s

        real*8 s3,s5,dar,de,xir,vd

        dar = text_d / (text_a*text_r)
        de  = text_lg1/text_lg2 - 2D0
        xir  = sqrt(65D0/8D0 * text_lg2/text_a)/ text_H / text_r
        vd  = dsqrt(0.4D0 * text_a / text_lhv)

        s3 = sqrt(3D0)
        s5 = sqrt(5D0)

        cos_a = cos(a)
        sin_a = sin(a)
        cos_b = cos(b)
        sin_b = sin(b)
        cos2b = cos(2*b)
        sin2b = sin(2*b)

        nr=-sin_b*cos_a
        nf=sin_b*sin_a
        nz=cos_b

        c=-0.25D0 !\cos\theta
        s=SQRT(15D0)/4D0 !\sin\theta

        ! H*Rij = cos(t) H + (1-cos(t))(H n) n - sin(t)(H x n)
        rzr=(1-c)*nz*nr-s*nf ! H*Rij/abs(H)
        rzf=(1-c)*nz*nf+s*nr
        rzz=c+(1-c)*nz**2

        E = 0
        Ea = 0
        Eb = 0
        Eda = 0
        Edb = 0

        ! magnetic free energy F_DH = -a * (n H)^2
        E = E + sin_b**2
        Eb = Eb + sin2b

        ! spin-orbit free energy
        help = 15D0 * text_ld / text_a / text_h**2
        E = E + help * (apsi * sin_b)**2
        Eb = Eb + help * sin2b * apsi**2

        ! flow free energy F_HV
        ! v_d = sqrt(2/5 a/lhv)
        E = E - 2/5 / (vd**2) * (rzr*vr+rzf*vf+rzz*vz)**2
        help = vr*(-(1-c)*cos2b*cos_a - s*cos_b*sin_a) &
             + vf*( (1-c)*cos2b*sin_a - s*cos_b*cos_a) &
             + vz*(-(1-c)*sin2b)
        Eb = Eb - 4/5 / (vd**2) * help * (rzr*vr+rzf*vf+rzz*vz)

        help = vr*((1-c)*sin_b*cos_b*sin_a - s*sin_b*cos_a) &
             + vf*((1-c)*sin_b*cos_b*cos_a + s*sin_b*sin_a)
        Ea = Ea - 4/5 / (vd**2) * help *(rzr*vr+rzf*vf+rzz*vz)

        ! vortex free energy F_LH
        E = E + text_lo * w/5 * (rzr*lr+rzf*lf+rzz*lz)**2

        help = lr*(-(1-c)*cos2b*cos_a - s*cos_b*sin_a) &
             + lf*( (1-c)*cos2b*sin_a - s*cos_b*cos_a) &
             + lz*(-(1-c)*sin2b)
        Eb = Eb + text_lo*w*2/5 * help * (rzr*lr+rzf*lf+rzz*lz)
        help = lr*((1D0-c)*sin_b*cos_b*sin_a - s*sin_b*cos_a) &
             + lf*((1D0-c)*sin_b*cos_b*cos_a + s*sin_b*sin_a)
        Ea = Ea + text_lo*w*2/5 * help * (rzr*lr+rzf*lf+rzz*lz)

        ! bending free energy F_G
        con1 = 4*(4+de)*xir**2/13

        E = E + con1*(db**2 + (sin_b**2)*da**2 + (sin_b**2)/r**2) ! (\nabla n)^2 (?)
        Eda = Eda + con1*2*da*sin_b**2
        Edb = Edb + con1*2*db
        Eb = Eb + con1 * 2*sin_b*cos_b*(da**2 + 1/r**2)

        con2 = -(2+de)*xir**2/26
        help=(s5*sin_a-s3*cos_b*cos_a)*db + &
             (s5*cos_b*cos_a+s3*sin_a)*sin_b*da + &
             (s5*cos_b*sin_a-s3*cos_a)*sin_b/r  ! s3 \div n + s5 n \rot n (?)
        E = E + con2 * help**2

        Eda = Eda + 2*con2*help * (s5*cos_b*cos_a + s3*sin_a)*sin_b
        Edb = Edb + 2*con2*help * (s5*sin_a - s3*cos_b*cos_a)

        Ea = Ea + 2*con2*help* &
          ( (s5*cos_a + s3*cos_b*sin_a)*db &
          - (s5*cos_b*sin_a - s3*cos_a)*sin_b*da &
          + (s5*cos_b*cos_a + s3*sin_a)*sin_b/r)
        Eb = Eb + 2*con2*help* &
          ( (s3*db - s5*sin_b*da)*sin_b*cos_a &
          + (s5*cos_b*cos_a + s3*sin_a)*cos_b*da &
          + (s5*cos_b*sin_a - s3*cos_a)*cos_b/r &
          - s5*sin_b*sin_a*sin_b/r)
      end


! Calculate E, dE/da, dE/db, dE/da', dE/db' at surface
! parameters used: dar,xir,lsg
      subroutine text1r_esurf(a,b,E,Ea,Eb)
        implicit none
        include 'text1r.fh'
        integer i
        real*8 a,b,E,Ea,Eb
        real*8 nr,nf,nz
        real*8 sin_a, sin_b, cos_a, cos_b, sin2b,cos2b

        real*8 s3,s5,de, xir,dar

        dar = text_d / (text_a*text_r)
        de  = text_lg1/text_lg2 - 2D0
        xir = sqrt(65D0/8D0 * text_lg2/text_a)/ text_H / text_r

        s3 = sqrt(3D0)
        s5 = sqrt(5D0)

        E=0D0
        Ea=0D0
        Eb=0D0
        cos_a = cos(a)
        sin_a = sin(a)
        cos_b = cos(b)
        sin_b = sin(b)
        sin2b = sin(2*b)
        cos2b = cos(2*b)
        nr=-sin_b*cos_a
        nf=sin_b*sin_a
        nz=cos_b

        ! dar is d/a/r
        E = E - 5D0*dar*(s5*nz*nr-s3*nf)**2/16D0
        Eb = Eb + 5D0*dar*(s5*nz*nr-s3*nf)*(s5*cos2b*cos_a+s3*cos_b*sin_a)/8D0
        Ea = Ea - 5D0*dar*(s5*nz*nr-s3*nf)*(s5*nz*nf + s3*nr)/8D0

        ! from bending free energy
        ! xir = xi_H = sqrt(65 lg2 /(8 a)) / H/R
        E = E + 4D0*(2D0+de)*xir**2*sin_b**2/13D0
        Eb = Eb + 4D0*(2D0+de)*xir**2*sin2b/13D0

        ! lsg ~ 3
        E = E - 2D0*text_lsg*xir**2*sin_b**2/13D0
        Eb = Eb - 2D0*text_lsg*xir**2*sin2b/13D0
      end

