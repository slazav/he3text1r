! Data block with parameter structure
      block data text1r_data_block
        implicit none
        include 'text1r.fh'
      end block data

! Check that number of points is less then array sizes
      function check_size()
        include 'text1r.fh'
        logical check_size
        check_size = text_n.gt.MAXN
        if (check_size) write (*,*) &
          'Error: number of points is larger then MAXN parameter'
      end

! Initialize parameter structure
      subroutine text1r_init(ttc,p,nu0,r, n, itype)
        implicit none
        real*8 ttc,p,nu0,r
        integer n, itype, i, zz,hh
        logical check_size
        include 'text1r.fh'
        include 'he3.f90h'

        text_a = he3_text_a(ttc,p)
        text_lg1 = he3_grad_lg1(ttc,p)
        text_lg2 = he3_grad_lg2(ttc,p)
        text_lhv = he3_text_lhv(ttc,p)
        text_lsg = 3D0
        text_ld  = he3_ld(ttc,p)
        text_lo  = 0D0
        text_d = he3_text_d(ttc,p)
        text_r = r
        text_h = const_2pi*nu0/he3_gyro
        text_n = n
        text_err = 0
        text_bnd = 0
        text_abnd = acos(0.5D0);       ! 60 deg
        text_bbnd = asin(sqrt(0.8D0)); ! 63.4 deg
        text_energy = 0D0
        if (check_size()) return

        ! initial conditions for alpha_n, beta_n
        hh = itype/2
        zz = (itype - hh*2)*2-1 ! parity of itype -1 or 1
        do i=1,text_n
          text_an(i) = -dble(zz)*acos(0.5D0)
          text_bn(i) = 2D0*acos(0D0)*dble(hh) + &
                       acos(-dble(zz)*1D0/sqrt(5D0)) * &
                       dble(i-1)/dble(text_n-1)
          if (text_bnd.ne.0) then
            text_rr(i)=text_r*dble(i-1)/dble(text_n)
          else
            text_rr(i)=text_r*dble(i-1)/dble(text_n-1)
          endif
          text_bm(i)=0D0
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
        logical check_size
        real*8 omega,omega_v,rr

        if (check_size()) return
        do i=1,text_n

          if (text_bnd.ne.0) then
            rr = dble(i-1)/dble(text_n)*text_r
          else
            rr = dble(i-1)/dble(text_n-1)*text_r
          endif

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
        logical check_size
        real*8 omega,omega_v,rr

        if (check_size()) return
        do i=1,text_n
          if (text_bnd.ne.0) then
            rr = dble(i-1)/dble(text_n)*text_r
          else
            rr = dble(i-1)/dble(text_n-1)*text_r
          endif
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
        logical check_size
        real*8 omega,kr,rr,r

        if (check_size()) return
        do i=1,text_n
          if (text_bnd.ne.0) then
            rr = dble(i-1)/dble(text_n)*text_r
          else
            rr = dble(i-1)/dble(text_n-1)*text_r
          endif
          r  = text_r
          text_vr(i)=0D0
          text_vf(i)=omega*rr-(kr**2/dlog(1D0+kr**2))*omega*rr/(1D0+(kr*rr/r)**2)
          text_vz(i)=-omega*r*(kr/(dlog(1D0+kr**2)*(1D0+(kr*rr/r)**2))-1D0/kr)
          text_lr(i)=0D0
          text_lf(i)=(kr*rr/r)/dsqrt(1D0+(kr*rr/r)**2)
          text_lz(i)=1D0/dsqrt(1D0+(kr*rr/r)**2)
          text_w(i)=2D0*omega*(kr**2/dlog(1D0+kr**2))*(1D0+(kr*rr/r)**2)**(-1.5D0)
        enddo
      end

! Print texture
      subroutine text1r_print(fname)
        implicit none
        integer fd,i
        logical check_size
        real*8 r2d
        character fname*(*)
        include 'text1r.fh'
  101   format (A,E12.4,A)
  102   format (A,I4)
        fd=100
        open(fd, file='result.dat')
        write (fd,101) '# Texture parameters: '
        write (fd,102) '#  Number of points n = ', text_n
        write (fd,101) '#  Cell radius r = ', text_r,   ' cm'
        write (fd,101) '#  Mag.field   H = ', text_H,   ' G'
        write (fd,101) '#          (nu_0 = ', text_H*20.0378D0/2D0/acos(0D0), ' kHz)'
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
        if (check_size()) return
  103   format ('# ' A6 '  ' A9 ' ' A9 '  ' A9                           &
      &              '  ' A9 ' ' A9 ' ' A9                               &
      &              '  ' A9 ' ' A9 ' ' A9                               &
      &              '  ' A9)
  104   format (F7.5 '  ' F10.5 ' ' F9.5 '  ' F9.5                       &
      &              '  ' F9.5 ' ' F9.5 ' ' F9.5                         &
      &              '  ' F9.5 ' ' F9.5 ' ' F9.5                         &
      &              '  ' F9.5)
        write(fd, 103) 'r,cm', 'a_n,deg', 'b_n,deg', 'bm,deg',           &
      &                'vr,cm/s', 'vz,cm/s', 'vf,cm/s',                  &
      &                'lr', 'lz', 'lf','w'
        do i=1,text_n
          write (fd,104) text_rr(i),                                     &
      &                  text_an(i)*r2d, text_bn(i)*r2d,                 &
      &                  text_bm(i)*r2d,                                 &
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
        integer maxnpar, lw, msglev
        parameter (maxnpar = 2*MAXN-2)
        parameter (lw = 14*maxnpar)
        real*8 x(maxnpar), ex(maxnpar)
        real*8 w(lw)
        logical check_size
        external text1r_mfunc
        integer nn,i

        ! additional parameters (default values, do not change)
        integer MAXIT, MAXFUN
        real*8 ETA,STEPMX,ACCRCY,XTOL,MCHPR1

        ! update text_rr (it can depend on text_bnd)
        do i=1,text_n
          if (text_bnd.ne.0) then
            text_rr(i)=text_r*dble(i-1)/dble(text_n)
          else
            text_rr(i)=text_r*dble(i-1)/dble(text_n-1)
          endif
        enddo

        nn=2*text_n-2

        MAXIT = nn/2
        IF (MAXIT .GT. 50) MAXIT = 50
        IF (MAXIT .LE. 0) MAXIT = 1
        MAXFUN = 150*nn
        ETA = .25D0
        STEPMX = 1.D1
        ACCRCY = 1.D2*MCHPR1()
        XTOL = DSQRT(ACCRCY)

        if (check_size()) return
        call text1r_text2x(text_n, text_an, text_bn, x)
        do i=1,20 ! we need several runs to catch small energy changes for flat textures
          call lmqn(text_err, nn,x,text_energy,ex,w,lw, &
                text1r_mfunc,msglev, &
                MAXIT, MAXFUN, ETA, STEPMX, ACCRCY, XTOL)
          if (text_err.eq.0) goto 10
        enddo
        write(*,*) "ERROR CODE = ", text_err, " after ", i, " iterations"
10      call text1r_x2text(text_n, text_an, text_bn, x)

!       derivatives
        text_db0 = text_bn(2)/text_rr(2)
        text_db1 = (text_bn(text_n)-text_bn(text_n-1))/ &
                   (text_rr(text_n)-text_rr(text_n-1))
        text_bmax = text_bn(text_n)

!       the l vector
        do i=2,text_n
          call n2l(text_an(i), text_bn(i), text_al(i), text_bl(i))
        enddo
        text_al(1)=text_al(2)
        text_bl(1)=0D0
      end

!      calculate the l vector from n
       subroutine n2l(an,bn,al,bl)
         implicit none
         real*8 an,bn,al,bl
         real*8 nr,nf,nz,lr,lf,lz
         nr=-sin(bn)*cos(an)
         nf=sin(bn)*sin(an)
         nz=cos(bn)
         lr = 1.25D0 * nz*nr - sqrt(15D0/16D0) * nf
         lf = 1.25D0 * nz*nf + sqrt(15D0/16D0) * nr
         lz = -0.25D0 + 1.25D0 * nz**2
         bl = acos(lz)
         al = atan2(lf,-nr)
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
!
!        if (text_n.gt.MAXN) then
!          write (*,*) 'Error: number of points is larger then MAXN parameter'
!          return
!        endif
!
!        call ab2x(text_n, text_an, text_bn, x)
!        call btnez(n,x,f,g, w, lw, text1r_mfunc, iflag)
!        call x2ab(text_n, text_an, text_bn, x)
!      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Various checks
      subroutine text1r_selfcheck(msglev)
        implicit none
        include 'text1r.fh'
        integer msglev, i

        logical check_size
        if (check_size()) return

        call text1r_conv_selfcheck( &
          text_n,text_an,text_bn, 1D-10)
        do i=2,text_n
          call text1r_ebulk_selfcheck( &
            text_an(i), text_bn(i), text_bm(i), &
                 text_vz(i), text_vr(i), text_vf(i), &
                 text_lz(i), text_lr(i), text_lf(i), &
                 text_w(i), 1D-3, 1D-4)
          call text1r_esurf_selfcheck( &
            text_an(i),text_bn(i), 1D-3, 1D-4)
          call text1r_egrad_selfcheck( &
            text_rr(i),text_an(i),text_bn(i), 1D-3, 1D-4)
        enddo
        call text1r_eint_selfcheck( &
          text_n,text_an,text_bn, 1D-3, 1D-2)
        call text1r_mfunc_selfcheck( &
          text_n,text_an,text_bn, 1D-3, 1D-2)

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Convert alpha and beta vectors to x vector
! x(1:2n-1) = [alpha(1:n) beta(2:n)]
! beta(1) is thrown away (it is 0)
      subroutine text1r_text2x(n,a,b, x)
        implicit none
        integer i,n
        real*8 a(n),b(n),x(2*n-2)
        do i=2,n
          x(i-1)       = -sin(b(i))*cos(a(i))/(1D0+cos(b(i)))
          x(n-1 + i-1) =  sin(b(i))*sin(a(i))/(1D0+cos(b(i)))
        enddo
      end

! Convert dE/dAi, dE/dBi to dE/dXi
      subroutine text1r_text2ex(n,a,b, ea,eb,ex)
        implicit none
        integer i,n
        real*8 a(n),b(n),ea(n),eb(n),ex(2*n-2)
        real*8 dau,dav,dbu,dbv,cb
        do i=2,n
          cb=cos(b(i))
          if (1D0-cb.gt.1D-30) then
            dau =  sin(b(i))*sin(a(i))/(1D0-cb)
            dav =  sin(b(i))*cos(a(i))/(1D0-cb)
            dbu = -sin(b(i))*cos(a(i))*sqrt((1D0+cb)/(1D0-cb))
            dbv =  sin(b(i))*sin(a(i))*sqrt((1D0+cb)/(1D0-cb))
          else
            dau =  sin(a(i)) * 2D0/1D-30
            dav =  cos(a(i)) * 2D0/1D-30
            dbu = -2D0*cos(a(i))
            dbv =  2D0*sin(a(i))
          endif
          ex(i-1)       = ea(i)*dau + eb(i)*dbu
          ex(n-1 + i-1) = ea(i)*dav + eb(i)*dbv
        enddo
      end

! Convert x to alpha, beta vectors
! beta(0) is set to 0
      subroutine text1r_x2text(n,a,b, x)
        implicit none
        integer i,n
        real*8 a(n),b(n),x(2*n-2),u,v
        do i=n,2,-1
          u = x(i-1)
          v = x(i-1+n-1)
          b(i) = acos((1D0-u**2-v**2)/(1D0+u**2+v**2))
          if (abs(b(i)).gt.1D-10) then
            a(i) = 4D0*atan(1D0)-atan2(v, u)
          else if (i.lt.n) then
            a(i)=a(i+1)
          else
            a(i)=0D0
          endif
!          write(*,*) a(i),b(i),u,v
        enddo
        a(1)=a(2)
        b(1)=0D0
      end

! Self test for convert functions
      subroutine text1r_conv_selfcheck(n, a,b, e)
        implicit none
        include 'text1r.fh'
        integer n,i
        real*8 a(n),b(n),x(2*MAXN-2),e
        real*8 an(MAXN),bn(MAXN)

        real*8 E1,Ea1(MAXN),Eb1(MAXN)
        real*8 E2,Ea2(MAXN),Eb2(MAXN)
        real*8 der1,der2

        ! check conversions a,b -> x - a,b
        call text1r_text2x(n, a, b, x)
        call text1r_x2text(n, an, bn, x)
        do i=1,n
          if (dabs(a(i) - an(i)) > e) then
            write (*,*) 'text1r_conv_selfcheck failed for A->X->A conversion:'
            write (*,*) 'i: ', i, ' A1: ', an(i), ' A2: ', a(i)
          endif
          if (dabs(b(i) - bn(i)) > e) then
            write (*,*) 'text1r_conv_selfcheck failed for B->X->B conversion:'
            write (*,*) 'i: ', i, ' B1: ', bn(i), ' B2: ', b(i)
          endif
        enddo
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Wrapper for egrad function for using in the TN.
! Calculate f and g from x values.
! x as array of both alpha and beta values
! g is array of both ga, gb
      subroutine text1r_mfunc(nx,x,e,ex)
        !! Wrapper for egrad function for using in the TN.
        implicit none
        include 'text1r.fh'
        integer i,nx,n
        real*8 e, x(nx),ex(nx)
        real*8 a(MAXN), b(MAXN)
        real*8 ea(MAXN), eb(MAXN)
        n=nx/2+1
        call text1r_x2text(n, a,b, x)
        call text1r_eint(n,a,b,e,ea,eb)
        call text1r_text2ex(n, a,b,ea,eb, ex)
      end

! Self test for mfunc derivatives
      subroutine text1r_mfunc_selfcheck(n, a,b, d, e)
        implicit none
        include 'text1r.fh'
        integer n,i
        real*8 a(n),b(n),d,e
        real*8 x(2*MAXN-2)
        real*8 E1,Ex1(2*MAXN-2)
        real*8 E2,Ex2(2*MAXN-2)
        real*8 der1,der2

        call text1r_text2x(n, a, b, x)
        call text1r_mfunc(2*n-2, x, e1, ex1)
        do i=1,2*n-2

          x(i) = x(i) + d
          call text1r_mfunc(2*n-2, x, e2, ex2)
          x(i) = x(i) - d
          der1=(e2-e1)/d
          der2=(ex2(i)+ex1(i))/2D0
          if ( dabs( der1/der2 - 1D0 ) > e ) then
            write (*,*) 'text1r_mfunc_selfcheck failed for derivative dE/dXi:'
            write (*,*) 'i: ', i, ' (e2-e1)/dx: ', der1, ' ex*dx: ', der2
          endif
        enddo
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Calculate total free energy E and derivatives dE/da(i), dE/db(i)
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
! We also need r*dr/2 eactor as in energy calculation
!
! Strightforward approach is to calculate sum for these 4 points
!   (Ea*sm*dr - Eda)*r/2 for i+sp, i!=n
!   (Ea*sp*dr - Eda)*r/2 for i+sm, i!=n
!   (Ea*sp*dr + Eda)*r/2 for i-sm, i!=0
!   (Ea*sm*dr + Eda)*r/2 for i-sp, i!=0
! but we can calculate E* only in two points instead of 4
! and add some terms to both dE/da(i) and dE/da(i+1)
      subroutine text1r_eint(n,a,b,ei,eai,ebi)

        implicit none
        include 'text1r.fh'
        integer i,n
        real*8 a(n), b(n), eai(n), ebi(n)
        real*8 rp,rm,bp,bm,ap,am,ei
        real*8 da,db
        real*8 E0,Ea,Eb,Eda,Edb
        real*8 E1, E2,E1a,E1b,E2a,E2b
        real*8 dx,sp,sm

        ! points for Gaussian quadrature
        sp = (3D0 + sqrt(3D0))/6D0
        sm = (3D0 - sqrt(3D0))/6D0
        ! grid step
        dx = 1D0/dble(n-1)

        do i=1,n
           eai(i)=0D0
           ebi(i)=0D0
        enddo
        ei=0D0

        do i=1,n-1
          ! gaussian quadrature points
          rp = (dble(i-1)+sp)*dx
          rm = (dble(i-1)+sm)*dx

          ! bulk energy at the i point
          if (i.eq.1) then
            call text1r_ebulk(a(i),b(i), text_bm(i), &
                 text_vz(i), text_vr(i), text_vf(i), &
                 text_lz(i), text_lr(i), text_lf(i), &
                 text_w(i), E1, E1a, E1b)
          else
            E1=E2
            E1a=E2a
            E1b=E2b
          endif
          ! bulk energy at the i+1 point
          call text1r_ebulk(a(i+1),b(i+1), text_bm(i+1), &
               text_vz(i+1), text_vr(i+1), text_vf(i+1), &
               text_lz(i+1), text_lr(i+1), text_lf(i+1), &
               text_w(i+1), E2, E2a, E2b)

          ! interpolate bulk energy and derivatives to rp,rm points
          ei = ei + rp*(sp*E2 + sm*E1)*0.5D0*dx
          eai(i) = eai(i) + E1a*sm*dx*rp*0.5D0
          ebi(i) = ebi(i) + E1b*sm*dx*rp*0.5D0
          eai(i+1) = eai(i+1) + E2a*sp*dx*rp*0.5D0
          ebi(i+1) = ebi(i+1) + E2b*sp*dx*rp*0.5D0

          ei = ei + rm*(sm*E2 + sp*E1)*0.5D0*dx
          eai(i) = eai(i) + E1a*sp*dx*rm*0.5D0
          ebi(i) = ebi(i) + E1b*sp*dx*rm*0.5D0
          eai(i+1) = eai(i+1) + E2a*sm*dx*rm*0.5D0
          ebi(i+1) = ebi(i+1) + E2b*sm*dx*rm*0.5D0

          ! calculate gradient energy in rp, rm points
          ap = sp*a(i+1)+sm*a(i)
          am = sm*a(i+1)+sp*a(i)
          bp = sp*b(i+1)+sm*b(i)
          bm = sm*b(i+1)+sp*b(i)
          da = (a(i+1)-a(i))/dx
          db = (b(i+1)-b(i))/dx

          call text1r_egrad(rp, ap,bp,da,db, E0,Ea,Eb,Eda,Edb)
          ei = ei + rp*E0*0.5D0*dx
          eai(i) = eai(i) + (Ea*sm*dx - Eda)*rp*0.5D0
          ebi(i) = ebi(i) + (Eb*sm*dx - Edb)*rp*0.5D0
          eai(i+1) = eai(i+1) + (Ea*sp*dx + Eda)*rp*0.5D0
          ebi(i+1) = ebi(i+1) + (Eb*sp*dx + Edb)*rp*0.5D0

          call text1r_egrad(rm, am,bm,da,db, E0,Ea,Eb,Eda,Edb)
          ei = ei + rm*E0*0.5D0*dx
          eai(i) = eai(i) + (Ea*sp*dx - Eda)*rm*0.5D0
          ebi(i) = ebi(i) + (Eb*sp*dx - Edb)*rm*0.5D0
          eai(i+1) = eai(i+1) + (Ea*sm*dx + Eda)*rm*0.5D0
          ebi(i+1) = ebi(i+1) + (Eb*sm*dx + Edb)*rm*0.5D0

        enddo


        !! strict boundary conditions: add one more point (n+1)
        !! with fixed angles text_abnd, text_bbnd
        if (text_bnd.ne.0) then
          ! gaussian quadrature points
          rp = (dble(n-1)+sp)*dx
          rm = (dble(n-1)+sm)*dx

          ! use already calculated values an n point
          E1=E2
          E1a=E2a
          E1b=E2b

          ! bulk energy at the n+1 point; use all parameters from the
          ! n point
          call text1r_ebulk(text_abnd, text_bbnd, text_bm(n), &
               text_vz(n), text_vr(n), text_vf(n), &
               text_lz(n), text_lr(n), text_lf(n), &
               text_w(n), E2, E2a, E2b)

          ! interpolate bulk energy and derivatives to rp,rm points
          ei = ei + rp*(sp*E2 + sm*E1)*0.5D0*dx
          eai(n) = eai(n) + E1a*sm*dx*rp*0.5D0
          ebi(n) = ebi(n) + E1b*sm*dx*rp*0.5D0

          ei = ei + rm*(sm*E2 + sp*E1)*0.5D0*dx
          eai(n) = eai(n) + E1a*sp*dx*rm*0.5D0
          ebi(n) = ebi(n) + E1b*sp*dx*rm*0.5D0

          ! calculate gradient energy in rp, rm points
          ap = sp*text_abnd+sm*a(n)
          am = sm*text_abnd+sp*a(n)
          bp = sp*text_bbnd+sm*b(n)
          bm = sm*text_bbnd+sp*b(n)
          da = (text_abnd-a(n))/dx
          db = (text_bbnd-b(n))/dx

          call text1r_egrad(rp, ap,bp,da,db, E0,Ea,Eb,Eda,Edb)
          ei = ei + rp*E0*0.5D0*dx
          eai(n) = eai(n) + (Ea*sm*dx - Eda)*rp*0.5D0
          ebi(n) = ebi(n) + (Eb*sm*dx - Edb)*rp*0.5D0

          call text1r_egrad(rm, am,bm,da,db, E0,Ea,Eb,Eda,Edb)
          ei = ei + rm*E0*0.5D0*dx
          eai(n) = eai(n) + (Ea*sp*dx - Eda)*rm*0.5D0
          ebi(n) = ebi(n) + (Eb*sp*dx - Edb)*rm*0.5D0
        else
          ! surface energy
          call text1r_esurf(a(n),b(n),E0,Ea,Eb)
          ei = ei + E0
          eai(n) = eai(n) + Ea
          ebi(n) = ebi(n) + Eb
        endif

      end


! Self test for eint derivatives
      subroutine text1r_eint_selfcheck(n, a,b, d, e)
        implicit none
        integer n,i
        real*8 a(n),b(n),d,e
        real*8 E1,Ea1(n),Eb1(n)
        real*8 E2,Ea2(n),Eb2(n)
        real*8 der1,der2

        call text1r_eint(n,a,b,e1,ea1,eb1)
        do i=2,n

          b(i) = b(i) + d
          call text1r_eint(n,a,b,e2,ea2,eb2)
          b(i) = b(i) - d
          der1=(e2-e1)/d
          der2=(eb2(i)+eb1(i))/2D0
          if ( dabs( der1/der2 - 1D0 ) > e ) then
            write (*,*) 'text1r_eint_selfcheck failed for derivative dE/dBi:'
            write (*,*) 'i: ', i, ' (e2-e1)/db: ', der1, ' ea1*db: ', der2
          endif

          a(i) = a(i) + d
          call text1r_eint(n,a,b,e2,ea2,eb2)
          a(i) = a(i) - d
          der1=(e2-e1)/d
          der2=(ea2(i)+ea1(i))/2D0
          if ( dabs( der1/der2 - 1D0 ) > e ) then
            write (*,*) 'text1r_eint_selfcheck failed for derivative dE/dAi:'
            write (*,*) 'i: ', i, ' (e2-e1)/da: ', der1, ' ea1: ', der2
          endif
        enddo
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Calculate He3 bulk energy E and derivatives dE/da, dE/db
!! in 1d radial coordinated as a function of r, n-vector angles a and b,
!! Energy is divided by (a r H^2).
!! parameters used:
!!   chia*(nub/nu0)^2, for non-zero apsi
!!   lo for non-zero rotation
!!   vd for non-zero flow
!!   de and xi
      subroutine text1r_ebulk(a,b, bm, vz,vr,vf, lz,lr,lf, w,&
                          E,Ea,Eb)
        implicit none
        include 'text1r.fh'
        real*8 a,b,E,Ea,Eb
        real*8 bm, vz,vr,vf, lz,lr,lf, w
        real*8 nz,nr,nf, rzz,rzr,rzf
        real*8 sin_a, sin_b, cos_a, cos_b, cos2b, sin2b
        real*8 con1, con2, help, c,s

        real*8 s3,s5,vd

        vd  = dsqrt(0.4D0 * text_a / text_lhv)

        s3 = sqrt(3D0)
        s5 = sqrt(5D0)

        cos_a = cos(a)
        sin_a = sin(a)
        cos_b = cos(b)
        sin_b = sin(b)
        cos2b = cos(2D0*b)
        sin2b = sin(2D0*b)

        nr=-sin_b*cos_a
        nf=sin_b*sin_a
        nz=cos_b

        c=-0.25D0 !\cos\theta
        s=SQRT(15D0)/4D0 !\sin\theta

        ! H*Rij = cos(t) H + (1-cos(t))(H n) n - sin(t)(H x n)
        rzr=(1D0-c)*nz*nr-s*nf ! H*Rij/abs(H)
        rzf=(1D0-c)*nz*nf+s*nr
        rzz=c+(1D0-c)*nz**2

        E = 0D0
        Ea = 0D0
        Eb = 0D0

        ! magnetic free energy F_DH = -a * (n H)^2
        E = E + sin_b**2
        Eb = Eb + sin2b

        ! spin-orbit free energy
        help = 15D0 * text_ld / text_a / text_h**2  * sin(bm/2D0)**2
        E = E + help * sin_b**2
        Eb = Eb + help * sin2b

        ! flow free energy F_HV
        ! v_d = sqrt(2/5 a/lhv)
        E = E - 2D0/5D0 / (vd**2) * (rzr*vr+rzf*vf+rzz*vz)**2
        help = vr*(-(1D0-c)*cos2b*cos_a - s*cos_b*sin_a) &
             + vf*( (1D0-c)*cos2b*sin_a - s*cos_b*cos_a) &
             + vz*(-(1D0-c)*sin2b)
        Eb = Eb - 4D0/5D0 / (vd**2) * help * (rzr*vr+rzf*vf+rzz*vz)

        help = vr*((1D0-c)*sin_b*cos_b*sin_a - s*sin_b*cos_a) &
             + vf*((1D0-c)*sin_b*cos_b*cos_a + s*sin_b*sin_a)
        Ea = Ea - 4D0/5D0 / (vd**2) * help *(rzr*vr+rzf*vf+rzz*vz)

        ! vortex free energy F_LH
        E = E + text_lo * w/5D0 * (rzr*lr+rzf*lf+rzz*lz)**2

        help = lr*(-(1D0-c)*cos2b*cos_a - s*cos_b*sin_a) &
             + lf*( (1D0-c)*cos2b*sin_a - s*cos_b*cos_a) &
             + lz*(-(1D0-c)*sin2b)
        Eb = Eb + text_lo*w*2D0/5D0 * help * (rzr*lr+rzf*lf+rzz*lz)
        help = lr*((1D0-c)*sin_b*cos_b*sin_a - s*sin_b*cos_a) &
             + lf*((1D0-c)*sin_b*cos_b*cos_a + s*sin_b*sin_a)
        Ea = Ea + text_lo*w*2D0/5D0 * help * (rzr*lr+rzf*lf+rzz*lz)
      end

! Self test for ebulk derivatives
      subroutine text1r_ebulk_selfcheck(a,b, bm, vz,vr,vf, lz,lr,lf, w, d, e)
        implicit none
        real*8 a,b,d,e
        real*8 E1,Ea1,Eb1
        real*8 E2,Ea2,Eb2
        real*8 der1,der2
        real*8 bm, vz,vr,vf, lz,lr,lf, w

        call text1r_ebulk(a,b, bm, vz,vr,vf, lz,lr,lf, w, E1,Ea1,Eb1)
        call text1r_ebulk(a+d,b, bm, vz,vr,vf, lz,lr,lf, w, E2,Ea2,Eb2)
        der1=(E2-E1)/d
        der2=(Ea2+Ea1)/2D0
        if ( dabs( der1/der2 - 1D0 ) > e ) then
          write(*,*) 'text1r_ebulk_selfcheck failed for dE/da:'
          write(*,*) ' a: ', a, ' b: ', b, ' (E2-E1)/da: ', der1, ' Ea: ', der2
        endif
        call text1r_ebulk(a,b+d, bm, vz,vr,vf, lz,lr,lf, w, E2,Ea2,Eb2)
        der1=(E2-E1)/d
        der2=(Eb2+Eb1)/2D0
        if ( dabs( der1/der2 - 1D0 ) > e ) then
          write(*,*) 'text1r_ebulk_selfcheck failed for dE/db:'
          write(*,*) ' a: ', a, ' b: ', b, ' (E2-E1)/db: ', der1, ' Eb: ', der2
        endif
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Calculate He3 gradient energy E and derivatives dE/da, dE/db, dE/da', dE/db'
!! in 1d radial coordinated as a function of r, n-vector angles a and b,
!! and derivatives a'=da/dr, b'=db/dr.
!! Energy is divided by (a H^2).
      subroutine text1r_egrad(r,a,b,da,db, E,Ea,Eb,Eda,Edb)
        implicit none
        include 'text1r.fh'
        real*8 r,a,b,da,db,E,Ea,Eb,Eda,Edb
        real*8 sin_a, sin_b, cos_a, cos_b
        real*8 con1, con2, con3, help
        real*8 s3,s5

        s3 = sqrt(3D0)
        s5 = sqrt(5D0)

        cos_a = cos(a)
        sin_a = sin(a)
        cos_b = cos(b)
        sin_b = sin(b)

        E = 0D0
        Ea = 0D0
        Eb = 0D0
        Eda = 0D0
        Edb = 0D0

        con1 = 5D0*(text_lg2 + text_lg1/2D0) &
                /text_a/text_H**2/text_r**2

        E = E + con1*(db**2 + (sin_b**2)*da**2 + (sin_b**2)/r**2) ! (\nabla n)^2 (?)
        Eda = Eda + con1*2D0*da*sin_b**2
        Edb = Edb + con1*2D0*db
        Eb = Eb + con1 * 2D0*sin_b*cos_b*(da**2 + 1D0/r**2)

        con2 = - 5D0/16D0*text_lg1 &
                 /text_a/text_H**2/text_r**2

        help=(s5*sin_a-s3*cos_b*cos_a)*db + &
             (s5*cos_b*cos_a+s3*sin_a)*sin_b*da + &
             (s5*cos_b*sin_a-s3*cos_a)*sin_b/r  ! s3 \div n + s5 n \rot n (?)
        E = E + con2 * help**2

        Eda = Eda + 2D0*con2*help * (s5*cos_b*cos_a + s3*sin_a)*sin_b
        Edb = Edb + 2D0*con2*help * (s5*sin_a - s3*cos_b*cos_a)

        Ea = Ea + 2D0*con2*help* &
          ( (s5*cos_a + s3*cos_b*sin_a)*db &
          - (s5*cos_b*sin_a - s3*cos_a)*sin_b*da &
          + (s5*cos_b*cos_a + s3*sin_a)*sin_b/r)
        Eb = Eb + 2D0*con2*help* &
          ( (s3*db - s5*sin_b*da)*sin_b*cos_a &
          + (s5*cos_b*cos_a + s3*sin_a)*cos_b*da &
          + (s5*cos_b*sin_a - s3*cos_a)*cos_b/r &
          - s5*sin_b*sin_a*sin_b/r)

        con3 = 5D0*text_lg1 &
               /text_a/text_H**2/text_r**2

        E   = E   + con3 * cos_b*sin_b*db/r
        Edb = Edb + con3 * cos_b*sin_b/r
        Eb  = Eb  + con3 * (cos_b**2 - sin_b**2)*db/r

      end

! Self test for egrad derivatives
      subroutine text1r_egrad_selfcheck(r,a,b, d, e)
        implicit none
        real*8 r,a,b,d,e,g0
        real*8 E1,Ea1,Eb1,Eda1,Edb1
        real*8 E2,Ea2,Eb2,Eda2,Edb2
        real*8 der1,der2

        g0=1D-3
        call text1r_egrad(r, a,  b,g0,g0,E1,Ea1,Eb1,Eda1,Edb1)
        call text1r_egrad(r, a+d,b,g0,g0,E2,Ea2,Eb2,Eda2,Edb2)
        der1=(E2-E1)/d
        der2=(Ea2+Ea1)/2D0
        if ( dabs( der1/der2 - 1D0 ) > e ) then
          write(*,*) 'text1r_egrad_selfcheck failed for dE/da:'
          write(*,*) ' a: ', a, ' b: ', b, ' (E2-E1)/da: ', der1, ' Ea: ', der2
        endif

        call text1r_egrad(r, a,b+d,g0,g0,E2,Ea2,Eb2,Eda2,Edb2)
        der1=(E2-E1)/d
        der2=(Eb2+Eb1)/2D0
        if ( dabs( der1/der2 - 1D0 ) > e ) then
          write(*,*) 'text1r_egrad_selfcheck failed for dE/db:'
          write(*,*) ' a: ', a, ' b: ', b, ' (E2-E1)/db: ', der1, ' Eb: ', der2
        endif

        call text1r_egrad(r, a,b,g0+d,g0,E2,Ea2,Eb2,Eda2,Edb2)
        der1=(E2-E1)/d
        der2=(Eda2+Eda1)/2D0
        if ( dabs( der1/der2 - 1D0 ) > e ) then
          write(*,*) 'text1r_egrad_selfcheck failed for dE/d(da):'
          write(*,*) ' a: ', a, ' b: ', b, ' (E2-E1)/d(da): ', der1, ' Eda: ', der2
        endif

        call text1r_egrad(r, a,b,g0,g0+d,E2,Ea2,Eb2,Eda2,Edb2)
        der1=(E2-E1)/d
        der2=(Edb2+Edb1)/2D0
        if ( dabs( der1/der2 - 1D0 ) > e ) then
          write(*,*) 'text1r_egrad_selfcheck failed for dE/d(db):'
          write(*,*) ' a: ', a, ' b: ', b, ' (E2-E1)/d(db): ', der1, ' Edb: ', der2
        endif
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
        sin2b = sin(2D0*b)
        cos2b = cos(2D0*b)
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

! Self test for esurf derivatives
      subroutine text1r_esurf_selfcheck(a,b, d, e)
        implicit none
        real*8 a,b,d,e
        real*8 E1,Ea1,Eb1
        real*8 E2,Ea2,Eb2
        real*8 der1,der2

        call text1r_esurf(a,b,E1,Ea1,Eb1)
        call text1r_esurf(a+d,b,E2,Ea2,Eb2)
        der1=(E2-E1)/d
        der2=(Ea2+Ea1)/2D0
        if ( dabs( der1/der2 - 1D0 ) > e ) then
          write(*,*) 'text1r_esurf_selfcheck failed for dE/da:'
          write(*,*) ' a: ', a, ' b: ', b, ' (E2-E1)/da: ', der1, ' Ea: ', der2
        endif
        call text1r_esurf(a,b+d,E2,Ea2,Eb2)
        der1=(E2-E1)/d
        der2=(Eb2+Eb1)/2D0
        if ( dabs( der1/der2 - 1D0 ) > e ) then
          write(*,*) 'text1r_esurf_selfcheck failed for dE/db:'
          write(*,*) ' a: ', a, ' b: ', b, ' (E2-E1)/db: ', der1, ' Eb: ', der2
        endif
      end

