      program main
        include '../text1r.fh'
        real*8 ttc, p, nu0, r, omega, omega_v, kr
        integer n, itype, msglev

        data   ttc /0.4D0/,    ! Temperature, T/Tc
     .         p   /15D0/,     ! Pressure, bar
     .         nu0 /826320D0/, ! Larmor frequency, Hz
     .         r   /0.3D0/,    ! Conteiner raduis, cm
     .         itype  /0/,     ! Initial conditions: 0 - normal, 1 - with 90deg peak...
     .         n      /50/,    ! Number of points
     .         msglev /-3/,    ! Message level: -3 silent ...
     .         omega   /1D0/,  ! Rotation velocity, rad/s
     .         omega_v /0.2D0/,! Rotation velocity of vortex cluster, rad/s
     .         kr      /1D0/   ! Parameter for a twisted vortex profile

!       Initialize texture calculation:
        call text1r_init(ttc, p, nu0, r, n, itype)

!       Set vortex and velocity profile if needed:
!       text_lo=5D0; ! lambda/omega is not set by default
!       call text1r_set_vortex_cluster(omega, omega_v);
!         call text1r_set_vortex_uniform(omega, omega_v);
!         call text1r_set_vortex_twisted(omega, kr);

!       Do minimization:
        call text1r_minimize(msglev)
!       Print results to file:
        call text1r_print('result.dat')
      end
