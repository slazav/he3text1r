      program main
        include '../text1r.fh'
        real*8 ttc=0.6D0
        call text1r_init(ttc, 15D0, 250D0, 0.3D0, 0)
        call text1r_minimize(0)
        call text1r_print('result.dat')
      end
