program main
  include '../text1r.fh'
  call text1r_init(0.6D0, 15D0, 250D0, 0.3D0, 0)
  call text1r_minimize(0)
  call text1r_print('result.dat')
end
