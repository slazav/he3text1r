#!/usr/bin/gnuplot

set style data lines
plot\
  "result.dat" using 1:2 title "alpha_n",\
  "result.dat" using 1:3 title "beta_n",\
  0
pause -1

plot\
  "result.dat" using 1:5 title "v_r",\
  "result.dat" using 1:6 title "v_z",\
  "result.dat" using 1:7 title "v_f",\
  "result.dat" using 1:8 title "l_r",\
  "result.dat" using 1:9 title "l_z",\
  "result.dat" using 1:10 title "l_f",\
  "result.dat" using 1:11 title "w",\
  0
pause -1
