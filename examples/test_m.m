function test_m()
  addpath ../matlab;
  dat = text1r_init(0.6, 15, 930000, 0.3, 0);
  dat = text1r_set_vortex_cluster(dat, 1.0, 0.5);
  dat.lo = 5; % lambda/omega
  dat = text1r_minimize(dat);
  figure; hold on;
  plot(dat.rr, 180/pi*dat.an, 'r-');
  plot(dat.rr, 180/pi*dat.bn, 'b-');
end
