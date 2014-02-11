function test_m2()
  addpath ../matlab;

  ttc = 0.13;     % Temperature, T/Tc
  p   = 0;      % Pressure, bar
  nu0 = 826320;  % Larmor frequency, Hz
  r   = 0.3;     % Conteiner raduis, cm
  itype   = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n       = 100; % Number of points
  msglev  = -3;  % Message level: -3 silent ...
  omega   = 1;   % Rotation velocity, rad/s
  omega_v = 0.2; % Rotation velocity of vortex cluster, rad/s
  kr      = 1;   % Parameter for a twisted vortex profile

  figure; hold on;

  % Initialize texture calculation:
  dat = text1r_init(ttc, p, nu0, r, n, itype);

  for e = 0:0.0001:0.001
    dat = text1r_qball_e(dat, e);
    plot(dat.rr, 180/pi*dat.an, 'r.-');
    plot(dat.rr, 180/pi*dat.bn, 'b.-');
    plot(dat.rr, 180/pi*abs(dat.bm), 'g.-');
  end
  xlim([0, r]);
end
