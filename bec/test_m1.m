function test_m1()
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

  ee=0;
  for i=1:100
    dat = text1r_minimize(dat);
    [dat e0] = text1r_wave(dat);
    dat = text1r_minimize(dat);
    if abs(ee-e0) < 1e-6; break; end
    ee=e0;
  end
  dat.apsi
  plot(dat.rr, 180/pi*dat.an, 'r.-');
  plot(dat.rr, 180/pi*dat.bn, 'b.-');
  plot(dat.rr, 180/pi*sqrt(abs(dat.apsi)), 'g.-');
  xlim([0, r]);
end
