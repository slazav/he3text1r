function test_m()
  addpath ../matlab;

  ttc = 0.4;     % Temperature, T/Tc
  p   = 15;      % Pressure, bar
  nu0 = 826320;  % Larmor frequency, Hz
  r   = 0.3;     % Conteiner raduis, cm
  itype   = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n       = 100; % Number of points
  msglev  = -3;  % Message level: -3 silent ...
  omega   = 1;   % Rotation velocity, rad/s
  omega_v = 0.2; % Rotation velocity of vortex cluster, rad/s
  kr      = 1;   % Parameter for a twisted vortex profile

  % Initialize texture calculation:
  dat = text1r_init(ttc, p, nu0, r, n, itype);
dat.H

  % Set vortex and velocity profile if needed:
%  dat.lo=5; % lambda/omega is not set by default
%  dat = text1r_set_vortex_cluster(dat, omega, omega_v);
% dat = text1r_set_vortex_uniform_(dat, omega, omega_v);
% dat = text1r_set_vortex_twisted_(dat, omega, kr);

  % Do minimization:
  dat = text1r_minimize(dat);

  % Plot result:
  figure; hold on;
  plot(dat.rr, 180/pi*dat.an, 'r-');
  plot(dat.rr, 180/pi*dat.bn, 'b-');
  xlim([0, r]);
end
