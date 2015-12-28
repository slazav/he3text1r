function test_m()
  addpath ../matlab;

  ttc = 0.1;     % Temperature, T/Tc
  p   = 15;      % Pressure, bar
  nu0 = 830000;  % Larmor frequency, Hz
  r   = 0.3;     % Conteiner raduis, cm
  itype   = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n       = 100; % Number of points
  msglev  = -3;  % Message level: -3 silent ...
  omega   = 0;   % Rotation velocity, rad/s
  omega_v = 0; % Rotation velocity of vortex cluster, rad/s
  kr      = 0;   % Parameter for a twisted vortex profile

  % Initialize texture calculation:
  dat = text1r_init(ttc, p, nu0, r, n, itype);

  % Set vortex and velocity profile if needed:
%  dat.lo=5; % lambda/omega is not set by default
%  dat = text1r_set_vortex_cluster(dat, omega, omega_v);
% dat = text1r_set_vortex_uniform_(dat, omega, omega_v);
% dat = text1r_set_vortex_twisted_(dat, omega, kr);

   figure; hold on;

  % Do minimization:
  l=dat.lsg;

  dat.d=dat.d*0.4;
  for dat.lsg = [l 0 -l]
    dat = text1r_minimize(dat);

    % Plot result:
    plot(dat.rr, 180/pi*dat.an, 'r-');
    plot(dat.rr, 180/pi*dat.bn, 'b-');

    plot(dat.rr, 180/pi*dat.al, 'm-');
    plot(dat.rr, 180/pi*dat.bl, 'c-');
  end



  xlim([0, r]);
  legend('\alpha_n', '\beta_n', '\alpha_l', '\beta_l',
         'location', 'northwest')
%  ylim([0 90])
  print -dfig -color a.fig
end
