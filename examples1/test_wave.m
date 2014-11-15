function test_wave()
% Calculate magnon wavefunction in a textural potential (no normalization)
  addpath ../matlab;

  ttc = 0.4;     % Temperature, T/Tc
  p   = 15;      % Pressure, bar
  nu0 = 826320;  % Larmor frequency, Hz
  r   = 0.3;     % Conteiner raduis, cm
  itype   = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n       = 100; % Number of points

  % Initialize texture calculation:
  dat = text1r_init(ttc, p, nu0, r, n, itype);

  % Do minimization:
  dat  = text1r_minimize(dat);

  % Calculate wavefunction Psi (no normalization)
  [psi, e0] = text1r_wave(dat, [0 2]);

  % Plot result:
  figure; hold on;
  plot(dat.rr, 180/pi*dat.an, 'r-');
  plot(dat.rr, 180/pi*dat.bn, 'b-');

  plot(dat.rr, psi/psi(1) * 60, 'g-');
  xlim([0, r]);
end
