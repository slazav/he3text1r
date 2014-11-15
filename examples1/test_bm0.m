function test_bm0()
% Calculate self-consistent texture and precessing magnetization
% distibution for a given beta_m(r=0) value.

  addpath ../matlab;

  ttc = 0.2;     % Temperature, T/Tc
  p   = 10;      % Pressure, bar
  nu0 = 826320;  % Larmor frequency, Hz
  r   = 0.3;     % Conteiner raduis, cm
  itype = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n     = 100; % Number of points

  bm0 = 10; % degrees

  % Initialize texture calculation:
  dat = text1r_init(ttc, p, nu0, r, n, itype);

  % Calculate ground state wavefunction Psi (no normalization)
  dat = text1r_qball_bm0(dat, pi/180*bm0);

  % Plot result:
  figure; hold on;
  plot(dat.rr, 180/pi*dat.an, 'r-');
  plot(dat.rr, 180/pi*dat.bn, 'b-');

  plot(dat.rr, dat.bm*180/pi, 'g-');
  xlim([0, r]);
  legend('beta_N, deg', 'alpha_N, deg', 'beta_M, deg',...
    'location', 'southeast');
  title(sprintf(...
    'texture for qball with beta_M=%.1f deg, ttc=%.2f, P=%.2f bar, f=%.1f Hz',...
    bm0, ttc,p,nu0))
end
