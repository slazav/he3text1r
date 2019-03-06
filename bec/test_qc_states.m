function test_m2()
  addpath ../matlab;
  addpath ../../he3lib/matlab;

  ttc = 0.13;     % Temperature, T/Tc
  p   = 29;      % Pressure, bar
  nu0 = 826320;  % Larmor frequency, Hz
  r   = 0.3;     % Conteiner raduis, cm
  itype   = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n       = 100; % Number of points
  msglev  = -3;  % Message level: -3 silent ...

  figure; hold on;

  % Initialize texture calculation:
  dat = text1r_init(ttc, p, nu0, r, n, itype);
  dat = text1r_minimize(dat);

  cper = he3_cperp(ttc, p);
  cpar = he3_cpar(ttc, p);
  nub  = he3_nu_b(ttc, p);

  mr = 1.652; % G/A/cm^2
  Imin = 2;
  f2r= -mr*Imin * he3_gyro/2/pi;

  df1 = text1r_qc_states(dat, cper, cpar, nu0, 0, nub, 0, 1000, 3);
  df2 = text1r_qc_states(dat, cper, cpar, nu0, f2r, nub, 0, 1000, 3);
  df3 = text1r_qc_states(dat, cper, cpar, nu0, 0, nub, 1, 1000, 3);
  df4 = text1r_qc_states(dat, cper, cpar, nu0, 0, nub, 2, 1000, 3);

  plot(0:length(df1)-2, df1(2:end)-df1(1:end-1), 'r.-');
  plot(0:length(df2)-2, df2(2:end)-df2(1:end-1), 'b.-');
  plot(0:length(df3)-2, df3(2:end)-df3(1:end-1), 'g.-');
  plot(0:length(df4)-2, df4(2:end)-df4(1:end-1), 'm.-');
  legend('H=const, Imin = 0',...
         'H=const, Imin = 2A',...
         'f=const',...
         'f=const, adjust texture',...
         'Location', 'southeast');
end
