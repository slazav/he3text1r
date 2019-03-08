function test_m2()
  addpath ../matlab;
  addpath ../../he3lib/matlab;

  ttc = 0;       % Temperature, T/Tc
  p   = 29;      % Pressure, bar
  nu0 = 833000;  % Larmor frequency, Hz
  r   = 0.3;     % Conteiner raduis, cm
  itype   = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n       = 180; % Number of points
  msglev  = -3;  % Message level: -3 silent ...

  figure; hold on;

  % Initialize texture calculation:
  dat = text1r_init(ttc, p, nu0, r, n, itype);
  dat = text1r_minimize(dat);

  cper = he3_cperp(ttc, p);
  cpar = he3_cpar(ttc, p);
  nub  = he3_nu_b(ttc, p);

  mr = 1.652; % G/A/cm^2
  Imin = 4;
  f2r= -mr*Imin * he3_gyro/2/pi * dat.rr.^2;

  df0 = text1r_qc_states(dat, cper, cpar, nu0, 0,   nub, -2, 300, 200, 3);
  df1 = text1r_qc_states(dat, cper, cpar, nu0, 0,   nub,  1, 300, 200, 3);
  df2 = text1r_qc_states(dat, cper, cpar, nu0, f2r, nub,  1, 300, 200, 3);
  df3 = text1r_qc_states(dat, cper, cpar, nu0, 0,   nub,  2, 300, 200, 3);
  df4 = text1r_qc_states(dat, cper, cpar, nu0, 0,   nub,  3, 300, 200, 3);
  [df5 rr5 psi5] = text1r_wv_states(dat, cper, cpar, nu0, 0, nub, -1, 800, 1000, 10);
  [df6 rr6 psi6] = text1r_wv_states(dat, cper, cpar, nu0, 0, nub,  1, 800, 1000, 10);
  [df7 rr7 psi7] = text1r_wv_states(dat, cper, cpar, nu0, 0, nub, -2, 800, 1000, 10);
  [df8 rr8 psi8] = text1r_wv_states(dat, cper, cpar, nu0, 0, nub,  2, 800, 1000, 10);

  plot(1:length(df0)-1, diff(df0), 'ks-');
  plot(1:length(df1)-1, diff(df1), 'g.-');
  plot(1:length(df2)-1, diff(df2), 'c.-');
  plot(1:length(df3)-1, diff(df3), 'k.-');
  plot(1:length(df4)-1, diff(df4), 'm.-');
  plot(1:length(df5)-1, diff(df5), 'rs-');
  plot(1:length(df6)-1, diff(df6), 'r.-');
  plot(1:length(df7)-1, diff(df7), 'bs-');
  plot(1:length(df8)-1, diff(df8), 'b.-');
  legend('qc, f=const, quadratic potential',...
         'qc, H=const',...
         'qc, H=const, Imin = 4A',...
         'qc, f=const',...
         'qc, f=const, adjust texture',...
         'wave, quadratic potential, f(R)=0',...
         'wave, f=const, f(R)=0',...
         'wave, quadratic potential, df(R)=0',...
         'wave, f=const, df(R)=0',...
         'Location', 'northwest');

  figure; hold on;
  for i = 1:length(psi6(1,:))
    plot(rr6, psi6(:,i)+i*0.1);
  end
end
