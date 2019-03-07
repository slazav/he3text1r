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
  f2r= -mr*Imin * he3_gyro/2/pi;

  df0 = text1r_qc_states(dat, cper, cpar, nu0, 0,   nub, -1, 300, 1000, 3);
  df1 = text1r_qc_states(dat, cper, cpar, nu0, 0,   nub,  0, 300, 1000, 3);
  df2 = text1r_qc_states(dat, cper, cpar, nu0, f2r, nub,  0, 300, 1000, 3);
  df3 = text1r_qc_states(dat, cper, cpar, nu0, 0,   nub,  1, 300, 1000, 3);
  df4 = text1r_qc_states(dat, cper, cpar, nu0, 0,   nub,  2, 300, 1000, 3);
  [df5 rr5 psi5] = text1r_wv_states(dat, cper, cpar, nu0, 0, nub, -1, 300, 1000, 10);
  [df6 rr6 psi6] = text1r_wv_states(dat, cper, cpar, nu0, 0, nub,  1, 300, 1000, 10);


  plot(0:length(df0)-2, df0(2:end)-df0(1:end-1), 'k.-');
  plot(0:length(df1)-2, df1(2:end)-df1(1:end-1), 'r.-');
  plot(0:length(df2)-2, df2(2:end)-df2(1:end-1), 'b.-');
  plot(0:length(df3)-2, df3(2:end)-df3(1:end-1), 'g.-');
  plot(0:length(df4)-2, df4(2:end)-df4(1:end-1), 'm.-');
  plot(0:length(df5)-2, df5(2:end)-df5(1:end-1), 'ms-');
  plot(0:length(df6)-2, df6(2:end)-df6(1:end-1), 'gs-');
  legend('qc, quadratic potential',...
         'qc, H=const, Imin = 0',...
         'qc, H=const, Imin = 4A',...
         'qc, f=const',...
         'qc, f=const, adjust texture',...
         'wave, quadratic potential',...
         'wave, f=const',...
         'Location', 'northwest');

  figure; hold on;

  for i = 1:length(psi6(1,:))
    plot(rr6, psi6(:,i)+i*0.1);
  end
end
