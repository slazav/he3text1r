function test_qc_states()
  addpath ../matlab;
  addpath ../../he3lib/matlab;

  ttc = 0;       % Temperature, T/Tc
  p   = 29;      % Pressure, bar
  nu0 = 833000;  % Larmor frequency, Hz
  r   = 0.3;     % Conteiner raduis, cm
  itype   = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n       = 200; % Number of points in the texture calculation
  msglev  = -3;  % Message level: -3 silent ...
  nq      = 1000 % Number of points in quasiclassical calculation
  nw      = 2000 % Number of points in wave calculation

  figure; hold on;

  % Initialize texture calculation:
  dat = text1r_init(ttc, p, nu0, r, n, itype);
  dat = text1r_minimize(dat);

  print_text(sprintf('calc_text_%.0f.txt', p), dat);

  cper = he3_cperp(ttc, p);
  cpar = he3_cpar(ttc, p);
  nub  = he3_nu_b(ttc, p);

  mr = 1.652; % G/A/cm^2
  Imin = 4;
  f2r= -mr*Imin * he3_gyro/2/pi * dat.rr.^2;

  fmax = nub^2/2/nu0*sin(dat.bmax)^2;

  % linear texture with beta0
  dfql = text1r_qc_states(dat, cper, cpar, nu0, 0,   nub, -2, nq, 1000, 10);
  [df1l rr1l psi1l] = text1r_wv_states(dat, cper, cpar, nu0, 0, nub, -1, nw, 1000, 10);
  [df2l rr2l psi2l] = text1r_wv_states(dat, cper, cpar, nu0, 0, nub, -2, nw, 1000, 10);

  plot(1:length(dfql)-1, diff(dfql), 'ks-');
  plot(1:length(df1l)-1, diff(df1l), 'rs-');
  plot(1:length(df2l)-1, diff(df2l), 'bs-');

  % constant field
  dfqh = text1r_qc_states(dat, cper, cpar, nu0, 0,   nub,  1, 300, 200, 3);
  nqh = interp1(dfqh,0:length(dfqh)-1, fmax);
  plot(1:length(dfqh)-1, diff(dfqh), 'k*-');
  print_states(sprintf('calc_states_qh_%.0f.txt', p), dfqh, fmax, nqh);

  % constant frequency
  dfqf = text1r_qc_states(dat, cper, cpar, nu0, 0,   nub,  2, nq, 1000, 10);
  [df1f rr1f psi1f] = text1r_wv_states(dat, cper, cpar, nu0, 0, nub,  1, nw, 1000, 10);
  [df2f rr2f psi2f] = text1r_wv_states(dat, cper, cpar, nu0, 0, nub,  2, nw, 1000, 10);

  print_modes("modes_psi1f", rr1f, psi1f);
  print_modes("modes_psi2f", rr2f, psi2f);

  nqf = interp1(dfqf,0:length(dfqf)-1, fmax);
  n1f = interp1(df1f,0:length(df1f)-1, fmax);
  n2f = interp1(df2f,0:length(df2f)-1, fmax);
  if 1
    plot(1:length(dfqf)-1, diff(dfqf), 'k.--');
    plot(1:length(df1f)-1, diff(df1f), 'r.--');
    plot(1:length(df2f)-1, diff(df2f), 'b.--');
    plot(nqf*[1 1], [550 750], 'k--')
    plot(n1f*[1 1], [550 750], 'r--')
    plot(n2f*[1 1], [550 750], 'b--')
  end
  print_states(sprintf('calc_states_qf_%.0f.txt', p), dfqf, fmax, nqf);
  print_states(sprintf('calc_states_1f_%.0f.txt', p), df1f, fmax, n1f);
  print_states(sprintf('calc_states_2f_%.0f.txt', p), df2f, fmax, n2f);

  % adaptive texture
  % we use quasiclassic states for reference
  if 1
  mm = max([length(dfqf), length(df1f), length(df2f)]);
  for i=1:mm
    if i<=length(dfqf)
      dat.H = 2*pi*(nu0-dfqf(i))/he3_gyro;
      dat = text1r_minimize(dat);
      tmp1 = text1r_qc_states(dat, cper, cpar, nu0, 0,   nub,  2, nq, 1000, 10);
      tmp2 = text1r_wv_states(dat, cper, cpar, nu0, 0, nub,  1, nw, 1000, 10);
      tmp3 = text1r_wv_states(dat, cper, cpar, nu0, 0, nub,  2, nw, 1000, 10);
    end
    if i<=length(tmp1); dfqfa(i) = tmp1(i); end
    if i<=length(tmp2); df1fa(i) = tmp2(i); end
    if i<=length(tmp3); df2fa(i) = tmp3(i); end
    fprintf('%d/%d\n', i,mm);
  end
  plot(1:length(dfqfa)-1, diff(dfqfa), 'k.-');
  plot(1:length(df1fa)-1, diff(df1fa), 'r.-');
  plot(1:length(df2fa)-1, diff(df2fa), 'b.-');
  nqf = interp1(dfqfa,0:length(dfqfa)-1, fmax);
  n1f = interp1(df1fa,0:length(df1fa)-1, fmax);
  n2f = interp1(df2fa,0:length(df2fa)-1, fmax);
  plot(nqf*[1 1], [550 750], 'k-')
  plot(n1f*[1 1], [550 750], 'r-')
  plot(n2f*[1 1], [550 750], 'b-')
  print_states(sprintf('calc_states_qfa_%.0f.txt', p), dfqfa, fmax, nqf);
  print_states(sprintf('calc_states_1fa_%.0f.txt', p), df1fa, fmax, n1f);
  print_states(sprintf('calc_states_2fa_%.0f.txt', p), df2fa, fmax, n2f);
  end

  legend(...
         'qc, f=const, quadratic potential',...
         'wave, quadratic potential, f(R)=0',...
         'wave, quadratic potential, df(R)=0',...
         'qc, H=const',...
         'qc, f=const',...
         'wave, f=const, f(R)=0',...
         'wave, f=const, df(R)=0',...
         'qc, f=const, adaptive texture',...
         'wave, f=const, f(R)=0, adaptive texture',...
         'wave, f=const, df(R)=0, adaptive texture',...
         'Location', 'northwest');

  figure; hold on;
  for i = 1:length(psi1f(1,:))
    plot(rr1f, psi1f(:,i)+i*0.1);
  end
end

function print_states(fname, df, fmax, nmax)
  fo = fopen(fname, 'w');
  fprintf(fo, '# N freq(Hz)\n');
  fprintf(fo, '%.2f %.2f\n', nmax, fmax);
  fprintf(fo, '\n');
  for i =1:length(df)
    fprintf(fo, '%d %.2f\n', i-1, df(i));
  end
  fclose(fo);
end

function print_text(fname, dat)
  fo = fopen(fname, 'w');
  fprintf(fo, '# r an bn\n');
  for i =1:length(dat.rr)
    fprintf(fo, '%f %f %f\n', dat.rr(i), dat.an(i), dat.bn(i));
  end
  fclose(fo);
end

function print_modes(fname, rr, psi)
  fo = fopen(fname, 'w');
  fprintf(fo, '# r an bn\n');
  for j = 1:length(psi(1,:))
    for i =1:length(rr)
      fprintf(fo, '%f %f\n', rr(i), psi(i,j));
    end
    fprintf(fo, '\n');
  end
  fclose(fo);
end
