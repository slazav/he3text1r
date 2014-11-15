function plot_db0()
  % beta_N derivative in the center vs pressure and field
  addpath ../matlab;

  ttc = 0.1;     % Temperature, T/Tc, low-temperature limit
  r   = 0.3;     % Conteiner raduis, cm
  itype   = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n       = 100; % Number of points

  p=[0 0.5 4 9.2 15 23.1 26.3 29.0];  % values used in our experiment
  f0 = 500:20:900; % kHz

  figure; clf; hold on;
  for np=1:length(p);
    for nf = 1:length(f0)
      % Initialize texture calculation:
      dat = text1r_init(ttc, p(np), 1000*f0(nf), r, n, itype);

      % Do minimization:
      dat  = text1r_minimize(dat);
      db0(nf) = dat.db0;
      db1(nf) = dat.db1;
      bmax(nf) = dat.bmax;
    end
    plot(f0, db0, 'r.-');
%    plot(f0, db1, 'b.-');
%    plot(f0, bmax, 'r.-');
  end
end
