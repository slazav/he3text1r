function test_m()
  addpath ../matlab;

  ttc = 0.13;     % Temperature, T/Tc
  p   = 0;      % Pressure, bar
  r   = 0.3;     % Conteiner raduis, cm
  itype   = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n       = 100; % Number of points
  msglev  = -3;  % Message level: -3 silent ...

  % Initialize texture calculation:

  figure; hold on;

  nu0=250000:10000:1500000;
  for i=1:length(nu0)
    dat = text1r_init(ttc, p, nu0(i), r, n, itype);
    % Do minimization:
    dat = text1r_minimize(dat);
    db(i)=(dat.bn(2)-dat.bn(1))/(dat.rr(2)-dat.rr(1));
    plot(dat.rr, 180/pi*dat.an, 'r-');
    plot(dat.rr, 180/pi*dat.bn, 'b-');
  end
  save db.txt nu0 db dat

  figure; hold on;
  plot(nu0/1000, db, 'r.-');

end
