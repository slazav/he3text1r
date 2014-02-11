function test_m()
  addpath ../matlab;

  ttc = 0.99;     % Temperature, T/Tc
  p   = 15;      % Pressure, bar
  nu0 = 926320;  % Larmor frequency, Hz
  r   = 0.3;     % Conteiner raduis, cm
  itype   = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n       = 100; % Number of points
  msglev  = -3;  % Message level: -3 silent ...
  omega   = 1;   % Rotation velocity, rad/s
  omega_v = 0.2; % Rotation velocity of vortex cluster, rad/s
  kr      = 1;   % Parameter for a twisted vortex profile

  % Initialize texture calculation:
  dat = text1r_init(ttc, p, nu0, r, n, itype);

  % Do minimization:
  dat = text1r_selfcheck(dat);
  dat = text1r_minimize(dat);
  dat = text1r_selfcheck(dat);


  if dat.err==3
    ee=0;
    for i=1:100
      for j=dat.n-1:-1:1
        if abs(dat.bn(j))<1e-2; dat.an(j)=dat.an(j+1); end
      end
      dat = text1r_minimize(dat);
      if (abs(ee-dat.energy)<1e-30); break; end
      ee=dat.energy;
    end
  end

  % Plot result:
  figure; hold on;
  plot(dat.rr, 180/pi*dat.an, 'r-');
  plot(dat.rr, 180/pi*dat.bn, 'b-');
end
