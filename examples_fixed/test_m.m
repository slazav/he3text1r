function test_m()
  addpath ../matlab;

  ttc = 0.4;     % Temperature, T/Tc
  p   = 15;      % Pressure, bar
  nu0 = 826320;  % Larmor frequency, Hz
  r   = 0.3;     % Conteiner raduis, cm
  itype   = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n       = 100; % Number of points
  msglev  = -3;  % Message level: -3 silent ...
  omega   = 1;   % Rotation velocity, rad/s
  omega_v = 0.2; % Rotation velocity of vortex cluster, rad/s
  kr      = 1;   % Parameter for a twisted vortex profile

  % Initialize texture calculation:

  % Plot result:
  figure; clf; hold on;

  a = 60;
  b = 63.4;
    dat = text1r_init(ttc, p, nu0, r, n, itype);
    dat.bnd=1;
    dat.an=pi+dat.an;
#    dat.bn=pi/2+dat.bn;
    dat.abnd=pi + a/180*pi;
#    dat.bbnd=pi/2+ b/180*pi;

    % Do minimization:
    dat = text1r_minimize(dat);

    plot(dat.rr, 180/pi*dat.an, 'r-');
    plot(dat.rr, 180/pi*dat.bn, 'b-');

    plot(dat.rr, 180/pi*dat.al, 'm-');
    plot(dat.rr, 180/pi*dat.bl, 'g-');

    dat = text1r_init(ttc, p, nu0, r, n, itype);
    dat.bnd=1;
    dat.abnd=a/180*pi;
    dat.bbnd=b/180*pi;

    % Do minimization:
    dat = text1r_minimize(dat);

    plot(dat.rr, 180/pi*dat.an, 'r-');
    plot(dat.rr, 180/pi*dat.bn, 'b-');

    plot(dat.rr, 180/pi*dat.al, 'm-');
    plot(dat.rr, 180/pi*dat.bl, 'g-');



  xlim([0, r]);
end
