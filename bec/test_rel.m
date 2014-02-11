function test_rel()
  addpath ../matlab;

  ttc = 0.13;     % Temperature, T/Tc
  p   = 0;      % Pressure, bar
  nu0 = 826320;  % Larmor frequency, Hz
  r   = 0.3;     % Conteiner raduis, cm
  itype   = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n       = 100; % Number of points
  msglev  = -3;  % Message level: -3 silent ...
  omega   = 1;   % Rotation velocity, rad/s
  omega_v = 0.2; % Rotation velocity of vortex cluster, rad/s
  kr      = 1;   % Parameter for a twisted vortex profile


  % Initialize texture calculation:
  dat = text1r_init(ttc, p, nu0, r, n, itype);

  e0=0.001;
  tau=10;
  time=1:50;
  e=e0*exp(-time/tau);

  figure; hold on;
  for i=1:length(time)
    [dat fre(i)] = text1r_qball_e(dat, e(i));
    amp(i)=trapz(dat.rr, 2*pi*dat.rr.*sin(abs(dat.bm)));
    plot(dat.rr, 180/pi*dat.an, 'r.-');
    plot(dat.rr, 180/pi*dat.bn, 'b.-');
    plot(dat.rr, 180/pi*abs(dat.bm), 'g.-');
  end

%  plot(time, e, 'r.-');
%  plot(time, fre, 'g.-');
%  plot(time, amp, 'b.-');

  figure; hold on;
  plot(fre, amp, 'r.-');

  figure; hold on;
  plot(time, amp, 'r.-');
  plot(time, fre, 'b.-');

end
