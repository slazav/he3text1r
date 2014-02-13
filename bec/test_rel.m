function test_rel()
  addpath ../matlab;
  addpath ../../lib/matlab;

  ttc = 0.13;     % Temperature, T/Tc
  p   = 0;      % Pressure, bar
  nu0 = 826320;  % Larmor frequency, Hz
  r   = 0.3;     % Conteiner raduis, cm
  itype   = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n       = 100; % Number of points
  msglev  = -3;  % Message level: -3 silent ...


  % Initialize texture calculation:
  dat = text1r_init(ttc, p, nu0, r, n, itype);

  e(1)=0.0005;
  tau=10;
  dt=1;
  time=1:dt:35;
  a=0.0001;

  lll=1;

  if lll
    load test_rel.dat
  else
    figure; hold on;
    for i=1:length(time)
      [dat fre(i)] = text1r_qball_e(dat, e(i));
      % gradient of bm

      amp(i)=trapz(dat.rr, 2*pi*dat.rr.*sin(abs(dat.bm)));
      grd(i)=text1r_qball_gm(dat)
      bm0(i)=abs(dat.bm(1));
      plot(dat.rr, 180/pi*dat.an, 'r.-');
      plot(dat.rr, 180/pi*dat.bn, 'b.-');
      plot(dat.rr, 180/pi*abs(dat.bm), 'g.-');

      if i<length(time);
        e(i+1) = e(i) - dt * a * grd(i);
        e(i+1)
      end
    end
    save test_rel.dat time e fre amp bm0 grd
  end

  nub=he3_nu_b(ttc,p);
  fre = fre * nub^2/(2*nu0);
%  plot(time, e, 'r.-');
%  plot(time, fre, 'g.-');
%  plot(time, amp, 'b.-');

  figure;
  subplot(2,2,1); hold on;
  plot(time, fre, 'r.-');
  title('freq, Hz vs time, s');

  subplot(2,2,2); hold on;
  plot(time, amp, 'r.-');
  title('amp = int(\sin\beta) vs time, s');

  subplot(2,2,3); hold on;
  plot(fre, amp, 'b.-');
  title('amp = int(\sin\beta) vs freq, Hz');

  subplot(2,2,4); hold on;
%  plot(time, grd, 'g.-');
%  plot(e, amp.^2 .* sqrt(fre), 'g.-');
%  plot(time, log(e), 'b.-');
  plot(time, log(amp), 'r.-');
%  plot(time, log(bm0)-5, 'g.-');

end
