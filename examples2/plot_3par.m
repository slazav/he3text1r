function plot_3par()
  % aproximation of the texture by three parameters,
  % db0, db1, bmax
  addpath ../matlab;

  ttc = 0.1;     % Temperature, T/Tc, low-temperature limit
  r   = 0.3;     % Conteiner raduis, cm
  itype   = 0;   % Initial conditions: 0 - normal, 1 - with 90deg peak...
  n       = 100; % Number of points

  p  = 0;   % bar
  f0 = 833; % kHz
  ttc = 0.1:0.1:0.9;

  figure; clf; hold on;
  for i=1:length(ttc);
    % Initialize texture calculation:
    dat = text1r_init(ttc(i), p, 1000*f0, r, n, itype);

    % Do minimization:
    dat  = text1r_minimize(dat);

    % three parameters
    bn1 = bn3par(dat.rr, dat.db0, dat.db1, dat.bmax);

    sh = (i-1)/5;
    plot(dat.rr, dat.bn + sh, 'r.-');
    plot(dat.rr, bn1 + sh, 'b.-');
    text(0.005, sh+0.1, sprintf('%.2f Tc', ttc(i)));
  end
  title('beta_N angle vs r and temperature (+ y shift)')
  legend('calculation',...
         '3-parametar approximation',...
         'location', 'northwest')
  xlim([0 dat.R])
end

function bl = bn2bl(bn)
  bl = 2*asin( sqrt(5/8) * sin(bn) );
end

function bn = bl2bn(bl)
  bn = asin( sqrt(1.6) * sin(bl/2) );
end

function bn = bn3par(rr, db0,db1,bm)
  R = rr(end);

  % 4-degree polinomial with d2/dr2 = 0 at r=0
  B = (4*(bm-db0*R) - (db1-db0)*R) / R^3;
  A = ((db1-db0)*R - 3*(bm-db0*R)) / R^4;
  bn = A*rr.^4 + B*rr.^3 + db0*rr;
end
