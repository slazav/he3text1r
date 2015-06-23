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
  h(1)=subplot(1,2,1); hold on;  title('beta_N angle vs r and temperature (+ y shift)');
  h(2)=subplot(1,2,2); hold on;  title('difference');
  for i=1:length(ttc);
    % Initialize texture calculation:
    dat = text1r_init(ttc(i), p, 1000*f0, r, n, itype);

    % Do minimization:
    dat  = text1r_minimize(dat);

    % three parameters
    bn1 = bn3par_poly(dat.rr, dat.db0, dat.db1, dat.bmax);
    bn2 =  bn3par_exp(dat.rr, dat.db0, dat.db1, dat.bmax);

    sh = (i-1)/5;
    plot(h(1), dat.rr, dat.bn + sh, 'r.');
    plot(h(1), dat.rr, bn1 + sh, 'b-');
    plot(h(1), dat.rr, bn2 + sh, 'g-');

    text(h(1), 0.005, sh+0.1, sprintf('%.2f Tc', ttc(i)));

    plot(h(2), dat.rr, (bn1-dat.bn)*10, 'b-');
    plot(h(2), dat.rr, (bn2-dat.bn)*10, 'g-');
  end
  legend(h(1), 'calculation',...
         'poly 3-parametar approximation',...
         'exp 3-parametar approximation',...
         'location', 'northwest')
  xlim([0 dat.R])
end

function bl = bn2bl(bn)
  bl = 2*asin( sqrt(5/8) * sin(bn) );
end

function bn = bl2bn(bl)
  bn = asin( sqrt(1.6) * sin(bl/2) );
end

function bn = bn3par_poly(rr, db0,db1,bm)
  R = rr(end);

  % 4-degree polinomial with d2/dr2 = 0 at r=0
  B = (4*(bm-db0*R) - (db1-db0)*R) / R^3;
  A = ((db1-db0)*R - 3*(bm-db0*R)) / R^4;
  bn = A*rr.^4 + B*rr.^3 + db0*rr;
end

function bn = bn3par_exp(rr, db0,db1,bm)
  R = rr(end);
  % linear + exponent
  xi = R^2 * (db1-db0)/(bm-db0*R);
  A  = (bm-db0*R)*exp(xi/R);
  bn = db0*rr +  A*exp(-xi./rr);
end
