% Calculate all magnon states in the texture using quasiclassical approach
%
% dat - Texture parameter structure. Do not forget to calculate texture before.
% cper, cpar - Spin-wave velocities (usually corresponds to gradient energy parameters lg1,lg2 in dat).
% f0 - Larmor field in cell center [Hz] / frequency (usually corresponds to dat.H).
% f2r  - d2(fL)/dr^2 across the cell (negative for usual trap) [Hz/cm^2].
% nu_b - Leggett frequency [Hz].
% mode - 0: constant field (default), 1: constant freq, 2: constant freq, recalculate texture.
% Nmax - Max number of states to return (default 1000).
% Nex  - Max number of non-local states to return (default 0).

function df = text1r_qc_states_hconst(dat, cper, cpar, f0, f2r, nu_b, mode, Nmax, Nex)
  format long;
  if (nargin < 7) mode = 0; end
  if (nargin < 8) Nmax = 1000; end
  if (nargin < 9) Nex  = 0; end

  % find lx,ly,lz
  nz = cos(dat.bn);
  nx = sin(dat.bn).*cos(dat.an);
  ny = sin(dat.bn).*sin(dat.an);

  ct=-1/4; st=sqrt(15)/4;
  lx=(1-ct)*nz.*nx - st*ny;
  ly=(1-ct)*nz.*ny + st*nx;
  lz= ct + (1-ct)*nz.^2;

  % velocity of transverse spin waves propagating
  % along cell radius
  ceff=f_ceff(dat, cpar, cper);

  % dipolar potential, Hz^2
  uD = f_uD(dat, nu_b);

  % For mode=0 f0 is larmor frequency in the center
  % frequency shifts are calculated from this level
  % For mode=1,2 f0 is a frequency.
  w0 = 2*pi*f0;
  wr = 2*pi*f2r*dat.rr.^2;
  H0 = dat.H; % for mode 2

  for N=0:Nmax
    % x - frequency from larmor
    if mode == 0
      % constant field mode
      k2func = @(x) ((w0+x)*(x-wr)-uD)./ceff.^2;
    else
      % constant frequency mode
      k2func = @(x) (w0*(x-wr)-uD)./ceff.^2;
    end
    zfunc1 = @(x) kint(N, dat.rr, k2func(x));

    if N>0 dwp = dw(N); else dwp = 0; end
    dw(N+1) = fzero(zfunc1, dwp);

    % in mode==2 recalculate the texture and update frequency shift
    if mode == 2
       dat.H = H0 - dw(N+1)/he3_gyro;
       dat = text1r_minimize(dat);
       % update uD and ceff:
       ceff=f_ceff(dat, cpar, cper);
       uD = f_uD(dat, nu_b);
       dw(N+1) = fzero(zfunc1, dw(N+1));
    end

    % check localization
    k2 = k2func(dw(N+1));
    if (k2(end) >= 0);
      if (Nex > 0) Nex = Nex-1;
      else break;
      end
    end
  end
  df = dw/2/pi;
end

% dipolar potential
function uD = f_uD(dat, nu_b)
  wB = 2*pi*nu_b;
  uD = 1/2 * wB^2 * sin(dat.bn).^2;
end


% velocity of transverse spin waves propagating
% along cell radius
function ceff = f_ceff(dat, cpar, cper)
  % find lx,ly,lz
  nz = cos(dat.bn);
  nx = sin(dat.bn).*cos(dat.an);
  ny = sin(dat.bn).*sin(dat.an);

  ct=-1/4; st=sqrt(15)/4;
  lx=(1-ct)*nz.*nx - st*ny;
  ly=(1-ct)*nz.*ny + st*nx;
  lz= ct + (1-ct)*nz.^2;

  ceff=sqrt(cper^2 + (cpar^2-cper^2)*lx.^2); % cm/s
end

% Zero function for finding state number N:
% 1/pi * integral(k) - 1/2 - N
% loc is set to 1 if the state is localized in the texture
function [res,loc] = kint(N, rr, k2)
  % localization of the wave
  ii=find(k2>0);
  if length(ii)<1; loc=0; res=-1; return; end

  % find zero crossing if any
  j = find(k2(2:end).*k2(1:end-1)<0, 1);
  if length(j);
    rc = rr(j) + (rr(j+1)-rr(j)) * k2(j)/(k2(j)-k2(j+1));
    kc = 0;
    rr= [rr(ii); rc];
    kk= [sqrt(k2(ii)); kc];

    % Near rc the integrand has form sqrt(1-r/rc).
    % We calculate it in polar coordinates k/kmax vs r/rc
    % int(f(x)dx) = 1/2 int(R^2 d phi)
    kmax=kk(1);
    ph = atan2(rr/rc, kk/kmax);
    R  = hypot(rr/rc, kk/kmax).^2;
    res = 0.5*kmax*rc*trapz(ph,R)/pi - 0.5 - N;
    % simple formila can be used, but it is less accurate:
    % res = trapz(rr,kk)/pi - 0.5 - N;
    loc = 1;
  else
    res = trapz(rr,sqrt(k2))/pi - 0.5 - N;
    loc = 0;
  end
end

