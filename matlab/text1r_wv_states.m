% Calculate all magnon states in the texture using wave approach
% original of this function (with r-z 2d version, tests, integrals)
% can be fount in he3misc/magnons_bphase
%
% Equation to solve:
%   %    n = 1/pi int k(r) dr - 1/2
%  [c_perp^2 - (c_par^2-c_perp^2) l_r] \nabla_r^2  = w(w-wL) - 1/2 wB^2 \sin^2\beta_n
%  c_eff^2/w \nabla_r^2 + 1/2 wB^2/w \sin^2\beta_n  = (w-wL)
%
% dat - Texture parameter structure. Do not forget to calculate texture before.
% cper, cpar - Spin-wave velocities (usually corresponds to gradient energy parameters lg1,lg2 in dat).
% f0 - Larmor field in cell center [Hz] / frequency (usually corresponds to dat.H).
% f2r  - d2(fL)/dr^2 across the cell (negative for usual trap) [Hz/cm^2].
% nu_b - Leggett frequency [Hz].
% mode - -1: simple mode: c=cper, beta_n = beta_n'*r
%         0: constant freq,
% NN   - interpolate texture to different grid size
% Nmax - Max number of states to return (default 1000).
% Nex  - Max number of non-local states to return (default 0).


% simplify wave_calc.m for 1D radial calculation
function [df, rr, psi] = text1r_wv_states(dat, cper, cpar, f0, f2r, fB, mode, NN, Nmax, Nex)

  rr = linspace(dat.rr(1),dat.rr(end), NN);
  an = interp1(dat.rr, dat.an, rr);
  bn = interp1(dat.rr, dat.bn, rr);
  dr=(rr(end)-rr(1))/(length(rr)-1);

  wB = 2*pi*fB;
  w0 = 2*pi*f0;
  w2r = 2*pi*f2r;

  % texture-dependent values:
  % velocity of transverse spin waves propagating
  % along cell radius, dipolar potential
  if mode == -1
    ceff=cper*ones(size(rr));
    uD = 0.5*wB^2/w0*(dat.db0*rr).^2 + w2r*rr.^2;
  else
    ceff=f_ceff(an,bn, cpar, cper);
    uD = 0.5*wB^2/w0*sin(bn).^2 + w2r*rr.^2;
  end
  %% for f=0 BC ar R: last point has F=0 and is not in the calculation
  %ceff = ceff(1:end-1);
  %uD   = uD(1:end-1);
  % N=length(rr)-1;

  % we are building matrix for rr(1:N)
  N=length(rr);
  Drr = sparse(N,N);
  uDr = sparse(N,N); % 1/r d/dr
  U   = sparse(N,N);
  L   = sparse(N,N);

  % Differential operators on the grid.
  % r>0:
  %   f = ar^2 + br + c
  %
  %   f(r)    = f0  = ar^2 + b r + c
  %   f(r-dr) = fm  = ar^2 - 2 a r dr + a dr^2 + br - b dr + c
  %   f(r+dr) = fp  = ar^2 + 2 a r dr + a dr^2 + br + b dr + c
  %
  %   fp-fm = 2 (2ar + b) dr
  %   fm+fp-2f0 = 2a dr^2
  %
  %   uDr: [1/r df/dr] = 2a + b/r = (fp-fm)/(2 r dr)
  %   Drr: [d2f/dr2] = 2a = (fm+fp-2f0)/dr^2
  %
  % r=0:
  %   fp=fm
  %
  %   f = f0 + ar^2
  %   fp = f0 + a dr^2
  %   [1/r df/dr] = [d2f/dr2] = 2a = 2(fp-f0)/dr^2
  % r=R:
  %   fp=fm
  %
  %   f = f0 + ar^2
  %   fm = f0 + a dr^2
  %   [1/r df/dr] = [d2f/dr2] = 2a = 2(fm-f0)/dr^2

  for ir=1:N % don't use Nr below!
      x0  = ir;
      if ir>1 xrm = ir-1; end
      if ir<N xrp = ir+1; end

      % Drr
      Drr(ir,ir) = -2/dr^2;
      if ir>1 Drr(ir,ir-1) = 1/dr^2; end
      if ir<N Drr(ir,ir+1) = 1/dr^2; end
      if ir==1
        Drr(ir,ir+1) = 2/dr^2; % Neumann boundary condition
      end
      if ir==N
        Drr(ir,ir-1) = 2/dr^2; % Neumann boundary condition
      end

      % uDr
      if ir>1
        uDr(ir,ir-1) = -1/(2*dr*rr(ir));
        if ir<N(1); uDr(ir,ir+1) =  1/(2*dr*rr(ir)); end
      elseif ir==1
        uDr(ir,ir) = -2/(dr^2);
        uDr(ir,ir+1) = 2/(dr^2);
      else
        uDr(ir,ir) = -2/(dr^2);
        uDr(ir,ir-1) = 2/(dr^2);
      end

      % U
      U(ir,ir)  = uD(ir);
  end

  % non-modified equation
  A = - diag(ceff.^2)*(Drr + uDr)/(2*pi*f0) + U;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % solve eigenvalue problem
  [v,l] = eig(A);
  df=diag(l)/(2*pi);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % sort states, add r=R point where F=0
  [df,ii] = sort(df);
  % for f=0 BC at R
  %  psi = [ real(v(:,ii)); zeros(size(v(1,:))) ];
  psi = real(v(:,ii));

  % wave is localized if derivative at r=R is large
  for i = 1:length(df)
    if psi(1,i) < 0; psi(:,i) = - psi(:,i); end
    loc(i) = abs( psi(end,i)/psi(1,i) );
    % f=0 BC: loc(i) = abs( dr*psi(end-1,i)/psi(1,i) );
  end
  loc = loc > 0.01;

  cut = min([find(loc, 1) + Nex, Nmax, length(loc)]);
  df  = df(1:cut);
  psi = psi(:,1:cut);

end

% velocity of transverse spin waves propagating
% along cell radius
function ceff = f_ceff(an, bn, cpar, cper)
  % find lx,ly,lz
  nz = cos(bn);
  nx = sin(bn).*cos(an);
  ny = sin(bn).*sin(an);

  ct=-1/4; st=sqrt(15)/4;
  lx=(1-ct)*nz.*nx - st*ny;
  ly=(1-ct)*nz.*ny + st*nx;
  lz= ct + (1-ct)*nz.^2;

  ceff=sqrt(cper^2 + (cpar^2-cper^2)*lx.^2); % cm/s
end
