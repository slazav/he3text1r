% Calculate all magnon states in the texture using wave approach
% original of this function (with r-z 2d version, tests, integrals)
% can be fount in he3misc/magnons_bphase
%
% Only constant frequency calculation.
%
% Equation to solve:
%   %    n = 1/pi int k(r) dr - 1/2
%  -[c_perp^2 - (c_par^2-c_perp^2) l_r] \nabla_r^2  = w(w-wL) - 1/2 wB^2 \sin^2\beta_n
%  -c_eff^2/w \nabla_r^2 + 1/2 wB^2/w \sin^2\beta_n  = (w-wL)
%
% dat - Texture parameter structure. Do not forget to calculate texture before.
% cper, cpar - Spin-wave velocities (usually corresponds to gradient energy parameters lg1,lg2 in dat).
% f0 - frequency (usually corresponds to dat.H).
% fprof - non-uniform field profile (in freq units, 0 in the center, negative for a usual traps).
% nu_b - Leggett frequency [Hz].
% mode - -2: simple mode: c=cper, beta_n = beta_n'*r, f'(r)=0
%        -1: simple mode: c=cper, beta_n = beta_n'*r, f(r)=0
%         1: f(r)=0
%         2: f'(r)=0
% NN   - interpolate texture to different grid size
% Nmax - Max number of states to return (default 1000).
% Nex  - Max number of non-local states to return (default 0).


% simplify wave_calc.m for 1D radial calculation
function [df, rr, psi] = text1r_wv_states(dat, cper, cpar, f0, fprof, fB, mode, NN, Nmax, Nex)

  wB = 2*pi*fB;
  w0 = 2*pi*f0;
  wprof = 2*pi*fprof;

  rr = linspace(dat.rr(1),dat.rr(end), NN);
  an = interp1(dat.rr, dat.an, rr);
  bn = interp1(dat.rr, dat.bn, rr);
  if length(wprof)>1
    wprof = interp1(dat.rr, wprof, rr);
  end
  dr=(rr(end)-rr(1))/(length(rr)-1);

  % texture-dependent values:
  % velocity of transverse spin waves propagating
  % along cell radius, dipolar potential

  simple = mode<0;
  ZBC = abs(mode)==1;

  if simple
    ceff=cper*ones(size(rr));
    uD = 0.5*wB^2/w0*(dat.db0*rr).^2 + wprof;
  else
    ceff=f_ceff(an,bn, cpar, cper);
    uD = 0.5*wB^2/w0*sin(bn).^2 + wprof;
  end
  % f=0 BC on 
  % last point has f=0 and is not in the calculation
  if ZBC
    ceff = ceff(1:end-1);
    uD   = uD(1:end-1);
    N=length(rr)-1;
  % f'=0 BC on 
  else
    N=length(rr);
  end

  % we are building matrix for rr(1:N)
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
  %   f = f0 + ar^2
  %   fp = fm = f0 + a dr^2
  %   [1/r df/dr] = [d2f/dr2] = 2a = 2(fp-f0)/dr^2
  % r=R, f'=0:
  %   fp = fm = f0 + a dr^2
  %   [1/r df/dr] = [d2f/dr2] = 2a = 2(fm-f0)/dr^2

  for ir=1:N % don't use NN below!

      % Drr
      Drr(ir,ir) = -2/dr^2;
      if ir>1 Drr(ir,ir-1) = 1/dr^2; end
      if ir<N Drr(ir,ir+1) = 1/dr^2; end
      if ir==1
        Drr(ir,ir+1) = 2/dr^2; % Neumann boundary condition
      end
      if ir==N && ~ZBC
        Drr(ir,ir-1) = 2/dr^2; % Neumann boundary condition
      end

      % uDr
      if ir>1 && ir<N
        uDr(ir,ir-1) = -1/(2*dr*rr(ir));
        uDr(ir,ir+1) =  1/(2*dr*rr(ir));
      end
      if ir==1
        uDr(ir,ir) = -2/(dr^2);
        uDr(ir,ir+1) = 2/(dr^2);
      end
      if ir==N && ~ZBC
        uDr(ir,ir) = -2/(dr^2);
        uDr(ir,ir-1) = 2/(dr^2);
      end

      % U
      U(ir,ir)  = uD(ir);
  end
  % non-modified equation
  A = - diag(ceff.^2)*(Drr + uDr)/w0 + U;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % solve eigenvalue problem
  [v,l] = eig(A);
  df=diag(l)/(2*pi);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % sort states, add r=R point where F=0
  [df,ii] = sort(df);
  if ZBC
    psi = [ real(v(:,ii)); zeros(size(v(1,:))) ];
  else
    psi = real(v(:,ii));
  end

  % wave is localized if derivative at r=R is large
  for i = 1:length(df)
    if psi(1,i) < 0; psi(:,i) = - psi(:,i); end
    if ZBC
      loc(i) = abs(psi(end-1,i)) / max(abs(diff(psi(:,i))));
    else
      loc(i) = abs(psi(end,i)/max(psi(:,i)) );
    end
  end
  loc = loc > 0.01;

  cut = min([find(loc, 1) + Nex, Nmax, N]);
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
