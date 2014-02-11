function [psi e0] = text1r_wave(dat)
  % Calculate ground state wavefunction Psi (no normalization)
  % for the potential sin(dat.bn).^2
  % psi = sqrt(2S/hbar) sin(bm/2)

  format long;

  N=dat.n-1;
  R=dat.R * 1e-2; % cm->m

  ksi=44.2*10^(-6); %dipolar length
  c=-48/65*(ksi)^2; %the constant in front of the laplace operator

  h   = R/(N+1);  % h=step size
  r   = (h:h:((N-1)*h+h))';   %discretization points starting from h'

  ttV=sin(dat.bn).^2;

  %Interpolating the potential within the signifficant radius
  V=spdiags(ttV(2:N+1),0,N,N);

  % Creating the discretization matrix L for the second derivative:
  % diagonal elements are "-2" and on the left and
  % on the right side of the twos there are "1":s
  tD = sparse(1:N,1:N,-2*ones(1,N),N,N);
  tE = sparse(2:N,1:N-1,ones(1,N-1),N,N);
  L = tE+tD+tE'; %'

  % Boundary conditions for L.
  % We demand S'(0)=0 (symmetry reasons)
  L(1,1)=-2/3;
  L(1,2)=2/3;
  %On the boundary S(R)=S(1)=0, no modification needed
  L=L*1/(h^2);     %correct scale

  %Creating the discretization matrix D for the first derivative : zeros on
  %the diagonal, "1" right from the diagonal and "-1" left from it

  ttE = sparse(2:N,1:N-1,ones(1,N-1),N,N);
  D=-ttE+ttE';%'

  %Boundary conditions for D

  %We demand S(R)=S(1)=0, again no modification needed
  %We want S'(0)=0, thus
  D(1,1)=-2/3;
  D(1,2)=2/3;
  D=D*1/(2*h); %correct scale

  %For the radial laplace we need 1/r * d/dr , not d/dr.
  %Creating a diagonal matrix T of "1/r" :s
  T=spdiags(r.^(-1),0,N,N);

  %Building the complete discr. matrix K
  K=(V+(L+T*D).*c); % c: see for the top

  %Organizing s (<=N) first eigenvectors into matrix "Eigenvectors" and corresponding eigenvectors in vector "E" 
  s=round(N/10);


  if 1==0 %method 1: fast but has convergence problems sometimes
    [Eigenvectors,Energies,flag]=eigs(K,s,'SR');%,opts);
    while flag ~= 0
      [Eigenvectors,Energies,flag]=eigs(K,s,'SR');
    end
    en=zeros(1,s);
    for n=1:s
     en(n)=Energies(n,n);
     eigv=eigv(:,n);
    end
  else %method 2: no convergence problems, slow
    [eigv, en]=eig(full(K));
    [en,ii]=sort(diag(en));
    en=en'; %'
    eigv=eigv(:,ii);
  end
  psi=[eigv(:,1); 0];
  e0=en(1);
end
