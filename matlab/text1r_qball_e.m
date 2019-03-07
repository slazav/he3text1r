% old script! use text1r_wv_states/text1r_qc_states instead of text1r_wave

function [dat e0] = text1r_qball_e(dat, E)
  % calculate self-consistent texture and precessing magnetization
  % distibution for a given integral of (1-\cos\beta)
  % (which is proportional to Zeeman energy)
  %
  % \Psi = sqrt(2S/hbar) sin(bm/2)
  % Zeeman energy is \chi H^2 (1-\cos\beta_M) = hbar omega_L |\Psi|^2

  hbar=1.054571726e-27;
  gyro=20378.0e0;

  ee=0;
  for i=1:100
    % calculate texture
    dat = text1r_minimize(dat);

    % calculate wave function
    [psi e0] = text1r_wave(dat);

    % renormalize wave-function and set beta_m
    int = trapz(dat.rr, 2*pi*dat.rr.*psi.^2);
    psi = psi*sqrt(E/int); % now integral(psi^2) = E

    % here psi=sqrt(2) sin(bm/2),
    %      |psi|^2 = 1-cos(bm)
    dat.bm = 2*asin(psi/sqrt(2));

    if abs(ee-e0) < 1e-10; break; end
    ee=e0;
  end

end
