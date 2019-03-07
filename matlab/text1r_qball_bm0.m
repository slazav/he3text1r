% old script! use text1r_wv_states/text1r_qc_states instead of text1r_wave

function dat = text1r_qball_bm0(dat, bm0)
  % Calculate self-consistent texture and precessing magnetization
  % distibution for a given beta_m(r=0) value (in radians).
  % psi = sqrt(2S/hbar) sin(bm/2)

  ee=0;
  for i=1:100
    % calculate texture
    dat = text1r_minimize(dat);

    % calculate wave function
    [psi e0] = text1r_wave(dat);

    % renormalize wave-function and set beta_m
    psi = psi/psi(1) * sin(bm0/2);
    dat.bm = 2*asin(psi);

    if abs(ee-e0) < 1e-10; break; end
    ee=e0;
  end

end
