function gm = text1r_qball_gm(dat)
  % calculate integral on (\nabla M_\perp)^2

  gbm=(dat.bm(2:end)-dat.bm(1:end-1)) / (dat.R/(dat.n-1));
  gbm = [gbm; gbm(end)];
  gm=trapz(dat.rr, 2*pi*dat.rr.*cos(dat.bm).^2.*gbm.^2);
end
