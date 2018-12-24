function mex_text1r()
  comp('text1r_init', 1);
  comp('text1r_set_vortex_cluster', 2);
  comp('text1r_set_vortex_uniform', 2);
  comp('text1r_set_vortex_twisted', 2);
  comp('text1r_minimize',  3);
  comp('text1r_selfcheck', 3);
end

function comp(fname, ftype)
  fprintf('>> compiling %s (type %d)\n', fname, ftype);
  mex( ['-DFUNC=' fname '_' ], ...
       ['-DFTYPE=' num2str(ftype) ], ...
       '-o', ['octave/' fname], '-ltext1r -lhe3',...
       ['-L' pwd ], ['-Wl,--rpath=' pwd ],...
       'mex_text1r.c');
end

