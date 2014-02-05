function compile()
  comp('text1r_init', 1);
  comp('text1r_set_vortex_cluster', 2);
  comp('text1r_set_vortex_uniform', 2);
  comp('text1r_set_vortex_twisted', 2);
  comp('text1r_minimize', 3);
end

function comp(fname, ftype)
  fprintf('>> compiling %s (%d args)\n', fname);
  mex( ['-DFUNC=' fname '_' ], ...
       ['-DFTYPE=' num2str(ftype) ], ...
       '-o', fname, '-lhe3', '-ltext1r', ...
       ['-L' pwd ], ['-Wl,--rpath=' pwd ], ...
       'mexfunc.c');
end

