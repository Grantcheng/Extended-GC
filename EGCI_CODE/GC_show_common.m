
od = 20;

gc0 = pos_nGrangerT2(X0, od);
gc  = pos_nGrangerT2(X , od);

fprintf('f X -> Y = %.5f (origin=%.5f)\n', gc(2,1), gc0(2,1));
fprintf('f Y -> X = %.5f (origin=%.5f)\n', gc(1,2), gc0(1,2));
fprintf('  (zero GC level: %.5f)\n', od/length(X));

