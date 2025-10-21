function catalog_save(simDir, T)
% CATALOG_SAVE  Atomic write of the unified catalog (MAT-only).
    f   = fullfile(simDir, 'catalog.mat');
    tmp = [f '.tmp'];
    save(tmp, 'T', '-v7');   % use -v7.3 if the table will exceed 2 GB
    movefile(tmp, f, 'f');
end
