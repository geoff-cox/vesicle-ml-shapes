function catalog_save(simDir, T)
    save(fullfile(simDir,'catalog.mat'),'T');
end