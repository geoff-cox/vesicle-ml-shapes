% cleanup.m
% -------------------------------------------------------------
% removes added paths
% -------------------------------------------------------------

function cleanup
    % Remove your source folders
    rmpath(genpath(fullfile(pwd, 'src')));
    rmpath(fullfile(pwd, 'bvp6c-solver'));
    rmpath(fullfile(pwd, 'InitialShapes'));
end