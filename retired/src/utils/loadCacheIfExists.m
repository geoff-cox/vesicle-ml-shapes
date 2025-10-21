% -------------------------------------------------------------------------
% EXTRACTED HELPER for "loadCacheIfExists"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function cache = loadCacheIfExists(S, cache)
    d = dir(fullfile(S.paths.root,'run_*')); % most recent prior run
    if isempty(d), return; end
    [~,i]=max([d.datenum]); prev = fullfile(d(i).folder, d(i).name, 'cache.mat');
    if exist(prev,'file')==2
        L = load(prev,'cache');
        % cat unique keys (avoid duplicates)
        cache = mergeCaches(cache, L.cache);
    end
end