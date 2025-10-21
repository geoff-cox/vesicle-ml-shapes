function [task, cache] = processQuadtree(cache, T)
% PROCESSQUADTREE  Scheduler for parameter points using a quadtree.
% 
% Contract (new):
%   [task, cache] = processQuadtree(cache, T)
%     - cache: persistent state; will contain cache.QT.{queue,cells} and cache.config
%     - T: catalog table from catalog_load(simDir) (3-col or legacy folded). Each row has:
%         T.hash, T.timestamp, T.entry (with .params and .meta fields)
%     - task: struct with at least .params (empty [] if no pending work)
%
% This function NEVER solves. It only:
%   1) pops a cell from the queue,
%   2) checks which corners are solved (by consulting T),
%   3) if an unsolved corner exists -> returns that .params as next task and requeues cell,
%   4) if all corners are solved -> evaluates uniformity and either finalizes or subdivides.
%
% Required utils available on path:
%   defaultArg.m, subdivideCell.m, uniformTest.m
%
% cache.config (defaults provided here if missing):
%   .bounds = [H0_1_min H0_1_max; H0_2_min H0_2_max]
%   .maxDepth  (default 6)
%   .maxCells  (default 2000)
%   .eTol      (default 5e-3)
%   .pTol      (default 5e-3)
%   .shapeTau  (default 0.08)
%
% Cell layout (stored in cache.QT.queue):
%   C.corners        : 4x2 [H0_1 H0_2] for (LL, LR, UL, UR)
%   C.depth          : integer
%   C.cornerSolved   : 4x1 logical
%   C.cornerLabel    : 4x1 string
%   C.cornerEnergy   : 4x1 double
%   C.cornerPressure : 4x1 double
%   C.isUniform      : logical
%   C.mixedEdges     : kx2 indices of edges with label transitions

function [task, cache] = processQuadtree(cache, T, MP)
    if nargin < 3, MP = struct(); end
    task = []; cache = ensure_defaults(cache);
    if ~isfield(cache,'QT') || isempty(cache.QT) || isempty(cache.QT.queue)
        if ~isfield(cache,'config') || ~isfield(cache.config,'bounds')
            cache.config.bounds = [-1 1; -1 1];
        end
        cache.QT.queue = {make_root_cell(cache.config.bounds)}; cache.QT.cells = [];
    end

    tol = 0;
    while ~isempty(cache.QT.queue)
        C = cache.QT.queue{1}; cache.QT.queue(1)=[];

        [C, anyUnknown] = refresh_corners_from_catalog(C, T, MP, tol);
        if anyUnknown
            k = find(~C.cornerSolved,1,'first');
            params = struct('H0_1',C.corners(k,1),'H0_2',C.corners(k,2));
            cache.QT.queue = [{C}, cache.QT.queue];
            task = struct('params',params);
            return
        end

        [uniform,mixedEdges] = uniformTest(C, cache.config.eTol, cache.config.pTol, cache.config.shapeTau);
        C.isUniform=uniform; C.mixedEdges=mixedEdges;
        if uniform || C.depth >= cache.config.maxDepth || numel(cache.QT.cells) >= cache.config.maxCells
            cache.QT.cells = [cache.QT.cells; C];
        else
            [C1,C2,C3,C4] = subdivideCell(C);
            cache.QT.queue(end+1:end+4) = {C1,C2,C3,C4};
        end
    end
end