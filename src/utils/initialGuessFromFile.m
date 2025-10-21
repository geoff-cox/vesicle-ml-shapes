function initSol = initialGuessFromFile(params, H0)
% INITIALGUESSFROMFILE  Load a seed solution from InitialShapes/.
% Usage:
%   initSol = initialGuessFromFile(params, [H0_1 H0_2])
% Looks for:
%   1) params.initialGuess (explicit filename, with or without path)
%   2) the nearest SIM_Node_* file in InitialShapes/ (by |H0_1 - angle| if present)
%
% Expects files saved with a struct like Version(1).Solution.

    % 0) project folders
    here  = fileparts(mfilename('fullpath'));
    root  = fileparts(fileparts(here));           % .../src/utils -> project root
    ishapes = fullfile(root,'InitialShapes');

    % 1) explicit filename provided
    if isfield(params,'initialGuess') && ~isempty(params.initialGuess)
        cand = params.initialGuess;
        if exist(cand,'file')~=2
            cand = fullfile(ishapes, cand);
        end
        assert(exist(cand,'file')==2, 'initialGuessFromFile: file not found: %s', cand);
        tmp = load(cand);
        initSol = tmp.Version(1).Solution;
        return
    end

    % 2) pick a reasonable default from InitialShapes/
    d = dir(fullfile(ishapes, 'SIM_Node_*.mat'));
    assert(~isempty(d), [...
        'initialGuessFromFile: ' ...
        'no SIM_Node_*.mat found in %s ' ishapes]);

    if isfield(params,'A') && isfield(params,'V') && isfield(params,'KG') ...
       && isfield(params,'KA') && isfield(params,'KB')
        base = sprintf('SIM_Node_%d_%d_%d_%d_%d_+00_+00.mat', ...
            round(params.A*100), round(params.V*100), ...
            params.KG, params.KA, params.KB);
        cand = fullfile(ishapes, base);
        if exist(cand,'file')==2
            tmp = load(cand);
            initSol = tmp.Version(1).Solution;
            return;
        end
    end

    % If H0 provided, choose nearest by first component; otherwise take the first
    if nargin >= 2 && ~isempty(H0) && numel(H0)>=1 && isfinite(H0(1))
        angles = nan(numel(d),1);
        for k=1:numel(d)
            % Try to parse angle from filename SIM_Node_XX_...
            tokens = regexp(d(k).name, 'SIM_Node_(\d+)_', 'tokens','once');
            if ~isempty(tokens), angles(k) = str2double(tokens{1}); end
        end
        [~,ix] = min(abs(angles - H0(1)), [], 'omitnan');
        if isempty(ix) || isnan(ix), ix = 1; end
    else
        ix = 1;
    end

    f = fullfile(d(ix).folder, d(ix).name);
    tmp = load(f);
    initSol = tmp.Version(1).Solution;
end
