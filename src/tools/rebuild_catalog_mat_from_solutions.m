function rebuild_catalog_mat_from_solutions()
    simDir = fullfile(pwd,'sim-results');
    solDir = fullfile(simDir,'solutions');
    addpath(fullfile(pwd,'src','utils'));
    T = catalog_load(simDir);
    files = dir(fullfile(solDir,'*.mat'));

    for i=1:numel(files)
        f = fullfile(files(i).folder, files(i).name);
        [~, base] = fileparts(f);
        hash = string(base);

        % Load meta/params if present
        params = struct(); meta = struct();
        try S = load(f); catch, warning('Unreadable: %s', f); continue; end
        if isfield(S,'meta'), meta = S.meta; end
        if isfield(S,'params'), params = S.params; end
        if isfield(meta,'hash') && strlength(string(meta.hash))==0, meta.hash = hash; end
        if ~isfield(meta,'version'), meta.version = "legacy"; end

        % Derive quick stats if missing
        if (~isfield(meta,'meshN') || isempty(meta.meshN)) && isfield(S,'result')
            meta.meshN = size(S.result.y, ndims(S.result.y));
        end
        if (~isfield(meta,'residualMax') || isempty(meta.residualMax)) && isfield(S,'result')
            r = []; if isfield(S.result,'bcResidual'), r=[r; S.result.bcResidual(:)]; end
            if isfield(S.result,'deResidual'), r=[r; S.result.deResidual(:)]; end
            if ~isempty(r), meta.residualMax = max(abs(r)); end
        end

        ts = datetime(files(i).datenum,'ConvertFrom','datenum', 'TimeZone','UTC', 'Format','yyyy-MM-dd HH:mm:ss');
        entry = struct('params',params,'meta',meta);
        T = catalog_append(simDir, hash, entry, ts);
    end
    fprintf('[rebuild] Final rows: %d\n', height(catalog_load(simDir)));
end
