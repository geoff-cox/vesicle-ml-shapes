function artifact_id = catalog_add_artifact_from_file(db, resultsRoot, solution_id, kind, format, filepath)
%CATALOG_ADD_ARTIFACT_FROM_FILE Register artifact and compute bytes+sha256.
%
%   artifact_id = catalog_add_artifact_from_file(db, resultsRoot, solution_id, kind, format, filepath)
%
% - resultsRoot is the root directory (e.g., .../SimResults_v2)
% - filepath is the absolute path (recommended) or a path under resultsRoot
% - stores relative path in DB for portability
%
% Example Usage:
%    h5abs = fullfile(resultsRoot,'artifacts','solution', solution_id + ".h5");
%    catalog_add_artifact_from_file(db, resultsRoot, solution_id, "solution_h5", "h5", h5abs);
%

    if nargin < 6
        error('catalog_add_artifact_from_file:args','Expected 6 inputs.');
    end

    resultsRoot = string(resultsRoot);
    filepath    = string(filepath);

    % Normalize and validate file exists
    f = char(filepath);
    if exist(f,'file') ~= 2
        error('catalog_add_artifact_from_file:notfound','File not found: %s', f);
    end

    bytes = file_bytes(f);
    sha   = file_sha256(f);

    relpath = make_relative_path(resultsRoot, filepath);

    artifact_id = catalog_add_artifact(db, solution_id, kind, format, relpath, sha, bytes);
end

function rel = make_relative_path(root, p)
%MAKE_RELATIVE_PATH Make p relative to root; error if not under root.

    root = char(string(root));
    p    = char(string(p));

    root = normalize_sep(root);
    p    = normalize_sep(p);

    % Ensure root ends with file separator
    if ~endsWith(root, filesep)
        root = [root filesep];
    end

    if startsWith(p, root)
        rel = string(p(numel(root)+1:end));
        rel = replace(rel, "\", "/"); % store in DB with forward slashes
    else
        % If user passed a relative path already, allow it if it exists under root
        if ~isabsolute(p)
            candidate = fullfile(root, p);
            if exist(candidate,'file')==2
                rel = string(p);
                rel = replace(rel, "\", "/");
                return;
            end
        end
        error('catalog_add_artifact_from_file:outsideRoot', ...
              'Artifact path must be under resultsRoot. resultsRoot=%s, filepath=%s', root, p);
    end
end

function tf = isabsolute(p)
    % Windows absolute: 'C:\' or '\\server\share'
    % Unix absolute: '/'
    if ispc
        tf = (~isempty(regexp(p,'^[A-Za-z]:[\\/]', 'once'))) || startsWith(p,"\\");
    else
        tf = startsWith(p, "/");
    end
end

function s = normalize_sep(s)
    % Use platform sep for matching.
    if ispc
        s = strrep(s,'/','\');
    else
        s = strrep(s,'\','/');
    end
end
