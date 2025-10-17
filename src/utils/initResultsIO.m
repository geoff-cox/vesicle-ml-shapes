% -------------------------------------------------------------------------
% EXTRACTED HELPER for "initResultsIO"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function S = initResultsIO(S)
    run_id = datestr(now,'yyyy-mm-dd_HHMMSS');
    S.paths.run    = fullfile(S.paths.root, ['run_', run_id], filesep);
    S.paths.sol_dir= fullfile(S.paths.run, 'solutions');
    if ~exist(S.paths.run,'dir'), mkdir(S.paths.run); end
    if ~exist(S.paths.sol_dir,'dir'), mkdir(S.paths.sol_dir); end

    % --- define BOTH catalogs explicitly
    S.paths.catalog_mat = fullfile(S.paths.run, 'catalog.mat');
    S.paths.catalog_csv = fullfile(S.paths.run, 'catalog.csv');

    % --- initialize MAT table with a canonical schema
    if ~exist(S.paths.catalog_mat,'file')
        T = table( ...
            strings(0,1), strings(0,1), ...           % timestamp, hash
            zeros(0,1), zeros(0,1), ...               % H0_1, H0_2
            zeros(0,1), zeros(0,1), ...               % A, V
            zeros(0,1), zeros(0,1), zeros(0,1), ...   % KA, KB, KG
            zeros(0,1,'int16'), ...                   % label
            zeros(0,1), zeros(0,1), ...               % E, P
            zeros(0,1), zeros(0,1), ...               % BCmax, DEmax
            zeros(0,1,'uint16'), ...                  % mesh
            'VariableNames', {'timestamp','hash','H0_1','H0_2','A','V','KA','KB','KG', ...
                              'label','E','P','BCmax','DEmax','mesh'});
        save(S.paths.catalog_mat,'T','-v7.3');
    end

    % --- initialize CSV with a header
    if ~exist(S.paths.catalog_csv,'file')
        hdr = {'timestamp','hash','H0_1','H0_2','A','V','KA','KB','KG','label','E','P','BCmax','DEmax','mesh'};
        fid = fopen(S.paths.catalog_csv,'wt');
        fprintf(fid, '%s\n', strjoin(hdr,','));
        fclose(fid);
    end
end