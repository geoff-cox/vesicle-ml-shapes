% -------------------------------------------------------------------------
% EXTRACTED HELPER for "writeCatalogRow"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function writeCatalogRow(S, hash, label, E, P, bcmax, demax, mesh)
    % Assemble a single logical row (cell array)
    ts  = datestr(now,'yyyy-mm-dd HH:MM:SS');
    row = {ts, hash, ...
           S.H0(1), S.H0(2), S.params.A, S.params.V, ...
           S.params.KA, S.params.KB, S.params.KG, ...
           label, E, P, bcmax, demax, mesh};

    % ---- 1) Append to MAT table (catalog.mat)
    if isfield(S.paths,'catalog_mat') && exist(S.paths.catalog_mat,'file')==2
        L = load(S.paths.catalog_mat,'T'); 
        T = L.T;
        % Ensure same variable order
        newRow = cell2table(row, 'VariableNames', T.Properties.VariableNames);
        T = [T; newRow]; %#ok<AGROW>
        save(S.paths.catalog_mat,'T','-v7.3');
    end

    % ---- 2) Append to CSV (catalog.csv)
    if isfield(S.paths,'catalog_csv')
        fid = fopen(S.paths.catalog_csv,'at');   % creates if missing
        fmt = '%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%d,%.16g,%.16g,%.3e,%.3e,%d\n';
        fprintf(fid, fmt, row{:});
        fclose(fid);
    end
end