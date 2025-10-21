% -------------------------------------------------------------------------
% EXTRACTED HELPER for "initFoldersAndLogging"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:11
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function sim = initFoldersAndLogging(sim)
    % INITFOLDERSANDLOGGING  Prepare output directories and logfile.

    switch lower(sim.mode)
        case 'new'
            base = fullfile('SimResults', [datestr(now,'dd-mmm-yyyy'), filesep]);
            sim.paths.root    = base;
            sim.paths.mats    = fullfile(base, 'MATS', filesep);
            sim.paths.matname = sim.paths.mats;   % original var name parity
            
            if ~exist(sim.paths.mats, 'dir'), mkdir(sim.paths.mats); end
            copyfile(fullfile('InitialShapes', sim.initialGuess), ...
                     fullfile(sim.paths.mats, sim.initialGuess));
            logName = fullfile(base, 'OPfile.txt');

            sim.paths.results = fullfile(base, 'results');      % hashed solutions
            sim.paths.catalog = fullfile(base, 'catalog.csv');  % index
            if ~exist(sim.paths.results,'dir'), mkdir(sim.paths.results); end

        case 'continue'
            baseOld = fullfile('SimResults', [sim.contDate, filesep]);
            sim.paths.root = baseOld;
            sim.paths.mats = fullfile(baseOld, 'Expand_MATS', filesep);
            if ~exist(sim.paths.mats,'dir'), mkdir(sim.paths.mats); end
            copyfile(fullfile(baseOld,'MATS'), sim.paths.mats);
            logName = fullfile(baseOld, 'Expand_OPfile.txt');

        case 'rescan'
            baseOld = fullfile('SimResults', [sim.rescanDate, filesep]);
            sim.paths.root = baseOld;
            % Find incrementing ReScan_MATS_N directory
            d = dir(baseOld);
            count = 1;
            for k = 1:numel(d)
                if numel(d(k).name) > 10 && strncmp(d(k).name,'ReScan_MATS',11)
                    count = count + 1;
                end
            end
            sim.paths.mats = fullfile(baseOld, sprintf('ReScan_MATS_%d',count), filesep);
            if ~exist(sim.paths.mats, 'dir'), mkdir(sim.paths.mats); end
            if count > 1
                copyfile(fullfile(baseOld, sprintf('ReScan_MATS_%d',count-1), filesep), sim.paths.mats);
            else
                copyfile(fullfile(baseOld,'MATS', filesep), sim.paths.mats);
            end
            logName = fullfile(baseOld, 'ReScan_OPfile.txt');

        otherwise
            error('Unknown mode: %s', sim.mode);
    end

    % Open log
    sim.log.path = logName;
    sim.log.fid  = fopen(logName, 'wt');
    logmsg(sim, '*** Node Simulation ***');
    logmsg(sim, 'Phase-plane: [%g,%g] x [%g,%g] | kappa_A=%g', ...
        sim.grid.A_lim, sim.grid.B_lim, sim.params.KA);
end