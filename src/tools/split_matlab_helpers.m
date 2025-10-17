function split_matlab_helpers(inputFile, outputDir, opts)
% SPLIT_MATLAB_HELPERS  Split a multi-function .m file into:
%   - Main file: first (top-level) function block
%   - One .m per *top-level* helper (keeps any nested sub-helpers inside)
%
% Only one level is separated: sub-helpers (nested functions) remain inside
% their parent helper file.
%
% Usage:
%   split_matlab_helpers('sim_driver_quad_tree.m', 'refactored')
%   split_matlab_helpers('sim_driver_quad_tree.m', 'refactored', ...
%       struct('helpersSubdir','helpers','writeCleanMain',true, ...
%              'addpathInMain',true,'overwrite',false, ...
%              'topLevelMustBeCol1',true))
%
% Key idea:
%   - Detect *top-level* function headers and their matching (column-1) ENDs.
%   - Anything between a helper's header and its top-level END stays together,
%     including any nested helper functions.
%
% Notes:
%   - By default, top-level headers must start in column 1 (no leading spaces).
%     If your file-scope helpers are indented, set opts.topLevelMustBeCol1=false.
%   - If a top-level function omits its END, we stop at the next top-level header.
%
% -------------------------------------------------------------------------

    arguments
        inputFile (1,1) string
        outputDir (1,1) string
        opts.helpersSubdir (1,1) string = "helpers"
        opts.writeCleanMain (1,1) logical = true
        opts.addpathInMain (1,1) logical = true
        opts.overwrite (1,1) logical = true
        opts.topLevelMustBeCol1 (1,1) logical = true
    end

    % --- Read file
    assert(isfile(inputFile), 'Input file not found: %s', inputFile);
    raw = fileread(inputFile);
    lines = splitlines(raw);
    n = numel(lines);

    % --- Helpers: quick predicates
    isCol1 = @(s) ~isempty(s) && (s(1) ~= ' ' && s(1) ~= sprintf('\t'));
    isHeaderLine = @(L) ~startsWith(strtrim(L), "%") && ~startsWith(strtrim(L), "%%") ...
                      && startsWith(strtrim(L), "function", 'IgnoreCase', true);
    isTopHeader = @(L) isHeaderLine(L) && ( ~opts.topLevelMustBeCol1 || isCol1Char(L) );

    function tf = isCol1Char(L)
        tf = ~isempty(L) && ~startsWith(L, " ") && ~startsWith(L, sprintf('\t'));
    end

    function tf = isEndAtCol1(L)
        s = strtrim(L);
        tf = isCol1Char(L) && strcmp(s, "end");
    end

    % --- Scan once to build top-level blocks
    % We treat a "top-level header" as a header seen when not currently inside
    % another top-level function; we close a top-level block only on a column-1 "end".
    blocks = struct('start', {}, 'stop', {}, 'name', {});
    openIdx = [];
    for i = 1:n
        L = lines{i};

        if isempty(openIdx)
            % Not inside a top-level block
            if isTopHeader(L)
                name = parseFuncName(L);
                if strlength(name) == 0
                    error('Could not parse function name at line %d.', i);
                end
                openIdx = i;
                blocks(end+1).start = i; %#ok<AGROW>
                blocks(end).name = name;
            end
        else
            % Inside a top-level block
            if isEndAtCol1(L)
                blocks(end).stop = i;      % close on column-1 end
                openIdx = [];
            end
        end
    end

    % If the last open block had no column-1 end, close it at EOF
    if ~isempty(openIdx)
        blocks(end).stop = n;
    end

    if isempty(blocks)
        error('No top-level function definitions found in: %s', inputFile);
    end

    % --- Post-process: ensure stops don't overlap next headers
    % If a block has no explicit col-1 end, trim it to the line before the next block start.
    for k = 1:numel(blocks)-1
        if blocks(k).stop >= blocks(k+1).start
            blocks(k).stop = blocks(k+1).start - 1;
        end
    end

    % --- Summaries
    mainName = blocks(1).name;
    helperNames = string({blocks(2:end).name});

    % --- Prepare output
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    helpersDir = fullfile(outputDir, opts.helpersSubdir);
    if ~exist(helpersDir, 'dir'), mkdir(helpersDir); end

    % --- Write main (cleaned)
    if opts.writeCleanMain
        mainPath = fullfile(outputDir, mainName + ".m");
        if ~opts.overwrite && isfile(mainPath)
            error('Refusing to overwrite existing main: %s (set opts.overwrite=true)', mainPath);
        end
        mainBlock = join(lines(blocks(1).start:blocks(1).stop), newline);
        if opts.addpathInMain && ~isempty(helperNames)
            mainBlock = injectAddpath(string(mainBlock), opts.helpersSubdir);
        end
        writeText(mainPath, addBanner(mainBlock, inputFile, mainName, true));
    end

    % --- Write helpers (each keeps its own nested sub-helpers intact)
    for k = 2:numel(blocks)
        fname = blocks(k).name;
        helperPath = fullfile(helpersDir, fname + ".m");
        if ~opts.overwrite && isfile(helperPath)
            error('Refusing to overwrite existing helper: %s (set opts.overwrite=true)', helperPath);
        end
        helperBlock = join(lines(blocks(k).start:blocks(k).stop), newline);
        writeText(helperPath, addBanner(helperBlock, inputFile, fname, false));
    end

    % --- Done
    fprintf('[OK] Wrote main "%s.m" and %d helper(s) to:\n  %s\n', ...
        mainName, numel(helperNames), outputDir);
    if ~isempty(helperNames)
        fprintf('Helpers folder: %s\n', helpersDir);
        for nm = helperNames, fprintf('  - %s.m\n', nm); end
    end
end

% =======================
% Local utilities
% =======================

function name = parseFuncName(headerLine)
    % Accepts:
    %   function y = foo(x)
    %   function [y1,y2] = foo(x)
    %   function foo(x)
    %   function foo
    L = char(strtrim(headerLine));
    pct = strfind(L, '%'); if ~isempty(pct), L = L(1:pct(1)-1); end
    expr = '^\s*function\s+(?:\[[^\]]*\]\s*=|\w+\s*=)?\s*([A-Za-z]\w*)';
    tok = regexp(L, expr, 'tokens', 'once');
    if ~isempty(tok), name = string(tok{1}); else, name = ""; end
end

function out = addBanner(blockText, sourceFile, funcName, isMain)
    timeStr = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    role = tern(isMain, 'CLEANED MAIN', 'EXTRACTED HELPER');
    banner = sprintf([ ...
        '%% -------------------------------------------------------------------------\n' ...
        '%% %s for "%s"\n' ...
        '%%   - Source: %s\n' ...
        '%%   - Extracted: %s\n' ...
        '%%   - Sub-helpers (nested functions) are retained in this file.\n' ...
        '%% -------------------------------------------------------------------------\n\n'], ...
        role, funcName, sourceFile, timeStr);
    out = banner + string(blockText);
end

function writeText(path, text)
    fid = fopen(path, 'wt'); assert(fid>=0, 'Could not open: %s', path);
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s', text);
end

function v = tern(cond, a, b), if cond, v=a; else, v=b; end, end

function txtOut = injectAddpath(mainBlock, helpersSubdir)
    lines = splitlines(string(mainBlock));
    if isempty(lines), txtOut = mainBlock; return; end
    hdrIdx = -1;
    for i = 1:numel(lines)
        Ls = strtrim(lines(i));
        if Ls ~= "" && startsWith(Ls, "function", 'IgnoreCase', true)
            hdrIdx = i; break;
        end
    end
    if hdrIdx < 0, txtOut = mainBlock; return; end
    insIdx = hdrIdx + 1;
    while insIdx <= numel(lines)
        Ls = strtrim(lines(insIdx));
        if Ls == "" || startsWith(Ls, "%") || startsWith(Ls, "%%")
            insIdx = insIdx + 1;
        else
            break;
        end
    end
    addLine = sprintf("addpath(fullfile(fileparts(mfilename('fullpath')), '%s'));", helpersSubdir);
    newLines = [lines(1:insIdx-1); addLine; lines(insIdx:end)];
    txtOut = strjoin(newLines, newline);
end
