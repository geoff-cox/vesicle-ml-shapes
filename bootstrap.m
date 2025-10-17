% bootstrap.m
% -------------------------------------------------------------
% One-stop setup for the vesicle simulation project.
% Run this once per MATLAB session before launching the driver.
% -------------------------------------------------------------

% Reset MATLAB path to a clean state
restoredefaultpath;
rehash toolboxcache;

% Add your source folders
addpath(fullfile(pwd, 'src'));
addpath(fullfile(pwd, 'src/helpers'));
addpath(fullfile(pwd, 'bvp6c-solver'));
addpath(fullfile(pwd, 'InitialShapes'));

% Ensure output folder exists
if ~exist('SimResults', 'dir')
    mkdir('SimResults');
end

addpath(fullfile(pwd, 'SimResults'));

% Optional: show MATLAB version and environment
disp('MATLAB environment initialized for Vesicle Simulation Project.');
ver

% Confirm paths (for debug)
disp('Active paths:');
disp(path);
