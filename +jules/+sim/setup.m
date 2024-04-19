function setup(direc)
arguments
    direc (1,:) char {mustBeFolder}
end
potential_fn = fullfile(direc,'ext','potential.h5');
if ~exist(potential_fn,'file')
    error('Please run jules.sim.replication first.')
end
mat_root = getenv('GEMINI_MAT_ROOT');
assert(~isempty(mat_root), ...
    ['Add environment variable GEMINI_MAT_ROOT directing to contents of ' ...
    'https://github.com/gemini3d/mat_gemini'])
addpath(mat_root)
jules.sim.run_ic(direc)
gemini3d.model.setup(direc,direc)
jules.sim.pbs(direc)
end