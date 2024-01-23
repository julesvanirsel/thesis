function make_sbatch_script(direc,opts)
arguments
    direc (1,:) char {mustBeFolder}
    opts.num_cpus (1,1) int32 {mustBePositive} = 64
    opts.num_nodes (1,1) int32 {mustBePositive} = 1
    opts.time_limit (1,1) int32 {mustBePositive} = 600
    opts.do_setup (1,1) logical = true
end

%% init
gem_root = getenv('GEMINI_ROOT');
mat_root = getenv('GEMINI_MAT_ROOT');
scr_root = getenv('GEMINI_SCR_ROOT');
sim_root = getenv('GEMINI_SIM_ROOT');

assert(~isempty(gem_root), ...
    'Add environment variable GEMINI_ROOT directing to gemini.bin')
assert(~isempty(mat_root), ...
    ['Add environment variable GEMINI_MAT_ROOT directing to contents of ' ...
    'https://github.com/gemini3d/mat_gemini'])
assert(~isempty(scr_root), ...
    ['Add environment variable GEMINI_SCR_ROOT directing to contents of ' ...
    'https://github.com/gemini3d/mat_gemini-scripts'])
assert(~isempty(sim_root), ...
    'Add environment variable GEMINI_SIM_ROOT directing to gemini simulations')

% setup for gemini matlab tools
addpath(mat_root)
addpath(scr_root)
setup

%% setup simulation
if any(direc(end)=='/\')
    direc = direc(1:end-1);
end
[~,sim_name] = fileparts(direc);
path = fullfile(sim_root,sim_name);
if opts.do_setup
    gemini3d.model.setup(path,path)
end

%% generate sbatch script
gem_bin = fullfile(gem_root,'build','gemini.bin');
fprintf('binary file: %s\n',gem_bin)
batch_cmd = '#SBATCH';
mpi_cmd = 'mpiexec';
jules_root = '//dartfs-hpc/rc/lab/L/LynchK/Jules/thesis';
matlab_cmd = sprintf('process(''%s'')',path);

cfg = gemini3d.read.config(path);
xg = gemini3d.read.grid(path);
lx2 = xg.lx(2);
lx3 = xg.lx(3);
if all(mod([lx2,lx3],opts.num_cpus)==0)
    num_cpus = opts.num_cpus;
else
    num_cpu_options = tools.common_factors(xg.lx(2),xg.lx(3));
    [~,num_cpu_id] = min(abs(num_cpu_options - opts.num_cpus));
    num_cpus = num_cpu_options(num_cpu_id);
    warning(['Cannot factor %i into both %i and %i. ' ...
        'Using %i workers instead.'],opts.num_cpus,lx2,lx3,num_cpus)
end

for d = cfg.E0_dir%[cfg.prec_dir,cfg.E0_dir]
    simsize_fn = fullfile(path,d,'simsize.h5');
    llon = h5read(simsize_fn,'/llon');
    llat = h5read(simsize_fn,'/llat');
    if all([llon,llat] == [lx2,lx3])
        fprintf('%s grid matches working grid size.\n',simsize_fn)
        h5write(simsize_fn,'/llon',-1)
        h5write(simsize_fn,'/llat',-1)
    end
end

script_fn = fullfile(path,'batch.script');
fprintf('write sbatch script: %s\n',script_fn)
fid = fopen(script_fn,'w');

fprintf(fid,'#!/bin/bash\n\n');

fprintf(fid,'# Job name\n');
fprintf(fid,'%s -J %s\n\n',batch_cmd,sim_name);

fprintf(fid,'# Output file\n');
fprintf(fid,'%s -o %s%s%%j.log\n\n',batch_cmd,path,filesep);

fprintf(fid,'# Error file\n');
fprintf(fid,'%s -e %s%s%%j.err\n\n',batch_cmd,path,filesep);

fprintf(fid,'# E-mail options\n');
fprintf(fid,'%s --mail-type=END,FAIL\n\n',batch_cmd);

fprintf(fid,'# Time limit (minutes)\n');
fprintf(fid,'%s --time=%i\n\n',batch_cmd,opts.time_limit);

fprintf(fid,'# Number of nodes\n');
fprintf(fid,'%s --nodes=%i\n\n',batch_cmd,opts.num_nodes);

fprintf(fid,'# Number of tasks per node\n');
fprintf(fid,'%s --ntasks-per-node=%i\n\n',batch_cmd,num_cpus);

fprintf(fid,'# Commands to run\n');
fprintf(fid,'%s %s %s;\\\n',mpi_cmd,gem_bin,path);
fprintf(fid,['matlab -nodisplay -nodesktop -r' ...
    ' "\\\naddpath(''%s'');\\\ninit;\\\n%s;\\\nexit\\\n"'], ...
    jules_root,matlab_cmd);

fclose(fid);
end