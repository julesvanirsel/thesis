while true
    user = input(['\nAccount IDs:\n' ...
        'Research Computing ... (0)\n' ...
        'Lynch Lab Laptop ..... (1)\n' ...
        'Jules ................ (2)\n' ...
        'Giselle .............. (3)\n' ...
        '\nAccount ID = ']);
    user = int16(user);
    if isinteger(user) & (user >= 0)
        break
    else
        warning('Invalid entry.')
    end
end

if user == 0
    addpath('//dartfs-hpc/rc/lab/L/LynchK/gemini/mat_gemini')
    addpath('//dartfs-hpc/rc/lab/L/LynchK/gemini/mat_gemini-scripts')
    addpath('//dartfs-hpc/rc/lab/L/LynchK/Jules/thesis')
    setenv('EDITOR','vi')
elseif user == 1
    addpath('C:\files\research\gemini\mat_gemini')
    addpath('C:\files\research\gemini\mat_gemini-scripts')
    addpath('C:\files\research\thesis')
elseif user == 2
    addpath('D:\Files\research\gemini\mat_gemini')
    addpath('D:\Files\research\gemini\mat_gemini-scripts')
    addpath('D:\Files\research\thesis')
elseif user == 3
    addpath('/Volumes/rc/lab/L/LynchK/gemini/mat_gemini')
    addpath('/Volumes/rc/lab/L/LynchK/gemini/mat_gemini-scripts')
    addpath('/Volumes/rc/lab/L/LynchK/Jules/thesis')
else
    error('Unknown user.')
end

setup
roots = ["GEMINI_ROOT","GEMINI_MAT_ROOT","GEMINI_SCR_ROOT","GEMINI_SIM_ROOT" ...
    ,"JULES_ROOT"];
fprintf('\nEnvironment variables:\n')
for r = roots
    fprintf('%s %s\n',pad(r+' ',30,'.'),getenv(r))
end

funcs = ["gemini3d.model.setup","gemini3d.read.grid","gemscr.functions.aurora"];
fprintf('\nCommon function locations\n')
for f = funcs
    fprintf('%s %s\n',pad(f+' ',30,'.'),which(f))
end

clear("user","roots","r","funcs","f")
% addpath('\\wsl$\Ubuntu-22.04\home\julesvanirsel\gemini\mat_gemini')
% addpath('\\wsl$\Ubuntu-22.04\home\julesvanirsel\gemini\mat_gemini-scripts')
% setenv('GEMINI_ROOT','\\wsl$\Ubuntu\home\julesvanirsel\gemini\gemini3d')
% setenv('CMAKE_PREFIX_PATH','\\wsl$\Ubuntu\home\julesvanirsel\gemini\libgem')