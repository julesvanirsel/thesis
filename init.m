% addpath('\\wsl$\Ubuntu-22.04\home\julesvanirsel\gemini\mat_gemini')
% addpath('\\wsl$\Ubuntu-22.04\home\julesvanirsel\gemini\mat_gemini-scripts')
if ispc
    addpath('D:\Files\research\gemini\mat_gemini')
    addpath('D:\Files\research\gemini\mat_gemini-scripts')
    addpath('D:\Files\research\thesis\')
    direc = '\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\aurora_null_02';
else
    addpath('//dartfs-hpc/rc/lab/L/LynchK/gemini/mat_gemini')
    addpath('//dartfs-hpc/rc/lab/L/LynchK/gemini/mat_gemini-scripts')
    addpath('//dartfs-hpc/rc/lab/L/LynchK/Jules/thesis/')
    setenv('EDITOR','vi')
    direc = '/dartfs-hpc/rc/lab/L/LynchK/public_html/Gemini3D/aurora_null_02';
end
setup