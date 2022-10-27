% addpath('\\wsl$\Ubuntu-22.04\home\julesvanirsel\gemini\mat_gemini')
% addpath('\\wsl$\Ubuntu-22.04\home\julesvanirsel\gemini\mat_gemini-scripts')
if ispc
    addpath('D:\Files\research\gemini\mat_gemini')
    addpath('D:\Files\research\gemini\mat_gemini-scripts')
    direc = '\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\aurora_null_01';
else
    addpath('//dartfs-hpc/rc/lab/L/LynchK/gemini/mat_gemini')
    addpath('//dartfs-hpc/rc/lab/L/LynchK/gemini/mat_gemini-scripts')
    setenv('EDITOR','vi')
    direc = '/dartfs-hpc/rc/lab/L/LynchK/public_html/Gemini3D/aurora_null_01';
end
setup