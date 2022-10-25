% addpath('\\wsl$\Ubuntu-22.04\home\julesvanirsel\gemini\mat_gemini')
% addpath('\\wsl$\Ubuntu-22.04\home\julesvanirsel\gemini\mat_gemini-scripts')
if ispc
    addpath('D:\Files\research\gemini\mat_gemini')
    addpath('D:\Files\research\gemini\mat_gemini-scripts')
else
    addpath('//dartfs-hpc/rc/lab/L/LynchK/gemini/mat_gemini')
    addpath('//dartfs-hpc/rc/lab/L/LynchK/gemini/mat_gemini-scripts')
    setenv('EDITOR','vi')
end
setup