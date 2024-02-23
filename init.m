% addpath('\\wsl$\Ubuntu-22.04\home\julesvanirsel\gemini\mat_gemini')
% addpath('\\wsl$\Ubuntu-22.04\home\julesvanirsel\gemini\mat_gemini-scripts')
% setenv('GEMINI_ROOT','\\wsl$\Ubuntu\home\julesvanirsel\gemini')
setenv('CMAKE_PREFIX_PATH','\\dartfs-hpc\rc\lab\L\LynchK\gemini\libgem\')
user = input('User: ');

if ispc
    if user==1
        addpath('C:\files\research\gemini\mat_gemini')
        addpath('C:\files\research\gemini\mat_gemini-scripts')
        addpath('C:\files\research\thesis\')
    elseif user==2
        addpath('D:\Files\research\gemini\mat_gemini')
        addpath('D:\Files\research\gemini\mat_gemini-scripts')
        addpath('D:\Files\research\thesis\')
    else
        fprintf('Unknown user.\n')
    end
else
    addpath('//dartfs-hpc/rc/lab/L/LynchK/gemini/mat_gemini')
    addpath('//dartfs-hpc/rc/lab/L/LynchK/gemini/mat_gemini-scripts')
    addpath('//dartfs-hpc/rc/lab/L/LynchK/Jules/thesis/')
    setenv('EDITOR','vi')
end
clear("user")
setup