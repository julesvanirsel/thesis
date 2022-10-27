% Run a series of Gemini simulations
%
% Example usage:
%   gemini_series('LynchK\public_html\Gemini3D\aurora_null_01','E_amp_h',[2e3,3e3,5e3,10e3])
%
% Arguments:
%   direc   base gemini run directory
%   par     function paramater to change
%   vals    list of values to change par to
%
% Dependencies:
%   matlab R2022a or higher
%
% Contact: jules.van.irsel.gr@dartmouth.edu

function gemini_series(direc,par,vals)
arguments
    direc (1,:) char {mustBeFolder}
    par (1,:) char {mustBeNonzeroLengthText}
    vals (1,:) double {mustBeNonempty}
end

lp = length(par);
fn = gemini3d.find.config(direc);

for i = 1:length(vals)
    val = vals(i);
    new_direc = [direc,'_',erase(par,'_'),'=',num2str(val,'%+6.0e')];
    if ~exist(new_direc,'dir')
        mkdir(new_direc);
    end
    new_fn = fullfile(new_direc,'config.nml');
    new_fid = fopen(new_fn,'w');
    fid = fopen(fn);
    found = false;
    while ~feof(fid)
        line = fgetl(fid);
        if length(line) > lp
            if strcmp(line(1:lp),par)
                found = true;
                line = [line(1:lp+3),num2str(vals(i)),' ! auto generated'];
            end
        end
        fprintf(new_fid,[line,newline]);
    end
    if ~found
        fclose all;
        error(['Parameter ',par,' not found in ',char(fn),'.'])
    end
    gemini3d.model.setup(new_direc)
end
fclose all;
end