% Run a series of Gemini simulations
%
% Example usage:
%   run_series('LynchK\public_html\Gemini3D\aurora_null_01'...
%               ,["E_amp_h","E_amp_l"]...
%               ,[[2e3,3e3,5e3,10e3];[2e3,3e3,5e3,10e3]])
%
% Arguments:
%   direc           base gemini run directory
%   pars            function paramaters to change
%   vals            list of values to change parameters to
%   setup = true    (option) run gemini3d.model.setup
%
% Dependencies:
%   matlab R2022a or higher
%
% Contact: jules.van.irsel.gr@dartmouth.edu

function run_series(direc,pars,vals,options)
arguments
    direc (1,:) char {mustBeFolder}
    pars (1,:) string {mustBeNonempty}
    vals (:,:) double {mustBeNonempty}
    options.setup (1,1) logical {mustBeNonempty} = true
end

npars = length(pars);
nruns = length(vals);
assert(npars==size(vals,1),'Number of parameters does not match the values array size.')
fn = gemini3d.find.config(direc);

for i = 1:nruns
    fid = fopen(fn);
    if npars==1
        new_direc = [direc,'_',erase(char(pars),'_'),'=',num2str(vals(i),'%+6.0e')];
    else
        new_direc = [direc,'_',char(64+i)];
    end
    if ~exist(new_direc,'dir')
        mkdir(new_direc);
    end
    new_fn = fullfile(new_direc,'config.nml');
    new_fid = fopen(new_fn,'w');
    found = false(size(pars));
    while ~feof(fid)
        line = fgetl(fid);
        for j = 1:npars
            par = char(pars(j));
            val = vals(j,i);
            lp = strlength(par);
            if length(line) > lp
                if strcmp(line(1:lp),par)
                    found(j) = true;
                    line = [line(1:lp+3),num2str(val),' ! auto generated'];
                end
            end
        end
        fprintf(new_fid,[line,newline]);
        %         if options.run
        %             gemini_bin = fullfile(getenv('GEMINI_ROOT'),'build','gemini.bin');
        %             system(['mpiexec -np ',num2str(options.np),' ',gemini_bin,' ',new_direc],'-echo')
        %         end
    end
    if not(all(found))
        fclose all;
        error(['One or more parameters not found in ',char(fn),'.'])
    end
    if options.setup
        gemini3d.model.setup(new_direc,new_direc)
    end
end
fclose all;
end