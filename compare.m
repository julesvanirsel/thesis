function compare(direc0,direc1,options)
arguments
    direc0 (1,:) char {mustBeFolder}
    direc1 (1,:) char {mustBeFolder}
    options.ignore (1,:) string = ["nml","times"]
end

if any(direc0(end)=='/\')
    direc0 = direc0(1:end-1);
end
if any(direc1(end)=='/\')
    direc1 = direc1(1:end-1);
end

[~,run0] = fileparts(direc0);
[~,run1] = fileparts(direc1);
cfg0 = gemini3d.read.config(direc0);
cfg1 = gemini3d.read.config(direc1);
fields0 = string(fieldnames(cfg0));
fields1 = string(fieldnames(cfg1));
n0 = length(fields0);
n1 = length(fields1);

for i = 1:n0
    field = fields0(i);
    if not(any(strcmp(field,options.ignore)))
        if any(strcmp(field,fields1))
            v0 = cfg0.(field);
            v1 = cfg1.(field);
            if v0 ~= v1
                disp([run0,': ',char(field),' = ',num2str(v0)])
                disp([run1,': ',char(field),' = ',num2str(v1),newline])
            end
        else
            warning off backtrace
            warning(['fieldname ',char(field),' does not exist in ',run1,newline])
        end
    end
end

for i = 1:n1
    field = fields1(i);
    if not(any(strcmp(field,fields0)))
        warning(['fieldname ',char(field),' does not exist in ',run0,newline])
    end
end

warning on backtrace
end