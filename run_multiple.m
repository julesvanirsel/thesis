function run_multiple(varargin,options)
arguments (Repeating)
    varargin (1,:) char {mustBeFolder}
end
arguments
    options.process (1,1) logical {mustBeNonempty} = false
    options.np (1,1) int16 {mustBePositive} = 36
end

for i = 1:nargin
    gemini_bin = fullfile(getenv('GEMINI_ROOT'),'build','gemini.bin');
    system(['mpiexec -np ',num2str(options.np),' ',gemini_bin,' ',varargin{i}],'-echo')
end
if options.process
    for i = 1:nargin
        process(varargin{i})
    end
end
end