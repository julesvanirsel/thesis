function process(direc,options)
arguments
    direc (1,:) char {mustBeFolder}
    options.start (1,1) double {mustBeNonempty} = -1
    options.cad (1,1) double {mustBeNonempty} = -1
    options.stop (1,1) double {mustBeNonempty} = -1
    options.mlon_ref (1,1) double {mustBeNonempty} = -1
    options.plot (1,1) logical {mustBeNonempty} = 1
    options.video (1,1) logical {mustBeNonempty} = 1
    options.vtk (1,1) logical {mustBeNonempty} = 0
end

tic

if options.plot
    plot_run(direc,"all",300e3...
        ,start = options.start...
        ,stop = options.stop...
        ,cad = options.cad...
        ,mlon_ref = options.mlon_ref...
        )
end

if options.video
    plots_dirs = dir(fullfile(direc,'plots'));
    for i = 3:length(plots_dirs)
        plots_name = plots_dirs(i).name;
        images2video(direc,fullfile('plots',plots_name))
    end
end

if options.vtk
    gemini3d.write.vtk(direc,'njv')
end

toc

end