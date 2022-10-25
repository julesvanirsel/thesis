function process(direc,options)
arguments
    direc (1,:) char {mustBeFolder}
    options.start (1,1) double {mustBeNonempty} = -1
    options.cad (1,1) double {mustBeNonempty} = -1
    options.stop (1,1) double {mustBeNonempty} = -1
    options.plot (1,1) logical {mustBeNonempty} = 1
    options.video (1,1) logical {mustBeNonempty} = 1
    options.vtk (1,1) logical {mustBeNonempty} = 0
end

tic

if options.plot
    plot_run(direc,"all",300e3,start=options.start,stop=options.stop,cad=options.cad)
end

if options.video
    images2video(direc,'plots_conductance');
    images2video(direc,'plots_continuity');
    images2video(direc,'plots_contour-auto');
    images2video(direc,'plots_contour-standard');
    images2video(direc,'plots_density');
    images2video(direc,'plots_fccp');
    images2video(direc,'plots_jouleheating');
    images2video(direc,'plots_multipanel-auto');
    images2video(direc,'plots_multipanel-standard');
end

if options.vtk
    gemini3d.write.vtk(direc,'njv')
end

toc

end