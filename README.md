# thesis
Ph. D. thesis work

## How to run a Gemini simulation on Andes:
1)  "mkdir LynchK/public_html/Gemini3D/<SIMULATION NAME>"
    -   Try and stick to a sensible simulation name such as <function>_<type>_<descriptor>_<version>, e.g. aurora_sharc_wide_02
2)  "cp LynchK/public_html/Gemini3D/<PREVIOUS SIMULATION NAME>/config.nml LynchK/public_html/Gemini3D/<SIMULATION NAME>"
    -   If you're building on an old simulation, copy that simulation's config.nml
3)  Edit config.nml in <SIMULATION NAME> to your liking.
    -   See below for details on config.nml entries
4)  If not already done so, open a MATLAB screen with "screen -S mat"
5)  In MATLAB screen, run "matlab -nodisplay"
6)  In MATLAB, navigate to <SIMULATION NAME> and run "gemini3d.model.setup('.','.')"
7)  If not already done so, open a Gemini screen with "screen -S gem"
8)  In Gemini screen, navigate to LynchK directory
9)  In Gemini screen, run "mpiexec -np 36 gemini/gemini3d/build/gemini.bin public_html/Gemini3D/<SIMULATION NAME>"
    -   Make sure that -np even divides nx2 x nx3, e.g. 144 x 216 / 36 = 864 cells per processor
10) In MATLAB screen, navigate to LynchK/Jules/thesis
11) In MATLAB screen, run "process('../../public_html/Gemini3D/<SIMULATION NAME>',<options>)". Options include
    -   plot=1 for plotting
    -   video=1 for making videos
    -   vtk=0 for making paraview files

## Configuration file details:

### base
- ymd = 2015,2,1 -------- year,month,day
- UTsec0 = 35850 -------- initial second of day
- tdur = 300 -------- run duration in seconds
- dtout = 10 -------- output cadance in seconds
- activ = 137.6,143.7,20 -------- F10.7 3-month avg, F10.7 daily, Ap daily
- tcfl = 0.9 -------- Courant number
- Teinf = 1500 -------- exospheric electron temperature

### setup
- glat = 65.8 -------- geographic center latitude
- glon = 207.7 -------- geographic center longitude
- xdist = 3400e3 -------- rough east-west span in meters
- ydist = 1200e3 -------- rough north-south span in meters
- alt_min = 80e3 -------- minimum altitude in meters
- alt_max = 950e3 -------- maximum altitude in meters
- alt_scale = 10e3, 8e3, 500e3, 150e3 -------- d1 + d2*tanh((alt - d3)/d4)
- x2parms = 400e3, 18.8e3, 50e3, 100e3 -------- d2 + d3*(1 + tanh((x - xdist/2 + d1)/d4))/2 mirrored
- x3parms = 400e3, 1.625e3, 18.5e3, 50e3 -------- d2 + d3*(1 + tanh((y - ydist/2 + d1)/d4))/2 mirrored
- lxp = 1 -------- number of east-west cells if uniform
- lyp = 1 -------- number of north-south cells if uniform
- Bincl = 90 -------- geomagnetic inclination in degrees
- eq_dir = '/dartfs-hpc/rc/lab/L/LynchK/Jules/initial_conditions/null_02' -------- directory to equalibrium run
- setup_functions = 'gemscr.functions.aurora' -------- callable matlab input function

### flags
- potsolve = 1 -------- 0 - no; 1 - electrostatic; 2 - inductive
- flagperiodic = 0 -------- 0 - not periodic in x3; 1 - periodic in x3
- flagoutput = 1 -------- what info in output:  1 - all; 2 - avg plasma parameters; 3 - ne only

### files
- file_format = 'h5' -------- data filename extensions
- indat_size = 'inputs/simsize.h5' -------- location of grid size file
- indat_grid = 'inputs/simgrid.h5' -------- location of grid file
- indat_file = 'inputs/initial_conditions.h5' -------- location of intial conditions file

### precip
- dtprec = 1 -------- precip input cadance
- prec_dir = 'inputs/particles' -------- location of precip files

### efield
- dtE0 = 10 -------- current/potential input cadance
- E0_dir = 'inputs/fields' -------- location of current/potential files

### aurora_parameters
- driftE = 0 -------- TBD
- driftN = 0
- ctr_spn = 1e5
- ctr_slp = 0
- ctr_pos = 0
- bar_frc = 0
- bar_pos = 4.5e5
- bar_vel = 0
- bar_gsl = 1e5
- dim_frc = 0
- dim_tim = 100
- dim_del = 100
- Q_amp_h = 3
- Q_amp_l = 3
- Q_wth_h = 2e4
- Q_wth_l = 2e5
- Q_gsl_h = 0.4
- Q_gsl_l = 0.1
- Q_off_h = -0.6
- Q_off_l = 5e4
- Q_floor = 0.3
- E_amp_h = 3e3
- E_amp_l = 3e3
- E_wth_h = 2e4
- E_wth_l = 2.2e5
- E_gsl_h = 0.4
- E_gsl_l = 0.1
- E_floor = 1.5e3
- K_amp = 0.3
- J_wth = 3e5
- J_gsl = 0.1
- F_amp = 2e3
- F_wth = 5e4
- F_gsl = 0.3
- flagdirich = 1
- ap_ExBg = 0
- ap_EyBg = 0
- ap_np = 16
- ap_cad2 = 4
- ap_cad3 = 1

## How to make an ssh key on andes
1)  ssh-keygen -t ed25519 -C "your_email@example.com"
2)  eval "$(ssh-agent -s)"
3)  ssh-add ~/.ssh/id_ed25519
4)  clip < ~/.ssh/id_ed25519.pub
5)  in github.com/julesvanirsel go to
    -   settings
    -   SSH and GPG keys
    -   New SSH key
    enter a Title and paste the key
6)  git remote set-url origin git@github.com:julesvanirsel/thesis.git
7)  git config --global user.name "your_username"
8)  git config --global user.email "your_email_address@example.com"
