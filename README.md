# thesis
Ph. D. thesis work

How to run a Gemini simulation on Andes:
1)  "mkdir LynchK/public_html/Gemini3D/<SIMULATION NAME>"
    -   Try and stick to a sensible simulation name such as <function>_<type>_<descriptor>_<version>, e.g. aurora_sharc_wide_02
2)  "cp LynchK/public_html/Gemini3D/<PREVIOUS SIMULATION NAME>/config.nml LynchK/public_html/Gemini3D/<SIMULATION NAME>"
    -   If you're building on an old simulation, copy that simulation's config.nml
3)  Edit config.nml in <SIMULATION NAME> to your liking.
    -   See below for details on config.nml


Configuration file details:

&base
ymd = 2015,2,1              year,month,day
UTsec0 = 35850              initial second of day
tdur = 300                  run duration in seconds
dtout = 10                  output cadance in seconds
activ = 137.6,143.7,20      F10.7 3-month avg, F10.7 daily, Ap daily
tcfl = 0.9                  Courant number
Teinf = 1500                exospheric electron temperature
/

&setup
glat = 65.8                             geographic center latitude
glon = 207.7                            geographic center longitude
xdist = 3400e3                          rough east-west span in meters
ydist = 1200e3                          rough north-south span in meters
alt_min = 80e3                          minimum altitude in meters
alt_max = 950e3                         maximum altitude in meters
alt_scale = 10e3, 8e3, 500e3, 150e3     d1 + d2*tanh((alt - d3)/d4)
x2parms = 400e3, 18.8e3, 50e3, 100e3    d2 + d3*(1 + tanh((x - xdist/2 + d1)/d4))/2 mirrored
x3parms = 400e3, 1.625e3, 18.5e3, 50e3  d2 + d3*(1 + tanh((y - ydist/2 + d1)/d4))/2 mirrored
lxp = 1
lyp = 1
Bincl = 90
eq_dir = '/dartfs-hpc/rc/lab/L/LynchK/Jules/initial_conditions/null_02'
setup_functions = 'gemscr.functions.aurora'
/

&flags
potsolve = 1
flagperiodic = 0
flagoutput = 1
/

&files
file_format = 'h5'
indat_size = 'inputs/simsize.h5'
indat_grid = 'inputs/simgrid.h5'
indat_file = 'inputs/initial_conditions.h5'
/

&precip
dtprec = 1
prec_dir = 'inputs/particles'
/

&efield
dtE0 = 10
E0_dir = 'inputs/fields'
/

&aurora_parameters
driftE = 0
driftN = 0
ctr_spn = 1e5
ctr_slp = 0
ctr_pos = 0
bar_frc = 0
bar_pos = 4.5e5
bar_vel = 0
bar_gsl = 1e5
dim_frc = 0
dim_tim = 100
dim_del = 100
Q_amp_h = 3
Q_amp_l = 3
Q_wth_h = 2e4
Q_wth_l = 2e5
Q_gsl_h = 0.4
Q_gsl_l = 0.1
Q_off_h = -0.6
Q_off_l = 5e4
Q_floor = 0.3
E_amp_h = 3e3
E_amp_l = 3e3
E_wth_h = 2e4
E_wth_l = 2.2e5
E_gsl_h = 0.4
E_gsl_l = 0.1
E_floor = 1.5e3
K_amp = 0.3
J_wth = 3e5
J_gsl = 0.1
F_amp = 2e3
F_wth = 5e4
F_gsl = 0.3
flagdirich = 1
ap_ExBg = 0
ap_EyBg = 0
ap_np = 16
ap_cad2 = 4
ap_cad3 = 1
/
