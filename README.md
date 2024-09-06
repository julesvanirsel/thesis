# thesis
Ph. D. thesis work

## How to install inversion infrastructure
Create a new environment
```sh
conda create --name "asi" python=3.10
```

Install 
```sh
conda activate asi
pip install scipy apexpy glob2 scikit-image matplotlib h5py
git clone https://github.com/almule12/asispectralinversion.git
pip install asispectralinversion/src/
```

## How to install GEMINI on WSL
### Install WSL distribution
1)  From Windows Powershell:
```sh
wsl --install Ubuntu
```
2)  Generate a username, \<firstnamelastname\>, and password.

### Install standard packages
3)  In WSL:
```sh
sudo apt update && sudo apt upgrade
sudo apt install build-essential # for c and c++ compiler: gcc and g++
sudo apt install gcc gfortran # for fortran compiler: gfortran
sudo apt install cmake # for cmake
sudo apt install libopenmpi-dev openmpi-bin # for mpi commands: e.g. mpiexec
sudo apt install libhdf5-dev # for h5c++, h5cc, and h5fc
sudo apt-get install liblapack-dev # for lapack
which gcc g++ gfortran cmake mpiexec h5c++ h5cc h5fc # make sure all of these commands exists
```

### Install external packages
4)  In WSL's /home/username:
```sh
mkdir gemini
cd gemini
git clone https://github.com/gemini3d/external.git
cmake -B external/build -S external -DCMAKE_INSTALL_PREFIX=~/gemini/libgem
cmake --build external/build
```
5)  At the end of ~/.bashrc add
```sh
export CMAKE_PREFIX_PATH=~/gemini/libgem
```
6)  Restart WSL.
```sh
echo $CMAKE_PREFIX_PATH
```

### Install GEMINI and run self tests
7)  In WSL's /home/username/gemini:
```sh
git clone https://github.com/gemini3d/gemini3d.git
cd gemini3d
cmake -B build -DCMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH
cmake --build build --parallel
ctest --test-dir build
```

### Install MATLAB packages
8)  In WSL's /home/username/gemini:
```sh
git clone --recurse-submodules https://github.com/gemini3d/mat_gemini
git clone --recurse-submodules https://github.com/gemini3d/mat_gemini-scripts
```

## Mount LynchK to WSL (Ubuntu)
```sh
sudo apt install cifs-utils
sudo mount -t cifs -o username=f123451,domain=KIEWIT.DARTMOUTH.EDU -o vers=3.0,file_mode=0660,dir_mode=0770,uid=julesvanirsel //dartfs-hpc.dartmouth.edu/rc/lab/L/LynchK /mnt/lynchk
```

## How to run a Gemini simulation on Andes:
1)  "mkdir LynchK/public_html/Gemini3D/\<SIMULATION NAME\>"
    -   Try and stick to a sensible simulation name such as \<function\>\_\<type\>\_\<descriptor\>\_\<version\>, e.g. aurora_sharc_wide_02
2)  "cp LynchK/public_html/Gemini3D/\<PREVIOUS SIMULATION NAME\>/config.nml LynchK/public_html/Gemini3D/\<SIMULATION NAME\>"
    -   If you're building on an old simulation, copy that simulation's config.nml
3)  Edit config.nml in \<SIMULATION NAME\> to your liking.
    -   See below for details on config.nml entries
4)  If not already done so, open a MATLAB screen with "screen -S mat"
5)  In MATLAB screen, run "matlab -nodisplay"
6)  In MATLAB, navigate to \<SIMULATION NAME\> and run "gemini3d.model.setup('.','.')"
7)  If not already done so, open a Gemini screen with "screen -S gem"
8)  In Gemini screen, navigate to LynchK directory
9)  In Gemini screen, run "mpiexec -np 36 gemini/gemini3d/build/gemini.bin public_html/Gemini3D/\<SIMULATION NAME\>"
    -   Make sure that -np even divides nx2 x nx3, e.g. 144 x 216 / 36 = 864 cells per processor
10) In MATLAB screen, navigate to LynchK/Jules/thesis
11) In MATLAB screen, run "process('../../public_html/Gemini3D/\<SIMULATION NAME\>',\<OPTIONS\>)". Options include
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

### precip_BG
- PhiWBG = 1e-5 -------- total energy flux (mW/m^2)
- W0BG = 3e3 -------- characteristic energy (eV)

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
```sh
ssh-keygen -t ed25519 -C "your_email@example.com"
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
clip < ~/.ssh/id_ed25519.pub
```
in github.com/julesvanirsel go to
  -   settings
  -   SSH and GPG keys
  -   New SSH key
enter a Title and paste the key
```sh
git remote set-url origin git@github.com:julesvanirsel/thesis.git
git config --global user.name "your_username"
git config --global user.email "your_email_address@example.com"
```
