%  ______                   _           ______ _           _                  
% (_____ \  _              | |         / _____) |         | |                 
%  _____) )| |_ _   _  ____| |__   ___( (____ | |__  _____| |_   _ _____  ___ 
% |  ____(_   _) | | |/ ___)  _ \ / _ \\____ \|  _ \| ___ | | | | | ___ |/___)
% | |      | |_| |_| ( (___| | | | |_| |____) ) | | | ____| |\ V /| ____|___ |
% |_|       \__)\__  |\____)_| |_|\___(______/|_| |_|_____)\_)\_/ |_____|___/ 
%              (____/                                                         
%% Reconstruction parameters
%  Edit only this section
%  References in ptycho_recons.m
%
% For 2 detector positions
%  Put multiple scan numbers e.g. p.scan_number = [5 6];
%  Multiple centers e.g.  p.ctr = [100 207;120 217]; 
  
% For wrapper
%  Change testmode = false
%  Put a for loop of scans around this whole code

% For OMNI automatic queue
% while 1==1  % around this whole code
% and check p.queue.name = 'filelist' 
% and check  p.queue.path = ['../../specES1/ptycho_reconstruct/'];
% it is also recommended to leave the p.scan_number=[] so that it crashes in 
% case it does not read a scan number from the queue file

% For using external C-code - tested on x12sa-cn-4
%  [linux]$ source setup-environment.sh        (load necessary modules and libraries in linux)

% In case of trouble with HDF5 matlab can be opened without its native HDF5 libraries
%   HDF5_DISABLE_VERSION_CHECK=1 matlab &   
% by now this is not needed anymore but I keep it here for documentation or future use

% local_path = fileparts(mfilename('fullpath')); 
% 
% addpath(fullpath(local_path, 'utils'))

addpath(core.find_base_package)

testmode = true;       % In test mode, lock files are not written and plots are generated every iteration

p = struct();
eng = struct();

%% General
p.   verbose_level = 4;%1+testmode;                            % verbosity for standard output (0-1 for loops, 2-3 for testing and adjustments, >= 4 for debugging)
p.   use_display = [];                                      % global switch for display, if [] then true for verbose > 1

p.   scan_number = [35];                                    % Multiple scan numbers for shared scans


% Geometry
p.   z = [5]; %added 5                                             % Distance from object to detector 
p.   asize = [256 256];                                     % Diffr. patt. array size   
p.   ctr = [515 786];                                       % Diffr. patt. center coordinates (y,x) (empty means middle of the array); e.g. [100 207;100+20 207+10];

p.   prop_regime = 'farfield';                              % propagation regime: nearfield, farfield (default), 
                                                                % nearfield info and tips: 
                                                                % 1) This is supported only by GPU engines. 
                                                                % 2) Give the p.focus_to_sample_distance as accurate as possible, positive numbers are sample downstream of focus. Negative numbers or sample upstream of focus is not tested and probably does not work based on a quick look at the code.
                                                                % 3) Dont use eng.estimate_NF_distance, set it = inf. This distance should only be optimized as last resort when you are completely lost in the experiment geometry.  
                                                                % 4) For nearfield the best input probe has no phase curvature, a good approximation can be obtained by just the sqrt of an empty dataframe, or better by taking a farfield probe, propagating to focus and then take the FFT. 
                                                                % 5) If using nearfield use the DM algorithm with many iterations, LSQML is not very stable for some reason. If you want to do ML after DM, and it gets unstable then hold the probe fixed during LSQML.
                                                                %    We have observed sometimes DM diverges, then increasing a bit eng. pfft_relaxation = 0.1, seemed to help
p.   focus_to_sample_distance = [];                         % sample to focus distance, parameter to be set for nearfield ptychography, otherwise it is ignored. Postive value for sample downstream of focus, for negative values originally crashed, we need to check to see if the implementation is correct in that case.
p.   FP_focal_distance = [];                                %  if nonempty -> assume Fourier ptychography configuration, FP_focal_distance = focal length of objective lens for Fourier Ptychography only,
p.   angular_correction_setup = '';                         % if src_positions=='orchestra' or 'nexusCXS' with orchestra files, choose angular correction for specific cSAXS experiment: 'flomni', 'omny', 'lamni', 'none', 
p.   energy = [];                                           % Energy (in keV), leave empty to use spec entry mokev
p.   sample_rotation_angles = [0,0,0];                      % Offaxis ptychography correction , 3x1 vector rotation around [X,Y,beam] axes in degrees , apply a correction accounting for tilted plane oR the sample and ewald sphere curvature (high NA correction)
                                                            % For example, for a grazing incidence experiment with sample oriented along the (x,z) plane,  rotation around x and incidence angle of rotx=-0.6, if we make this [-90-0.6,0,0] it is neither needed to correct the sign of the scanning z direction, nor to flip the diffraction data upside down. 

%p.   affine_angle = 0;                                     % Not used by ptycho_recons at all. This allows you to define a variable for the affine matrix below and keep it in p for future record. This is used later by the affine_matrix_search.m script
p.   affine_matrix = [];                                    % Applies affine transformation (e.g. rotation, stretching) to the positions (ignore by = []). Convention [yn;xn] = M*[y;x]. For flOMNI we found in June 2019: = [1 , 0.0003583 ; 5.811e-05 , 1 ]; for OMNY we found in October 2018: = [1 0;tan(0.4*pi/180) 1]; laMNI in Jul 2019 [1 , -0.006421 ; 0.003766 , 1 ] 

% Scan meta data
p.   src_metadata = 'nexusCXS';                             % source of the meta data, following options are supported: 'spec', 'none' , 'artificial' - or add new to +scan/+meta/

% Scan queue
p.   queue.name = '';                                       % specify file queue; currently only 'filelist' is supported
p.   queue.path=['../../specES1/ptycho_reconstruct/'];      % Folder where the queue of files is defined, note the content of files can overwrite some parameters in p-structure
p.   queue.max_attempts = 5;                                % Max number of attempts to reconstruct a scan.
p.   queue.file_queue_timeout = 10;                         % Time to wait when queue is empty before checking it again 
p.   queue.remote_recons = false;                           % divide the reconstruction into primary/replica processes to reconstruction on a remote server
p.   queue.recon_latest_first = 1;                          % When using 'p.queue_path', (1) reconstruct the latest measurement first or (0) reconstruct in lexicographical order
p.   queue.remote_path = '';                                % Queue list for remote reconstructions. Needs to be accessible for primary and replica processes
p.   queue.tmp_dir_remote = '';                             % shared directory for storing the remote reconstruction

p.   queue.lockfile = false;                                % If true writes a lock file, if lock file exists skips recontruction
p.   spec.waitforscanfinish = true;                         % Checks spec file for the scan end flag 'X#'
p.   spec.check_nextscan_started = true;                    % Waits until the next scan starts to begin reconstructing this one. It is important for OMNY scans with orchestra
p.   spec.isptycho = {};                                    % Use only when SPEC is used: = {'round_roi','cont_line','ura_mesh'}  ( = {} to skip)  List of ptycho spec commands for valid ptycho scans

% Data preparation
p.   detector.name = 'eiger1p5M';                           % see +detectors/ folder 
p.   detector.check_2_detpos = [];                          % = []; (ignores)   = 270; compares to dettrx to see if p.ctr should be reversed (for OMNY shared scans 1221122), make equal to the middle point of dettrx between the 2 detector positions
p.   detector.data_prefix = '';                             % Default using current eaccount e.g. e14169_1_
p.   detector.binning = false;                              % = true to perform 2x2 binning of detector pixels, for binning = N do 2^Nx2^N binning
p.   detector.upsampling = false;                           % upsample the measured data by 2^data_upsampling, (transposed operator to the binning), it can be used for superresolution in nearfield ptychography or to account for undersampling in a far-field dataset
p.   detector.burst_frames = 1;                             % number of frames collected per scan position

p.   prepare.auto_prepare_data = true;                      % if true: prepare dataset from raw measurements if the prepared data does not exist
p.   prepare.force_preparation_data = true;                 % Prepare dataset even if it exists, it will overwrite the file % Default: @prepare_data_2d
p.   prepare.store_prepared_data = true;                    % store the loaded data to h5 even for non-external engines (i.e. other than c_solver)
p.   prepare.prepare_data_function = '';                    % (used only if data should be prepared) custom data preparation function handle;
p.   prepare.auto_center_data = false;                      % if matlab data preparator is used, try to automatically center the diffraction pattern to keep center of mass in center of diffraction

p.   prealign_FP = false;                                   % use prealignment routines for Fourier Ptychography
p.   prealign.asize = [1000 1000];                          % array size for the alignment procedure
p.   prealign.crop_dft = 100;                               % crop the dftregistration input
p.   prealign.prealign_data = true;                         % recalculate the alignment
p.   prealign.axis = 1;                                     % alignment axis
p.   prealign.type = {'round'};                             % alignment routine
p.   prealign.numiter = 5;                                  % number of alignment iterations
p.   prealign.rad_filt_min = 25e-6;                         % discard positions < rad_filt_min radius
p.   prealign.rad_filt_max = 80e-6;                         % discard positions > rad_filt_max radius
p.   prealign.load_alignment = true;                        % load alignment from an alignment_file
p.   prealign.alignment_file = 'alignment_S00249.mat';      % alignment file
p.   prealign.mag_est = 160;                                % estimated magnification; used as an initial guess for the distortion correction matrix
p.   prealign.use_distortion_corr = true;                   % use distortion correction; if distortion_corr is empty, it will calculate a new correction based on the shifts retrieved from the alignment
p.   prealign.distortion_corr = [];                         % distortion correction matrix; [161.3003, 3.4321, -6.7294, 0.0000, 0.9675, 2.0220, 0.0540];
p.   prealign.static_mask = true;                           % use a static mask; helpful for masking the direct beam / beam stop
p.   prealign.static_mask_filename = './staticMaskFP.mat';  % saved static mask; if the file does not exist, a new one will be created


% Scan positions
p.   src_positions = 'nexusCXS';                           % 'spec', 'orchestra', 'load_from_file' or add new position loaders to +scan/+positions/

p.   orchestra.positions_file = ['../../specES1/scan_positions/scan_%05d.dat'];    % Filename pattern for position files, Example: ['../../specES1/scan_positions/scan_%05d.dat']; (the scan number will be automatically filled in) leave empty '' if no orchestra file
p.   spec.motor.fine_motors = {};                           % Y and X motor name for positions, leave empty for defaults
p.   spec.motor.fine_motors_scale = [];                     % ptycho expects real positions in m; 
p.   spec.motor.coarse_motors = {};                         % Coarse sample position for shared object, use {X-motor, Y-motor} 
p.   spec.motor.coarse_motors_scale = [];                   % Scale of the coarse motors (to scale the provided values to meters)

% I/O
p.   prefix = '';                                           % For automatic output filenames. If empty: scan number
p.   suffix = 'run_1';                                      % Optional suffix for reconstruction 
p.   scan_string_format = 'S%05d';                          % format for scan string generation, it is used e.g for plotting and data saving 

p.   base_path = '/users/stud/haffnerm/Ptychoshelves';                                  % base path : used for automatic generation of other paths 
p.   specfile = '';                                         % Name of spec file to get motor positions and check end of scan, defaut is p.spec_file == p.base_path;
p.   nexus_directory = 'data';                              % subdirectory of raw_data_path
p.   ptycho_package_path = '/users/stud/haffnerm/Ptychoshelves/fold_slice'; % ptycho package path. If empty it uses fullfile(p.base_path, 'cxs_software/base/')
p.   base_package_path = '/users/stud/haffnerm/Ptychoshelves/fold_slice'; % base package path. If empty it uses fullfile(p.base_path, 'cxs_software/ptycho/')
p.   raw_data_path{1} = '';                                 % Default using compile_x12sa_filename, used only if data should be prepared automatically
p.   prepare_data_path = '';                                % Default: base_path + 'analysis'. Other example: '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/'; also supports %u to insert the scan number at a later point (e.g. '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/S%.5u')
p.   prepare_data_filename = [];                            % Leave empty for default file name generation, otherwise use [sprintf('S%05d_data_%03dx%03d',p.scan_number(1), p.asize(1), p.asize(2)) p.prep_data_suffix '.h5'] as default 
p.   save_path{1} = '';                                     % Default: base_path + 'analysis'. Other example: '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/'; also supports %u to insert the scan number at a later point (e.g. '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/S%.5u')
p.   io.default_mask_file = '';                             % load detector mask defined in this file instead of the mask in the detector packages, (used only if data should be prepared) 
p.   io.default_mask_type = 'binary';                       % (used only if data should be prepared) ['binary', 'indices']. Default: 'binary' 
p.   io.file_compression = 0;                               % reconstruction file compression for HDF5 files; 0 for no compression
p.   io.data_compression = 1;                               % prepared data file compression for HDF5 files; 0 for no compression
p.   io.load_prep_pos = false;                              % load positions from prepared data file and ignore positions provided by metadata

p.   io.phone_number = [];                                  % phone number for sending messages
p.   io.send_failed_scans_SMS = true;                       % send message if p.queue_max_attempts is exceeded
p.   io.send_finished_recon_SMS = false;                    % send message after the reconstruction is completed
p.   io.send_crashed_recon_SMS = false;                     % send message if the reconstruction crashes
p.   io.SMS_sleep = 1800;                                   % max 1 message per SMS_sleep seconds

p.   artificial_data_file = 'template_artificial_data';     % artificial data parameters, set p.src_metadata = 'artificial' to use this template

%% Reconstruction

% Initial iterate object
p.   model_object = true;                                   % Use model object, if false load it from file 
p.   model.object_type = 'rand';                            % specify how the object shall be created; use 'rand' for a random initial guess; use 'amplitude' for an initial guess based on the prepared data

p.   initial_iterate_object_file{1} = '';                   %  use this mat-file as initial guess of object, it is possible to use wild characters and pattern filling, example: '../analysis/S%05i/wrap_*_1024x1024_1_recons*'



% Initial iterate probe
p.   model_probe = true;                                    % Use model probe, if false load it from file 
p.   model.probe_is_focused = true;                         % Model probe is focused (false: just a pinhole)
p.   model.probe_central_stop = true;                       % Model central stop
p.   model.probe_diameter = 170e-6;                         % Model probe pupil diameter
p.   model.probe_central_stop_diameter = 50e-6;             % Model central stop diameter
p.   model.probe_zone_plate_diameter = 170e-6;              % Model probe zone plate diameter
p.   model.probe_outer_zone_width = [];                     % Model probe zone plate outermost zone width (not used if not a focused probe) 
p.   model.probe_propagation_dist = 1.2e-3;                 % Model probe propagation distance (pinhole <-> sample for unfocused, focal-plane <-> sample for focused)
p.   model.probe_focal_length = 51e-3;                      % Model probe focal length (used only if model_is_focused is true
                                                            %   AND model_outer_zone_width is empty)
p.   model.probe_upsample = 10;                             % Model probe upsample factor (for focused probes)

p.   initial_probe_file = '';                               % Use probe from this h5 reconstruction file (not used if model_probe is true)
p.   probe_file_propagation = 0;                            % Distance for propagating the probe from file in meters, = 0 to ignore

p.   model.probe_zone_plate_diameter_for_propagation = 170e-6; % Parameter to calculate the propagation with energy
p.   prop_probe_with_energy = false;                         % Calculate the probe propagation if current and initial probe energy is different, used for spectroscopy scan.
%  If you have spectral scans, without the FZP moving, this propagation distance compensates for the change of focal length with energy. Useful in the 100-200 eV range
%  The calculation needs p.model.probe_zone_plate_diameter_for_propagation and p.model.probe_outer_zone_width
%  To enable set p.prop_probe_with_energy = true;

% Shared scans - Currently working only for sharing probe and object
p.   share_probe = 1;                                       % Share probe between scans. Can be either a number/boolean or a list of numbers, specifying the probe index; e.g. [1 2 2] to share the probes between the second and third scan. 
p.   share_object = 0;                                      % Share object between scans. Can be either a number/boolean or a list of numbers, specifying the object index; e.g. [1 2 2] to share the objects between the second and third scan. 

% Modes
p.   probe_modes  = 1;                                      % Number of coherent modes for probe
p.   object_modes = 1;                                      % Number of coherent modes for object
% Mode starting guess
p.   mode_start_pow = [0.02];                               % Normalized intensity on probe modes > 1. Can be a number (all higher modes equal) or a vector
p.   mode_start = 'herm';                                   % (for probe) = 'rand', = 'herm' (Hermitian-like base), = 'hermver' (vertical modes only), = 'hermhor' (horizontal modes only)
p.   ortho_probes = true;                                   % orthogonalize probes after each engine




%% Plot, save and analyze

p.   plot.prepared_data = testmode;                         % plot prepared data
p.   plot.interval = [];                                    % plot each interval-th iteration, does not work for c_solver code
p.   plot.log_scale = [0 0];                                % Plot on log scale for x and y
p.   plot.realaxes = true;                                  % Plots show scale in microns
p.   plot.remove_phase_ramp = false;                        % Remove phase ramp from the plotted / saved phase figures 
p.   plot.fov_box = true;                                   % Plot the scanning FOV box on the object (both phase and amplitude)
p.   plot.fov_box_color = 'r';                              % Color of the scanning FOV box
p.   plot.positions = true;                                 % Plot the scanning positions
p.   plot.mask_bool = true;                                 % Mask the noisy contour of the reconstructed object in plots
p.   plot.windowautopos = true;                             % First plotting will auto position windows
p.   plot.obj_apod = false;                                 % Apply apodization to the reconstructed object;
p.   plot.prop_obj = 0;                                     % Distance to propagate reconstructed object before plotting [m]
p.   plot.show_layers = true;                               % show each layer in multilayer reconstruction 
p.   plot.show_layers_stack = false;                        % show each layer in multilayer reconstruction by imagesc3D
p.   plot.object_spectrum = true;                             % Plot propagated object (FFT for conventional ptycho); if empty then default is false if verbose_level < 3 and true otherwise
p.   plot.probe_spectrum = true;                              % Plot propagated probe (FFT for conventional ptycho); if empty then default is false if verbose_level < 3 and true otherwise
p.   plot.conjugate = false;                                % plot complex conjugate of the reconstruction 
p.   plot.horz_fact = 2.5;                                  % Scales the space that the ptycho figures take horizontally
p.   plot.FP_maskdim = 180e-6;                              % Filter the backpropagation (Fourier Ptychography)
p.   plot.calc_FSC = false;                                 % Calculate the Fourier Shell correlation for 2 scans or compare with model in case of artificial data tests 
p.   plot.show_FSC = false;                                 % Show the FSC plots, including the cropped FOV
p.   plot.residua = false;                                  % highlight phase-residua in the image of the reconstructed phase

p.   save.external = ~testmode;                             % Use a new Matlab session to run save final figures (saves ~6s per reconstruction). Please be aware that this might lead to an accumulation of Matlab sessions if your single reconstruction is very fast.
p.   save.store_images = true;                              % Write preview images containing the final reconstructions in [p.base_path,'analysis/online/ptycho/'] if p.use_display = 0 then the figures are opened invisible in order to create the nice layout. It writes images in analysis/online/ptycho
p.   save.store_images_intermediate = false;                % save images to disk after each engine
p.   save.store_images_ids = 1:4;                           % identifiers  of the figure to be stored, 1=obj. amplitude, 2=obj. phase, 3=probes, 4=errors, 5=probes spectrum, 6=object spectrum
p.   save.store_images_format = 'png';                      % data type of the stored images jpg or png 
p.   save.store_images_dpi = 150;                           % DPI of the stored bitmap images 
p.   save.exclude = {'fmag', 'fmask', 'illum_sum'};         % exclude variables to reduce the file size on disk
p.   save.save_reconstructions_intermediate = false;        % save final object and probes after each engine
p.   save.save_reconstructions = true;                      % save reconstructed object and probe when full reconstruction is finished 
p.   save.output_file = 'h5';                               % data type of reconstruction file; 'h5' or 'mat'

% Was missing for YJ's version
p.fout = '/users/stud/haffnerm/Ptychoshelves/temp_out';


%% ENGINES
% External C++ code
if 0
    % Please notice that you have to force data preparation (force_prepare_h5_files=true) if you have made any changes to 
    % the already prepared data (fmag, fmask, positions, sharing ...). 
    eng.  name = 'presolver';                   % A wrapper around c_solver engine which uses DM to presolve reconstruction first in low resolution. It can be used to refine the lowest spatial frequencies in a computationally cheaper way
    eng.  asize_presolve = [256 256];           % Diffraction size (p.asize value) for the downsampled dataset, IMPORTANT: this is asize before binning
    eng.  number_iterations = 500;              % Total number of iterations of DM method 
    eng.  probe_regularization = .1;            % DM: Weight factor for the probe update; 
    eng.  probe_change_start = 1;               % DM: Start updating probe at this iteration number
    eng.  probe_support_radius = 0.8;           % DM:   Normalized radius of circular support, = 1 for radius touching the window    
    eng.  pfft_relaxation = .05;                % DM: Relaxation in the Fourier domain projection, = 0  for full projection    
    
    % advanced setting unrelated to reconstruction quality 
    eng.  single_prec = true;                   % single or double precision
    eng.  threads = 20;                         % number of threads for OMP
    eng.  beamline_nodes = [];                  % cSAXS only option: beamline nodes for the MPI/OMP hybrid, e.g. ['x12sa-cn-2'; 'x12sa-cn-3'];
    eng.  use_gpu = false;                      % use GPU if available
    eng.  ra_nodes = 2;                         % PSI only option:   number of nodes on ra cluster for the MPI/OMP hybrid; set to 0 for current node
    eng.  caller_suffix = '';                   % suffix for the external reconstruction program
    eng.  reconstruction_program = '';          % specify external reconstruction program that overwrites previous settings, e.g. 'OMP_NUM_THREADS=20 ./ptycho_single_OMP';
    eng.  check_cpu_load = true;                % check if specified nodes are already in use (only x12sa). Disable check if you are sure that the nodes are free.
    eng.  initial_conditions_path = '';         % path of the initial conditions file; default if empty (== prepare_data_path)
    eng.  initial_conditions_file = '';    		% Name of the initial conditions file, default if empty. 
    eng.  measurements_file = '';				% Name of the measurements file, default if empty. 
    eng.  solution_file = '';                   % Name of the solution file, default if empty. 
    eng.  force_prepare_h5_files = false;       % If true before running the C-code the data h5 file is created and the h5 file with initial object and probe too, regardless of whether it exists. It will use the matlab data preparator. 
    [p, ~] = core.append_engine(p, eng);        % Adds this engine to the reconstruction process
end

if 1
    % Please notice that you have to force data preparation (force_prepare_h5_files=true) if you have made any changes to 
    % the already prepared data (fmag, fmask, positions, sharing ...). 
    
    eng.  name = 'c_solver';                    % Fast parallelized CPU code that combines DM and ML engines 
    eng.  number_iterations = 300;              % Iterations for Difference maps
    eng.  opt_iter = 100;                       % Iterations for maximum likelihood optimization, ML is automatically ended when no further improvement is found 
    eng.  probe_regularization = .1;            % DM only: Weight factor for the probe update; 
    eng.  probe_change_start = 1;               % DM only: Start updating probe at this iteration number
    eng.  probe_support_radius = 0.8;           % DM+ML:   Normalized radius of circular support, = 1 for radius touching the window    
    eng.  pfft_relaxation = .05;                % DM only: Relaxation in the Fourier domain projection, = 0  for full projection    
    eng.  background = 0;                       % ML only: [PARTIALLY IMPLEMENTED (not fully optimized)] Add background to the ML model in form:  |Psi|^2+B, B is in average counts per frame and pixel
    eng.  probe_support_fft = false;            % ML only: Apply probe support in Fourier space,   uses p.model.probe_outer_zone_width  to estimate the support size 

    % multilayer ptychography
    Nlayers = 1; 
    eng.  delta_z = 0e-6 * ones(1, Nlayers-1); % Separation between object slices 
    if  Nlayers>1
        p.suffix = [p.suffix '_N' num2str(Nlayers)];
        eng.  number_iterations = 0; % highly recommended to skip the DM engine
    end
    eng. preshift_ML_probe = true;             % multilayer extension: if true, assume that the provided probe is reconstructed in center of the sample and the layers are centered around this position 

    % advanced setting unrelated to reconstruction quality 
    eng.  single_prec = true;                   % single or double precision
    eng.  threads = 20;                         % number of threads for OMP
    eng.  beamline_nodes = {};                  % beamline nodes for the MPI/OMP hybrid, e.g. {'x12sa-cn-2', 'x12sa-cn-3'};
    eng.  use_gpu = true;                       % use GPU if available
    eng.  gpu_id = 1;                           % GPU device number; can use   utils.find_unused_gpu_card(4);   to automatically check; or check with nvidia-smi (+1 for Matlab indexing)
    eng.  num_gpus = 1;                         % number of GPUs
    eng.  ra_nodes = 1;                         % number of nodes on ra cluster for the MPI/OMP hybrid; set to 0 for current node
    eng.  ra_reservation = '';                  % reservation flag for reserved ra nodes, example 'p18532_12Oct'
    eng.  slurm_partition = '';                 % Normally uses 'day' or 'week' depending on number of iterations, but you can specify another one, e.g. 'rhel77_test'
    eng.  caller_suffix = '';                   % suffix for the external reconstruction program
    eng.  reconstruction_program = '';          % specify external reconstruction program that overwrites previous settings, e.g. 'OMP_NUM_THREADS=20 ./ptycho_single_OMP';
    eng.  check_cpu_load = true;                % check if specified nodes are already in use (only x12sa). Disable check if you are sure that the nodes are free.
    eng.  initial_conditions_path = '';         % path of the initial conditions file; default if empty (== prepare_data_path)
    eng.  initial_conditions_file = '';    		% Name of the initial conditions file, default if empty. 
    eng.  measurements_file = '';				% Name of the measurements file, default if empty. 
    eng.  solution_file = '';                   % Name of the solution file, default if empty. 
    eng.  force_prepare_h5_files = false;       % If true before running the C-code the data h5 file is created and the h5 file with initial object and probe too, regardless of whether it exists. It will use the matlab data preparator. 
    [p, ~] = core.append_engine(p, eng);        % Adds this engine to the reconstruction process
end

% Difference Map (Matlab and MEX)               See for more details: Thibault, Pierre, et al. "High-resolution scanning x-ray diffraction microscopy." Science 321.5887 (2008): 379-382.
if 0
    eng. name = 'DM';
    eng. number_iterations = 300;             % Total number of iterations
    eng. probe_change_start = 1;              % Start updating probe at this iteration number
    eng. average_start = 300;                 % Start averaging at this iteration number
    eng. average_interval = 5;                % Number of iterations between reconstruction estimates for average 
    eng. count_bound = 4e-2;                  % Relaxed Fourier projection parameter - average photons of change per pixel (= 0 no relaxation) 
    eng. pfft_relaxation = 0.05               % Relaxation in the Fourier domain projection, = 0  for full projection 
    eng. probe_regularization = .1;           % Weigth factor for the probe update
    eng. probe_mask_bool = true;              % If true, impose a support constraint to the probe
    eng. probe_mask_area = .9;                % Area ratio of the mask
    eng. probe_mask_use_auto = false;         % Use autocorrelation for probe_mask (if false: circular circle)
    eng. object_flat_region = [];             % Mask for enforcing a flat region in the object (to reduce artifacts)
    eng. remove_scaling_ambiguity = true;     % Remove ambiguity of the probe times object scalling by probe normalization
    eng. clip_object = true;                  % Clip the object transmission function
    eng. clip_max = 1.0;                      % Upper bound
    eng. clip_min = 0.0;                      % Lower bound
    eng. compute_rfact = false;               % If set to true, R-factor is computed at every iteration (large overhead!!!)

    [p, ~] = core.append_engine(p, eng);      % Adds this engine to the reconstruction process
end

% Maximum Likelihood (Matlab)                 See for more details:  Thibault, P., and M. Guizar-Sicairos, New Journal of Physics 14.6 (2012): 063004.
if 0
    eng. name = 'ML';
    eng. opt_errmetric = 'L1';                % Error metric for max likelihood = 'poisson', 'L1' (approx poisson), 'L2' (uniform gaussian noise)
    eng. opt_flags = [1 1];                   % Optimize [object probe]
    eng. opt_iter = 100;                      % Iterations for optimization
    eng. opt_ftol = 1e-10;                    % Tolerance on error metric for optimization
    eng. opt_xtol = 1e-7;                     % Tolerance on optimizable parameters
    eng. probe_mask_bool = true;              % If true, impose a support constraint to the probe
    eng. probe_mask_area = .9;                % Area ratio of the mask
    eng. probe_mask_use_auto = false;         % Use autocorrelation for probe_mask (if false: circular circle)
    eng. scale_gradient = false;              % Preconditioning by scaling probe gradient - Reported useful for weak objects
    eng. inv_intensity = false;               % Make error metric insensitive to intensity fluctuations
    eng. use_probe_support = false;           % Use the support on the probe that was used in the difference-map
    eng. reg_mu = 0.01;                       % Regularization constant ( = 0 for no regularization)
    eng. smooth_gradient = true;              % Sieves preconditioning, =false no smoothing, = true uses Hanning, otherwise specify a small matrix making sure its sum = 1
    [p, ~] = core.append_engine(p, eng);      % Adds this engine to the reconstruction process
end

% 3D Maximum Likelihood (Matlab)             See for more details: Tsai EH, et al, Optics express. 2016 Dec 12;24(25):29089-108.
if 0
    eng. name = 'ML_MS';
    eng. ms_opt_iter = 200; 
    eng. N_layer = 2;
    eng. delta_z = [50]*1e-6 * ones(1, eng.N_layer-1);
    eng. ms_init_ob_fraction = [];            % Default: 1/N_layer    
    eng. ms_opt_flags = [1 1 0];              % Optimize [object, probe, and separation (delat_z)]
    eng. ms_opt_z_param = [200 50];           % [Every this iteratsion, run this many iterations to optimize z (likely to converge earlier anyway)]
    eng. ms_grado_roi = [];                   % Consider using scan_roi   

    
    eng. opt_errmetric = 'L1';                % Error metric for max likelihood = 'poisson', 'L1' (approx poisson), 'L2' (uniform gaussian noise)
    eng. opt_ftol = 1e-10;                    % Tolerance on error metric for optimization
    eng. opt_xtol = 1e-7;                     % Tolerance on optimizable parameters
    eng. probe_mask_bool = true;              % If true, impose a support constraint to the probe
    eng. probe_mask_area = .9;                % Area ratio of the mask
    eng. probe_mask_use_auto = false;         % Use autocorrelation for probe_mask (if false: circular circle)
    eng. scale_gradient = false;              % Preconditioning by scaling probe gradient - Reported useful for weak objects
    eng. inv_intensity = false;               % Make error metric insensitive to intensity fluctuations
    eng. use_probe_support = false;           % Use the support on the probe that was used in the difference-map
    eng. reg_mu = 0.01;  %0.01                % Regularization constant ( = 0 for no regularization)
    eng. smooth_gradient = true;              % Sieves preconditioning, =false no smoothing, = true uses Hanning, otherwise specify a small matrix making sure its sum = 1
    [p, ~] = core.append_engine(p, eng);      % Adds this engine to the reconstruction process
end

% 3D Difference Map (Matlab)                See for more details: Tsai EH, et al, Optics express. 2016 Dec 12;24(25):29089-108.
if 0
    eng. name = 'DM_MS';
    eng. number_iterations = 20;               % Total number of iterations
    eng. N_layer = 2;
    eng. delta_z = [50]*1e-6 * ones(1, eng.N_layer-1);
    eng. ms_init_ob_fraction = [];  % Default: 1/N_layer    
    % p.suffix = [p.suffix '_N' num2str(eng. N_layer)]; 
    
    eng. ms_reverse_order_iter = [];  % Run the additional reverse order update for these iterations
    eng. ratio_reverse = 0.5;         
    eng. ms_scan_roi = []*1e-6;  % (Not implemented yet!)     
    eng. use_mex = zeros(1,3);
    eng. opt_errmetric = 'L1';
    eng. inv_intensity = false; 
    
    eng. probe_change_start = 1;              % Start updating probe at this iteration number
    eng. average_start = 300;                 % Start averaging at this iteration number
    eng. average_interval = 5;                % Number of iterations between reconstruction estimates for average 
    eng. count_bound = 4e-2;                  % Relaxed Fourier projection parameter - average photons of change per pixel (= 0 no relaxation) 
    eng. probe_regularization = .1;           % Weigth factor for the probe update
    eng. probe_mask_bool = true;              % If true, impose a support constraint to the probe
    eng. probe_mask_area = .9;                % Area ratio of the mask
    eng. probe_mask_use_auto = false;         % Use autocorrelation for probe_mask (if false: circular circle)
    eng. object_flat_region = [];             % Mask for enforcing a flat region in the object (to reduce artifacts)
    [p, ~] = core.append_engine(p, eng);      % Adds this engine to the reconstruction process
end

% --------- GPU engines  -------------   See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
if 0 
    eng = struct();                        % reset settings for this engine 
    eng. name = 'GPU';    
    eng. use_gpu = true;                   % if false, run CPU code, but it will get very slow 
    eng. keep_on_gpu = true;               % keep data + projections on GPU, false is useful for large data if DM is used
    eng. compress_data = true;             % use automatic online memory compression to limit need of GPU memory
    eng. gpu_id = [];                      % default GPU id, [] means choosen by matlab
    eng. check_gpu_load = true;            % check available GPU memory before starting GPU engines 
    
    % general 
    eng. number_iterations = 300;          % number of iterations for selected method 
    %eng. asize_presolve = [196 196];      % crop data to "asize_presolve" size to get low resolution estimate that can be used in the next engine as a good initial guess 
    %eng. share_probe = 1;                 % Share probe between scans. Can be either a number/boolean or a list of numbers, specifying the probe index; e.g. [1 2 2] to share the probes between the second and third scan. 
    %eng. share_object = 0;                % Share object between scans. Can be either a number/boolean or a list of numbers, specifying the object index; e.g. [1 2 2] to share the objects between the second and third scan. 
    eng. align_shared_objects = false;     % before merging multiple unshared objects into one shared, the object will be aligned and the probes shifted by the same distance -> use for alignement and shared reconstruction of drifting scans  

    eng. method = 'MLc';                   % choose GPU solver: DM, ePIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
    eng. opt_errmetric = 'L1' ;            % optimization likelihood - poisson, L1
    eng. grouping = inf;                   % size of processed blocks, larger blocks need more memory but they use GPU more effeciently, !!! grouping == inf means use as large as possible to fit into memory 
                                           % * for hPIE, ePIE, MLs methods smaller blocks lead to faster convergence, 
                                           % * for MLc the convergence is similar 
                                           % * for DM is has no effect on convergence
    %eng. probe_modes  = 1;                % Number of coherent modes for probe
    eng. object_change_start = 1;          % Start updating object at this iteration number
    eng. probe_change_start = 1;           % Start updating probe at this iteration number

    % regularizations
    eng. reg_mu = 0;                       % Regularization (smooting) constant ( reg_mu = 0 for no regularization)
    eng. delta = 0;                        % press values to zero out of the illumination area in th object, usually 1e-2 is enough 
    eng. positivity_constraint_object = 0; % enforce weak (relaxed) positivity in object, ie O = O*(1-a)+a*|O|, usually a=1e-2 is already enough. Useful in conbination with OPRP or probe_fourier_shift_search  

    eng. apply_multimodal_update = false;   % apply all incoherent modes to object. Sometimes the higher probe modes account for some problems in the measurements, such as background, in those cases is important to set to false so that these modes do not interact with the object. 
    eng. probe_backpropagate = 0;           % backpropagation distance the probe mask, 0 == apply in the object plane. Useful for pinhole imaging where the support can be applied  at the pinhole plane
    eng. probe_support_radius = [];         % Normalized radius of circular support, = 1 for radius touching the window    
    eng. probe_support_fft = false;         % assume that there is not illumination intensity out of the central FZP cone
    eng. probe_support_fft_type = 'circ';   % choose the type of Fourier support; select 'circ' for a circular support based on FZP parameters, 'gaps' for detector gaps and 'mask' for the detector mask
    eng. probe_support_fft_float = false;   % apply only a weak probe support, a relaxed projection, instead of setting the values to 0. If = 0 or false then probe support projections are normal and amplitude outside support is set to zero. If e.g. = 0.1 then leaves 10 percent of amplitude outside the support.
    eng. probe_support_fft_relax = 0.9;     % probe support relaxation, 1 for full projection


    % basic recontruction parameters 
    % PIE / ML methods                    % See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
    eng. beta_object = 1;                 % object step size, larger == faster convergence, smaller == more robust, should not exceed 1
    eng. beta_probe = 1;                  % probe step size, larger == faster convergence, smaller == more robust, should not exceed 1
    eng. delta_p = 0.1;                   % LSQ dumping constant, 0 == no preconditioner, 0.1 is usually safe, Preconditioner accelerates convergence and ML methods become approximations of the second order solvers 
    eng. momentum = 0;                    % add momentum acceleration term to the MLc method, useful if the probe guess is very poor or for acceleration of multilayer solver, but it is quite computationally expensive to be used in conventional ptycho without any refinement. The momentum method works usually well even with the accelerated_gradients option.  eng.momentum = multiplication gain for velocity, eng.momentum == 0 -> no acceleration, eng.momentum == 0.5 is a good value
    eng. accelerated_gradients_start = 2; % iteration number from which the Nesterov gradient acceleration should be applied, this option is supported only for MLc method. It is very computationally cheap way of convergence acceleration. 

    % DM
    eng. pfft_relaxation = 0.05;          % Relaxation in the Fourier domain projection, = 0  for full projection 
    eng. probe_regularization = 0.1;      % Weight factor for the probe update (inertia)

    
    % ADVANCED OPTIONS                     See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
    % position refinement 
    eng. apply_subpix_shift = false;       % apply FFT-based subpixel shift, it is automatically allowed for position refinement
    eng. probe_position_search = inf;      % iteration number from which the engine will reconstruct probe positions, from iteration == probe_position_search, assume they have to match geometry model with error less than probe_position_error_max
    eng. probe_geometry_model = {'scale', 'asymmetry', 'rotation', 'shear'};  % list of free parameters in the geometry model, choose from: {'scale', 'asymmetry', 'rotation', 'shear'}
    eng. probe_position_error_max = 20e-9; % maximal expected random position errors, probe prositions are confined in a circle with radius defined by probe_position_error_max and with center defined by original positions scaled by probe_geometry_model

    % multilayer extension 
    eng. delta_z = [];                     % if not empty, use multilayer ptycho extension , see ML_MS code for example of use, [] == common single layer ptychography , note that delta_z provides only relative propagation distance from the previous layer, ie delta_z can be either positive or negative. If preshift_ML_probe == false, the first layer is defined by position of initial probe plane. It is useful to use eng.momentum for convergence acceleration 
    eng. regularize_layers = 0;            % multilayer extension: 0<R<<1 -> apply regularization on the reconstructed object layers, 0 == no regularization, 0.01 == weak regularization that will slowly symmetrize information content between layers 
    eng. preshift_ML_probe = true;         % multilayer extension: if true, assume that the provided probe is reconstructed in center of the sample and the layers are centered around this position 

    % other extensions 
    eng. background = 0.001;               % average background scattering level, for OMNI values around 0.3 for 100ms, for flOMNI <0.1 per 100ms exposure, see for more details: Odstrcil, M., et al., Optics letters 40.23 (2015): 5574-5577.
    eng. background_width = inf;           % width of the background function in pixels,  inf == flat background, background function is then convolved with the average diffraction pattern in order to account for beam diversion 
    eng. clean_residua = false;            % remove phase residua from reconstruction by iterative unwrapping, it will result in low spatial freq. artefacts -> object can be used as an residua-free initial guess for netx engine


    % wavefront & camera geometry refinement     See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
    eng. probe_fourier_shift_search = inf; % iteration number from which the engine will: refine farfield position of the beam (ie angle) from iteration == probe_fourier_shift_search
    eng. estimate_NF_distance = inf;       % iteration number from which the engine will: try to estimate the nearfield propagation distance using gradient descent optimization  
    eng. detector_rotation_search = inf;   % iteration number from which the engine will: search for optimal detector rotation, preferably use with option mirror_scan = true , rotation of the detector axis with respect to the sample axis, similar as rotation option in the position refinement geometry model but works also for 0/180deg rotation shared scans 
    eng. detector_scale_search = inf;      % iteration number from which the engine will: refine pixel scale of the detector, can be used to refine propagation distance in ptycho 
    eng. variable_probe = false;           % Use SVD to account for variable illumination during a single (coupled) scan, see for more details:  Odstrcil, M. et al. Optics express 24.8 (2016): 8360-8369.
    eng. variable_probe_modes = 1;         % OPRP settings , number of SVD modes using to describe the probe evolution. 
    eng. variable_probe_smooth = 0;        % OPRP settings , enforce of smooth evolution of the OPRP modes -> N is order of polynomial fit used for smoothing, 0 == do not apply any smoothing. Smoothing is useful if only a smooth drift is assumed during the ptycho acquisition 
    eng. variable_intensity = false;       % account to changes in probe intensity

    % extra analysis
    eng. get_fsc_score = false;            % measure evolution of the Fourier ring correlation during convergence 
    eng. mirror_objects = false;           % mirror objects, useful for 0/180deg scan sharing -> geometry refinement for tomography, works only if 2 scans are provided 

    % custom data adjustments, useful for offaxis ptychography
    eng.auto_center_data = false;           % autoestimate the center of mass from data and shift the diffraction patterns so that the average center of mass corresponds to center of mass of the provided probe 
    eng.auto_center_probe = false;          % center the probe position in real space before reconstruction is started 
    eng.custom_data_flip = [0,0,0];         % apply custom flip of the data [fliplr, flipud, transpose]  - can be used for quick testing of reconstruction with various flips or for reflection ptychography 
    eng.apply_tilted_plane_correction = ''; % if any(p.sample_rotation_angles([1,2]) ~= 0),  this option will apply tilted plane correction. (a) 'diffraction' apply correction into the data, note that it is valid only for "low NA" illumination  Gardner, D. et al., Optics express 20.17 (2012): 19050-19059. (b) 'propagation' - use tilted plane propagation, (c) '' - will not apply any correction 

    
    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process

end
    
    


%% Run the reconstruction
caller = dbstack;
if length(caller)==1
    tic
    out = core.ptycho_recons(p);
    toc
end

%end
% 2011-11-24
% Parameter to autoposition windows on first display - p.windowautopos
% Replaced powerbound with countbound. countbound represents the mean
  % number of photons in a change below which no projection is taken. It
  % scales automatically with exposure time (number of photons in
  % measurement)
% Real axes option to show plots in microns
% Read parameters from spec
% Implement user suplied object_flat_region
% Implemented option for reconstructing when having 2 repeated scans in the prepared data file

% 2011-11-29
% Template seemed extracted from an AFS run, I modified directories for
% direct use on ../../
% Implemented test mode
% Added cutoff value at beginning
% Added auto settings for prepare data, scan numbers
% Implemented reading from spec. Note it will use the values from the first
% scan
% Added option for repeated scan, should be enabled for 2 detector positions

% 2012-08-23
% Replaced default prepare data function to prepare_data_2d
% In I/O section: added option for a sufix 
% Added default option for raw data path based on compile_x12sa_filename
% Added options to autoprepare data, with cutoff and burstmode detected if
  % the file does not exist. Also added the possiblity to override and
  % force a repreparation of data
% Added a data prefix option (for eaccount_1_) and defaults using
  % identify_eaccount

% 2012-10-29
% Added option for binning and some checks for OMNY detector position scans

% 2012-10-31
% Added options to use the external C-code for testing

% 2015-05-13
% Added option to queue file tasks from OMNI.
% For this I moved the default checks and generation of default names and
% paths to ptycho_recons. Que Dios se apiade de nosotros.

% 2016-02-11
% Removed old option for dump files
% Added p.store_images, if this flag is on and p.use_display it will open 
% figures in the background and write nice jpegs of the reconstruction and error metric anyway

% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the PtychoShelves 
%   computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
%
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the PtychoShelves package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite 
% K. Wakonig, H.-C. Stadler, M. Odstrčil, E.H.R. Tsai, A. Diaz, M. Holler, I. Usov, J. Raabe, A. Menzel, M. Guizar-Sicairos, PtychoShelves, a versatile 
% high-level framework for high-performance analysis of ptychographic data, J. Appl. Cryst. 53(2) (2020). (doi: 10.1107/S1600576720001776)
% and for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for LSQ-ML:
% M. Odstrčil, A. Menzel, and M. Guizar-Sicairos, Iterative least-squares solver for generalized maximum-likelihood ptychography, Opt. Express 26(3), 3108 (2018). 
%   (doi: 10.1364/OE.26.003108),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089),
% and/or for OPRP:
% M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  Ptychographic coherent diffractive imaging with orthogonal probe relaxation. 
% Opt. Express 24.8 (8360-8369) 2016. (doi: 10.1364/OE.24.008360).
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland.