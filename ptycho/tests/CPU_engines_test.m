%% TEST TEMPLATE FOR FUNTIONALITY CHECK OF CPU ENGINES 
% 1) call standard template to get fresh settings defaults 
% 2) generate artificial data that should serve as a standart test "sample"
% 3) call GPU engine with different basic functionalities and test if all still works
%   !! THESE TEST ARE ONLY USEFUL TO FIND CRASHES IN THE CODE, QUALITY OF THE RECONSTRUCTION IS NOT EVALUATED !!


%% set shared parameters for all test scripts 
run(fullfile( fileparts(mfilename('fullpath')), 'init_test.m'))

%% general settings
p.   artificial_data_file = 'tests/test_data.m';     % artificial data parameters 
p.   asize = [192 192];                              % size of the reconstruction probe
p.   z = 5; 



%% ENGINES
% External C++ code
if isunix
    % Please notice that you have to force data preparation (force_prepare_h5_files=true) if you have made any changes to 
    % the already prepared data (fmag, fmask, positions, sharing ...). 
    eng = struct();
    eng.  name = 'presolver';
    eng.  method = 'DM';
    eng.  asize_presolve = [192 192]; 
    eng.  number_iterations = 10;              % Total number of iterations
    eng.  probe_regularization = .1;            % Weigth factor for the probe update; 
    eng.  probe_change_start = 1;               % Start updating probe at this iteration number
    eng.  probe_support_radius = 0.8;           % Normalized radius of circular support, = 1 for radius touching the window    
    eng.  pfft_relaxation = .05;                % Relaxation in the Fourier domain projection, = 0  for full projection    
    eng.  single_prec = true;                   % single or double precision
    eng.  threads = 20;                         % number of threads for OMP
    eng.  beamline_nodes = [];                  % beamline nodes for the MPI/OMP hybrid, e.g. ['x12sa-cn-2'; 'x12sa-cn-3'];
    eng.  ra_nodes = 0;                         % number of nodes on ra cluster for the MPI/OMP hybrid; set to 0 for current node
    eng.  caller_suffix = '';                   % suffix for the external reconstruction program
    eng.  reconstruction_program = '';          % specify external reconstruction program that overwrites previous settings, e.g. 'OMP_NUM_THREADS=20 ./ptycho_single_OMP';
    eng.  check_cpu_load = false;                % check if specified nodes are already in use (only x12sa). Disable check if you are sure that the nodes are free.
    eng.  initial_conditions_path = '';         % path of the initial conditions file; default if empty (== prepare_data_path)
    eng.  initial_conditions_file = '';    		% Name of the initial conditions file, default if empty. Do not use ~ in the path
    eng.  measurements_file = '';				% Name of the measurements file, default if empty. Do not use ~ in the path
    eng.  solution_file = '';                   % Name of the solution file, default if empty. Do not use ~ in the path
    eng.  force_prepare_h5_files = 0;           % If true before running the C-code the data h5 file is created and the h5 file with initial object and probe too, regardless of whether it exists. It will use the matlab data preparator. 
    [p, ~] = core.append_engine(p, eng);        % Adds this engine to the reconstruction process
end




if isunix
    % Please notice that you have to force data preparation (force_prepare_h5_files=true) if you have made any changes to 
    % the already prepared data (fmag, fmask, positions, sharing ...). 
    eng = struct();
    eng.  name = 'c_solver';
    eng.  method = 'DM+ML';
    eng.  number_iterations = 10;              % Total number of iterations
    eng.  opt_iter = 10;                       % Iterations for optimization     
    eng.  probe_regularization = .1;            % Weigth factor for the probe update; 
    eng.  probe_change_start = 1;               % Start updating probe at this iteration number
    eng.  probe_support_radius = 0.8;           % Normalized radius of circular support, = 1 for radius touching the window    
    eng.  pfft_relaxation = .05;                % Relaxation in the Fourier domain projection, = 0  for full projection    
    eng.  background = 0;                       % [PARTIALLY IMPLEMENTED (not fully optimized)] Add background to the ML model in form:  |Psi|^2+B, B is in average counts per frame and pixel
    eng.  probe_support_fft = false;            % [PARTIALLY IMPLEMENTED (not fully optimized)] Apply probe support in Fourier space, ! uses model zoneplate settings to estimate support size 

    eng.  N_layer = 1;                          % Number of virtual object layers (slices)
    eng.  delta_z = 0e-6 * ones(1, eng.N_layer-1); % Separation between object slices 
    %eng.  ms_init_ob_fraction = [1 0];
    if eng.  N_layer>1
        p.suffix = [p.suffix '_N' num2str(eng. N_layer)];
        eng.  number_iterations = 0; % highly recommended
    end
    
    eng.  single_prec = true;                   % single or double precision
    eng.  threads = 20;                         % number of threads for OMP
    eng.  beamline_nodes = [];                  % beamline nodes for the MPI/OMP hybrid, e.g. ['x12sa-cn-2'; 'x12sa-cn-3'];
    eng.  ra_nodes = 0;                         % number of nodes on ra cluster for the MPI/OMP hybrid; set to 0 for current node
    eng.  caller_suffix = '';                   % suffix for the external reconstruction program
    eng.  reconstruction_program = '';          % specify external reconstruction program that overwrites previous settings, e.g. 'OMP_NUM_THREADS=20 ./ptycho_single_OMP';
    eng.  check_cpu_load = true;                % check if specified nodes are already in use (only x12sa). Disable check if you are sure that the nodes are free.
    eng.  initial_conditions_path = '';         % path of the initial conditions file; default if empty (== prepare_data_path)
    eng.  initial_conditions_file = '';    		% Name of the initial conditions file, default if empty. Do not use ~ in the path
    eng.  measurements_file = '';				% Name of the measurements file, default if empty. Do not use ~ in the path
    eng.  solution_file = '';                   % Name of the solution file, default if empty. Do not use ~ in the path
    eng.  force_prepare_h5_files = 1;           % If true before running the C-code the data h5 file is created and the h5 file with initial object and probe too, regardless of whether it exists. It will use the matlab data preparator. 
    [p, ~] = core.append_engine(p, eng);        % Adds this engine to the reconstruction process
end

% --------- Matlab engines  ------------- 
if 1
    eng = struct();                 % reset settings for this engine 
    eng. name = 'GPU';    
    eng. use_gpu = false;            % if false, run CPU code, but it will get very slow 
    eng. keep_on_gpu = true;        % keep data + projections on GPU, false is useful for large data if DM is used
    eng. compress_data = true;      % use automatic online memory compression to limit meed of GPU memory
    eng. gpu_id = [];               % default GPU id, [] means choosen by matlab
    eng. check_gpu_load = true;     % check available GPU memory before starting GPU engines 
    
    % general 
    eng. number_iterations = 10;   % number of iterations for selected method 
    %eng. asize_presolve = [196 196]; % crop data to "asize_presolve" size to get low resolution estimate 
    %eng. share_probe = 1;           % Share probe between scans. Can be either a number/boolean or a list of numbers, specifying the probe index; e.g. [1 2 2] to share the probes between the second and third scan. 
    %eng. share_object = 0;          % Share object between scans. Can be either a number/boolean or a list of numbers, specifying the object index; e.g. [1 2 2] to share the objects between the second and third scan. 
    eng. method = '';            % choose GPU solver: DM,  ePIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
    eng. opt_errmetric = 'L1' ;     % optimization likelihood - poisson, L1
    eng. grouping = 100;            % size of processed blocks, larger blocks need more memory but they use GPU more effeciently
                                    % for hPIE, ePIE, MLs methods smaller blocks lead to faster convergence, 
                                    % for MLc the convergence is similar 
                                    % for DM, RAAR is has no effect on convergence
    %eng. probe_modes  = 1;          % Number of coherent modes for probe
    eng. object_change_start = 1;    % Start updating object at this iteration number
    eng. probe_change_start = 1;     % Start updating probe at this iteration number

    % regularizations
    eng. reg_mu = 0;                % Regularization constant ( = 0 for no regularization)
    eng. delta = 0;                 % press values to zero out of the illumination area, usually 1e-2 is enough 
    eng. positivity_constraint_object = 0; % enforce weak positivity in object, usually 1e-2 is already enough 

    eng. apply_multimodal_update = false; % apply all incoherent modes to object, it can cause isses if the modes collect some crap 
    eng. probe_backpropagate = 0;         % backpropagate the probe mask, inf == farfield 
    eng. probe_support_radius = [];      % Normalized radius of circular support, = 1 for radius touching the window    
    eng. probe_support_fft = false;        % assume that there is not illumination intensity out of the central FZP cone 

    % basic recontruction parameters 
    % PIE / ML methods
    eng. beta_object = 1;           % object step size, larger == faster convergence, smaller == more robust, should not exceed 1
    eng. beta_probe = 1;            % probe step size, larger == faster convergence, smaller == more robust, should not exceed 1
    eng. delta_p = 0.1;             % LSQ dumping constant, 0 == no preconditioner, 0.1 is usually safe, 
    eng. momentum = 0.5;              % add momentum term to the MLc method, eng.momentum = multiplication gain for velocity
    % eng. delta_z = [50e-6, 50e-6];  % multilayer ptycho extension 

    % DM
    eng. pfft_relaxation = 0.05;     % Relaxation in the Fourier domain projection, = 0  for full projection 
    eng. probe_regularization = 0.1; % Weight factor for the probe update (inertia)

    
    % ADVANCED OPTIONS   
    % position refinement 
    eng. apply_subpix_shift = false;       % apply FFT-based subpixel shift, important for good position refinement but it is slow
    eng. probe_position_search = inf;      % reconstruct probe positions, from iteration == probe_position_search, assume they have to match geometry model with error less than probe_position_error_max
    eng. probe_position_error_max = 20e-9; % max expected random position error of the stages 
    
    % other extensions 
    eng. background = 0.001;                % average background scattering level, for OMNI values around 0.3 for 100ms, for flOMNI <0.1 per 100ms exposure 
    eng. clean_residua = false;            % remove residua from reconstruction by iterative unwrapping, may result in low spatial freq. artefacts 
    eng. regularize_layers = 0.01;         % 0<R<<1 -> apply regularization on the reconstructed object layers, 0 == no regularization 
    
    % wavefront refinement                  
    eng. probe_fourier_shift_search = inf; % refine farfield position of the beam (ie angle) from iteration == probe_fourier_shift_search
    eng. estimate_NF_distance = inf;       % try to estimate the nearfield propagation distance  
    eng. variable_probe = false;           % Use SVD to account for variable illumination during a single (coupled) scan
    eng. variable_probe_modes = 3;         % OPRP settings , number of SVD modes, apply only for PIE methods 
    eng. variable_probe_smooth = 1;        % OPRP settings , apply assumption of smooth evolution of the OPRP modes -> N is order of polynomial fit used for smoothing, 0 == n
    eng. variable_intensity = false;       % account to changes in probe intensity
    
    % extra analysis
    eng. get_fsc_score = false;         % measure evolution of the Fourier ring correlation during convergence 
    eng. mirror_objects = false;        % mirror objects, useful for 0/180deg scan sharing 

    eng. method = 'DM';            % choose GPU solver: DM,  ePIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
    eng. method = 'MLc';            % choose GPU solver: DM,  ePIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
    eng. method = 'MLs';            % choose GPU solver: DM,  ePIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
    eng. method = 'hPIE';            % choose GPU solver: DM,  ePIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
    eng. method = 'ePIE';            % choose GPU solver: DM,  ePIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process

end
 

% Difference Map (Matlab and MEX)


if false
    eng = struct();
    eng. name = 'DM';
    eng. method = 'mex'; 
    eng. number_iterations = 10;               % Total number of iterations
    eng. probe_change_start = 1;              % Start updating probe at this iteration number
    eng. average_start = 300;                 % Start averaging at this iteration number
    eng. average_interval = 5;                % Number of iterations between reconstruction estimates for average 
    eng. count_bound = 4e-2;                  % Relaxed Fourier projection parameter - average photons of change per pixel (= 0 no relaxation) 
    eng. pfft_relaxation = 0.05 ;               % Relaxation in the Fourier domain projection, = 0  for full projection 
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

    eng. use_mex = [1,1,1]; 
    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
end

if 1
    eng = struct();
    eng. name = 'DM';
    eng. method = 'matlab'; 
    eng. number_iterations = 10;               % Total number of iterations
    eng. probe_change_start = 1;              % Start updating probe at this iteration number
    eng. average_start = 300;                 % Start averaging at this iteration number
    eng. average_interval = 5;                % Number of iterations between reconstruction estimates for average 
    eng. count_bound = 4e-2;                  % Relaxed Fourier projection parameter - average photons of change per pixel (= 0 no relaxation) 
    eng. pfft_relaxation = 0.05 ;               % Relaxation in the Fourier domain projection, = 0  for full projection 
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
    eng. use_mex = [0,0,0]; 
    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
end

% Maximum Likelihood (Matlab)
if 1
    eng = struct();
    eng. name = 'ML';
    eng. method = 'matlab'; 
    eng. opt_errmetric = 'L1';                % Error metric for max likelihood = 'poisson', 'L1' (approx poisson), 'L2' (uniform gaussian noise)
    eng. opt_flags = [1 1];                   % Optimize [object probe]
    eng. opt_iter = 5;                        % Iterations for optimization
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
    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
end

    
run(fullfile( ptycho_path, 'tests/run_test.m'))



% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
%   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
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
%       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089).
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
