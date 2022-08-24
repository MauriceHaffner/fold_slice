%ARTIFICIAL_DATA parameters
%   Specify simulation parameters and overwrite fields for a more realistic
%   simulation.

%% overwrite parameters in p 
p.   src_positions = 'matlab_pos';                                    % 'spec', 'omny' or 'matlab_pos' (scan params are defined below)
p.   data_prefix = 'artif_data'; 


%% simulation specific paramters
% p.simulation.material{1} = 'H2O';                                     % Simulated material
% p.simulation.material{2} = 'CH3';                                     % Simulated material
% p.simulation.material_density = [-1, 0.95];
% p.simulation.dataset = 7522;                                         % Simulated sample; see +core/+simulations/create_object for more information
% p.simulation.objheight = 0.5e-6;                                    % Object height of the simulated sample

p.simulation.energy = 6.2;                                          % Simulated energy
p.simulation.delta = [];                                          % Simulated refractive index (n = 1-delta+1j*beta); leave empty for CXRO values
p.simulation.beta = [];                                           % Simulated refractive index (n = 1-delta+1j*beta); leave empty for CXRO values
p.simulation.material_density = [-1];
p.simulation.material = {'Au'};                                   % Simulated material
p.simulation.dataset = {5}; % {[1,8], [8,1]};                     % Simulated sample; see +core/+simulations/create_object for more information
p.simulation.rotation_angle = [];                                 % apply rotation angle [degrees] on the 3D phantom  dataset, [] == ignore
p.simulation.flip_objects_180deg = [];                            % Simulated 0 vs 180 deg projection 
p.simulation.thickness = 0e-6;                                    % keep 0 for thin samples simulation , Total thickness of the simulated sample if multiple layers are provided
p.simulation.objheight = 1e-6;                                    % Object height of the simulated sample

p.simulation.incoherence_blur=0;                                  % Incoherence blur
p.simulation.photons_per_pixel = inf; % 50;                             % Photon count; if photons_per_pixel = inf -> no noise 
p.simulation.position_uncertainty = 0;                            % Position uncertainty
p.simulation.apply_sub_px_shifts = true;                                     % Allow sub-pixel shifts
p.simulation.delta_z_total = 300e-6; 
p.simulation.objwidth = 30e-6;
p.simulation.subscans = 1;                                        % multiple "scans" per diffraction dataset 
p.simulation.simulate_detector = true;                            % Simulate detector with proper bad pixel mask
p.simulation.prop_from_focus = 0;                                 % Propagate views before FFT
p.simulation.apply_sub_px_shifts = true;                          % Allow sub-pixel shifts

p.simulation.affine_matrix = [];                                 % Apply affine matrix on the positions before producing data 


% Offaxis ptychography  
p.simulation.sample_rotation_angles = [0,0,0];                   % 3x1 vector rotation around [X,Y,beam] axes in degrees , apply a correction accounting for tilted plane oR the sample and ewald sphere curvature (high NA correction)



% Fourier Ptychography
p.simulation.illumination = 30e-6;                                % Illumination size

% Detector settings
p.simulation.det_pixel_size = 75e-6;                                  %  default pixel size for virtual detector 
                             
                   
% Scan parameters
p.   scan.type = 'fermat';                                       % {'round', 'raster', 'round_roi', 'custom', 'fermat'}
p.   scan.radius_in = 0;                                              % round scan: interior radius of the round scan
p.   scan.radius_out = 5e-6;                                          % round scan: exterior radius of the round scan
p.   scan.nr = 10;                                                    % round scan: number of intervals (# of shells - 1)
p.   scan.nth = 3;                                                    % round scan: number of points in the first shell
p.   scan.lx = 15e-6;                                                 % round_roi scan: width of the roi
p.   scan.ly = 15e-6;                                                 % round_roi scan: height of the roi
p.   scan.dr = 1.5e-6;                                                % round_roi scan: shell step size
p.   scan.nx = 10;                                                    % raster scan: number of steps in x
p.   scan.ny = 10;                                                    % raster scan: number of steps in y
p.   scan.step_size_x = 0.5e-6;                                         % raster scan: step size (grid spacing)
p.   scan.step_size_y = 0.5e-6;                                         % raster scan: step size (grid spacing)
p.   scan.b = 0;                                                      % fermat: angular offset
p.   scan.n_max = 1e4;                                                % fermat: maximal number of points generated 
p.   scan.step = 1e-6;                                                % fermat: step size 
p.   scan.cenxy = [0,0];                                              % fermat: position of center offset 
p.   scan.scan_roi = [];                                              % Region of interest in the object [xmin xmax ymin ymax] in meters. Points outside this region are not used for reconstruction.


% Initial iterate probe
p.   model_probe = false;                                           % Use model probe
% p.   model.probe_is_focused = true;                               % Model probe is focused (false: just a pinhole)
% p.   model.probe_central_stop = true;                            % Model central stop
% p.   model.probe_diameter = 170e-6;                               % Model probe pupil diameter
% p.   model.probe_central_stop_diameter = 50e-6;                   % Model central stop diameter
% p.   model.probe_zone_plate_diameter = 170e-6;                    % Model probe zone plate diameter
% p.   model.probe_outer_zone_width = [];                           % Model probe zone plate outermost zone width (not used if not a focused probe) 
% p.   model.probe_propagation_dist = 1.2e-3;                       % Model probe propagation distance (pinhole <-> sample for unfocused, focal-plane <-> sample for focused)
% p.   model.probe_focal_length = 51e-3;                            % Model probe focal length (used only if model_is_focused is true
%                                                                   %   AND model_outer_zone_width is empty)
% p.   model.probe_upsample = 10;                                   % Model probe upsample factor (for focused probes)
% p.   model.probe_structured_illum_power = 1;                      % Intensity of phase variation in the lens 


p.   initial_probe_file = 'utils/imgs/probe_PSI.mat';     % Use probe from this mat-file (not used if model_probe is true)
p.   probe_file_propagation = [];                                % Distance for propagating the probe from file in meters, = [] to ignore





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
