%% TEST ALL ENGINES IN PTYCHOSHELVES TOOLKIT
% run this script to automatically test 
% all common features of the ptychoshelves toolkit 
% - GPU code will be tested only if a Matlab supported GPU is detected 
% - Fast c_solver will be tested only on unix systems 


close all
clear all

mkdir ./tomo
addpath ./tomo

%% %%%%%%%%% SET BASIC PATHS %%%%%%%%%%%%%%%%
base_package_path = '/users/stud/haffnerm/Ptychoshelves/fold_slice';  % use automatic search for the cSAXS_matoab_base package 
ptychoshelves_path = fileparts(mfilename('fullpath'));   % path to this ptychoshelves package 
% temporal_data_path = './temp/';  % directory to store temporal data created during tests 
base_path ='./' ;   % working path 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%% SET PARAMS %%%%%%%%%%%%%%%%
plot_results = true; 
check_simulated_data = true;  % only tests quality using simulated datasets, no automatic detector of quality is implemented
check_real_data = true;       % compares the achievable resolution with previous reconstructions; can be very time consuming
verbose_level = -1;           % default = -1 -> no verbosity, only in case of crash 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% add the relevant paths 
addpath(base_package_path)
addpath(ptychoshelves_path)
addpath(fullfile(ptychoshelves_path, 'tests'))
import utils.*



%% run different test templates 
if check_simulated_data
disp( 'TESTING ENGINES ON CPU')
run('CPU_engines_test.m')
disp( 'TESTING MULTILAYER ENGINES ON CPU')
run('CPU_multilayer_test.m')
disp( 'TESTING ENGINES WITH BINNING')
run('CPU_binning_test.m')
disp( 'TESTING CPU-BASED LOW COHERENCE PTYCHO ')
run('CPU_engines_test_low_coherence.m')


if gpuDeviceCount > 0
    disp('TESTING GPU-BASED ENGINES')
    run('GPU_engines_test.m')
    disp( 'TESTING GPU-BASED MULTILAYER EXTENSION')
    run('GPU_multilayer_test.m')
    disp( 'TESTING GPU-BASED UPSAMPLED NEARFIELD PTYCHO ')
    run('GPU_nearfield_test.m')
    disp( 'TESTING GPU-BASED LOW COHERENCE PTYCHO ')
    run('GPU_engines_test_low_coherence.m')
    disp( 'TESTING OFFAXIS PTYCHO ')
    run('GPU_offaxis_ptycho_test.m')
    disp( 'TESTING GPU-BASED PTYCHO WITH SCAN MIRRORING, SHARING')
    run('GPU_engines_test_mirror_scans.m')
    disp( 'TESTING GPU-BASED PTYCHO WITH SHIFTED SCAN SHARING')
    run('GPU_engines_test_shifted_scan_sharing.m')
    disp( 'TESTING GPU-BASED PTYCHO GEOMETRY REFINEMENT WITH SCAN MIRRORING, SHARING')
    run('GPU_engines_refine_geometry.m')
    disp( 'TESTING GPU-BASED PTYCHO CAMERA POSITION REFINEMENT WITH SCAN MIRRORING, SHARING')
    run('GPU_engines_camera_geometry.m')
end
end

[~, host] = system('hostname');
if strcmpi(host(1:2), 'ra') && check_real_data
    disp( 'TESTING REAL DATASET (NeXus EIGER1p5M H5 - OMNY) ')
    run('NeXusEiger1p5M_omny_test.m');
    disp( 'TESTING REAL DATASET (EIGER1p5M H5 - OMNY) ')
    run('Eiger1p5M_omny_test.m');
%    disp( 'TESTING REAL DATASET (PILATUS CBF - OMNY) ')
%    run('Pilatus_omny_test.m');
elseif check_real_data
    warning('Tests of the real datasets is implemented only for the ra-cluster in PSI')
end

disp( 'TESTING PLOTTING ')
run('test_plotting.m')

if exist('+ptychotomo/', 'dir') && strcmpi(input('Do you want to test iterative multilayer-ptychotomography? [y/N]','s'), 'y')
    disp('TESTING PTYCHOTOMOGRAPHY, IT CAN TAKE A WHILE ... ')
    run('GPU_multilayer_3D_test.m')
else
    disp('Skipping')
end

    
disp( 'TESTING FINISHED')





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
