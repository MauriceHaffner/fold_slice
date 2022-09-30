% resolution = fourier_shell_corr_3D_2e(image,param)
% Computes the Fourier shell correlation for a single image by splitting
% into even-even, odd-odd, even-odd and odd-even parity sites. Images can be complex-valued.
% Can handle non-cube arrays but assumes the voxel is isotropic
% Modified by YJ for electron ptychography
%
% Inputs:
%     **fftkernel         fast fourier transform of convolution kernel
%     **image             image
%     **param             Structure containing parameters
% *optional*: 
%     **SNRt = 0.5        Power SNR for threshold, popular options:
%                         SNRt = 0.5;      1 bit threshold for average
%                         SNRt = 0.2071;   1/2 bit threshold for average
%     **thickring         Normally the pixels get assigned to the closest integer pixel ring in Fourier domain. 
%                         With thickring the thickness of the rings is increased by
%                         thickring, so each ring gets more pixels and more statistics
%     **auto_thickring    do not calculate overlaps if thickring > 1 is used
%     **freq_thr =0.05    mimimal freq value above which the resolution is detected   
% 
% returns: 
%     ++resolution        [min, max] resolution estimated from FSC curve

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%|                                                                       |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language this copyright should be retained and the authors 
%   and institution should be acknowledged in written form. Additionally 
%   you should cite the publication most relevant for the implementation 
%   of this code, namely
%   Vila-Comamala et al. "Characterization of high-resolution diffractive 
%   X-ray optics by ptychographic coherent diffractive imaging," Opt. 
%   Express 19, 21333-21344 (2011).
%   
%   Note however that the most relevant citation for the theoretical
%   foundation of the FSC criteria we use here is 
%   M. van Heela, and M. Schatzb, "Fourier shell correlation threshold
%   criteria," Journal of Structural Biology 151, 250-262 (2005).
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resolution= FSC_3d_e(fftkernel,image,param)
import math.isint
import utils.*
import engines.debluring.*

%%%%%%%%%%%%%%%%%%%%% PROCESS PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    param = struct();
end

parser = inputParser;
parser.addParameter('SNRt', 0.5 , @isnumeric )                       % SNRt = 0.2071 for 1/2 bit threshold for average of 2 images
                                                                        % SNRt = 0.5 for 1 bit threshold for average of 2 images
parser.addParameter('thickring',  0 , @isnumeric )                      % thick ring in Fourier domain
parser.addParameter('auto_binning',  false , @islogical )               % bin FRC before calculating rings, it makes calculations faster 
parser.addParameter('max_rings',  200 , @isnumeric )                    % maximal number of rings if autobinning is used 
parser.addParameter('freq_thr',  0 , @isnumeric )                    % mimimal freq value where resolution is detected  
parser.addParameter('pixel_size',  []  )                                % size of pixel in angstrom
parser.addParameter('mask',  [], @(x)(isnumeric(x) || islogical(x)) )   % array, equal to 0 for ignored pixels of the fft space and 1 for rest 

r = parser.Results;

% load all to the param structure 
for name = fieldnames(r)'
    if ~isfield(param, name{1})  % prefer values in param structure 
        param.(name{1}) = r.(name{1});
    end
end

if isempty(param.pixel_size)
    warning('Pixel size not specified. Please use param.pixel_size. \n');
    param.pixel_size = nan; 
end

image_real = real(image);
image_imag = imag(image);
deconv_image_imag = deconvlucy(image_imag, fftkernel,param.lucy_iters,param.SNRt);
deconv_image_total = image_real + 1j * deconv_image_imag;
sub_images = image_split(deconv_image_total);
split_masks = image_split(param.mask);

FSC = 0;
for idx = 1 : 1 : 2

    img1 = sub_images{idx};
    img2 = sub_images{idx+1};

    mask1 = split_masks{idx};
    mask2 = split_masks{idx+1};

    if any(size(img1) ~= size(img2))
        error('Images must be the same size')
    end

    % remove masked values from consideration (i.e. for laminography)
    F1 = fftn(img1);
    F2 = fftn(img2);   
    if ~isempty(param.mask)
        F1 = bsxfun(@times,F1 , mask1+eps);
        F2 = bsxfun(@times,F2 , mask2+eps);
    end
    F1cF2  = F1 .* conj(F2);
    F1 = abs(F1).^2; 
    F2 = abs(F2).^2;
    
    [ny,nx,nz] = size(img1);
    nmin = min(size(img1));

    thickring = param.thickring; 

    if param.auto_binning 
        % bin the correlation values to speed up the following calculations 
        % find optimal binning to make the volumes roughly cubic 
        bin = ceil(thickring/4) * floor(size(img1)/ nmin);        
        % avoid too large number of rings
        bin = max(bin, floor(nmin ./ param.max_rings));
        if any(bin > 1)
            thickring = ceil(thickring / min(bin));
            % fftshift and crop the arrays to make their size dividable by binning number 
            if ismatrix(img1); bin(3) = 1; end
            % force the binning to be centered 
            subgrid = {fftshift(ceil(bin(1)/2):(floor(ny/bin(1))*bin(1)-floor(bin(1)/2)-1)), ...
                   fftshift(ceil(bin(2)/2):(floor(nx/bin(2))*bin(2)-floor(bin(2)/2)-1)), ...
                   fftshift(ceil(bin(3)/2):(floor(nz/bin(3))*bin(3)-floor(bin(3)/2)-1))}; 
            if ismatrix(img1); subgrid(3) = [] ; end 
            % binning makes the shell / ring calculations much faster
            F1 = ifftshift(utils.binning_3D(F1(subgrid{:}), bin));
            F2 = ifftshift(utils.binning_3D(F2(subgrid{:}), bin));
            F1cF2 = ifftshift(utils.binning_3D(F1cF2(subgrid{:}), bin));
        end
    else
        bin = 1;
    end

    [ny,nx,nz] = size(F1);
    nmax = max([nx ny nz]);
    nmin = min(size(img1));

    % empirically tested that thickring should be >=3 along the smallest axis to avoid FRC undesampling 
    thickring = max(thickring, ceil(nmax/nmin)); 
    %param.thickring  = thickring;

    rnyquist = floor(nmax/2);
    freq = [0:rnyquist];

    x = ifftshift([-fix(nx/2):ceil(nx/2)-1])*floor(nmax/2)/floor(nx/2);
    y = ifftshift([-fix(ny/2):ceil(ny/2)-1])*floor(nmax/2)/floor(ny/2);
    if nz ~= 1
        z = ifftshift([-fix(nz/2):ceil(nz/2)-1])*floor(nmax/2)/floor(nz/2);
    else
        z = 0;
    end

    % deal with asymmetric pixel size in case of 2D FRC
    if  length(param.pixel_size) == 2
        if param.pixel_size(1) > param.pixel_size(2)
            y = y .* param.pixel_size(2) / param.pixel_size(1); 
        else
            x = x .* param.pixel_size(1) / param.pixel_size(2); 
        end
        param.pixel_size = min(param.pixel_size);  % FSC will be now calculated up to the maximal radius given by the smallest pixel size
    end
                                                               
    [X,Y,Z] = meshgrid(single(x),single(y),single(z));
    index = (sqrt(X.^2+Y.^2+Z.^2));

    clear X Y Z

    Nr = length(freq);    
    for ii = 1:Nr
        r = freq(ii);
        % calculate always thickring, min ring thickness is given by the smallest axis
        ind = index>=r-thickring/2 & index<=r+thickring/2;
        ind = find(ind);  % find seems to be faster then indexing 
        auxF1 = F1(ind); 
        auxF2 = F2(ind); 
        auxF1cF2 = F1cF2(ind);       
        C(ii)  = sum(auxF1cF2);
        C1(ii) = sum(auxF1);
        C2(ii) = sum(auxF2);
        n(ii) = numel(ind);  % Number of points
    end
    FSC = FSC + abs(C)./(sqrt(C1.*C2));
    
end
FSC = FSC / 2;
  
n = n*prod(bin);  % account for larger number of elements in the binned voxels

T = (  param.SNRt + 2*sqrt(param.SNRt)./sqrt(n+eps) + 1./sqrt(n)  )./...
    (  param.SNRt + 2*sqrt(param.SNRt)./sqrt(n+eps) + 1  );

freq_fine = 0:1e-3:max(freq);
freq_fine_normal = freq_fine/max(freq);

FSC_fine = max(0,interpn(freq, FSC, freq_fine, 'spline')); % spline, linear
T_fine = interpn(freq, T, freq_fine, 'spline');

idx_intersect = abs(FSC_fine-T_fine)< param.correlation_threshold;
%intersect_array = FSC_fine(idx_intersect);
range = freq_fine_normal(idx_intersect);
if length(range)<1
    range = [0 1];
    %intersect_array = [1 1];
end

%%%%%% CALCULATE STATISTICS %%%%%%%%%%%%%%
pixel = param.pixel_size; % angstrom

%range_start = range(find(range>param.freq_thr, 1, 'first'));
%if isempty(range_start)
%    range_start = range(1);
%end

resolution = pixel/range(1); %[pixel/range_start, pixel/range(end)];

end