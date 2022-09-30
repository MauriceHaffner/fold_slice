function object_deconv = apply_deconvolution(p,object)

    import engines.debluring.*
    import math.*
    import utils.*
    import io.*

    objectROI = object(p.object_ROI{:},:);
    objectROI = make_sizes_even(objectROI);
    output_cell = optimize_kernel(p,object);
    p.optimal_kernel = output_cell{1};
    p.optimal_kernel_params = output_cell{2};
    objectROI_real = real(objectROI);
    objectROI_imag = imag(objectROI);
    objectROI_imag_deconv = deconvlucy(objectROI_imag, p.optimal_kernel,p.deconvlucy_iters,p.SNRt);
    object_deconv = objectROI_real + 1j * objectROI_imag_deconv;

end    

