function [hr_image] = srrecon(lr_image,kernel_type,kernel_width,ground_truth)
%SRRECON Performs a super-resolution reconstruction
%   Takes a 1D image vector and, using a specified kernel_type, performed a
%   super-resolution reconstruction
%
%   input_image - 1D image
%   kernel_type - kernel type, e.g. 'gaussian'
%   kernel_width - the nominal width of the kernel, e.g. FWHM, in pixels
%                   (units of slice spacing)

    % Set number of iterations for iterative back projection
    max_iter = 1000;
    
    % Switch to turn plotting errors on and off
    plot_errors = 0;

    % Create forward and backward projection kernels
    m = 3; % Multiplier used to set how big vector representing kernel is
    n_kernel_pts = ceil(kernel_width)*m+1;
    % Make sure kernel and image either both have an odd number or both
    % have an even number of pixels. This avoids conv function introducing 
    % a shift
    if mod(n_kernel_pts,2)~=mod(size(lr_image,2),2)
        n_kernel_pts = n_kernel_pts+1;
    end
    fp_kernel = zeros(n_kernel_pts,1);
    kernel_pts = linspace(-(n_kernel_pts-1)/2,(n_kernel_pts-1)/2,n_kernel_pts);
    switch kernel_type
        case 'gaussian'
            sigma = kernel_width/(2*sqrt(2*log(2)));
            fp_kernel = exp(-((kernel_pts/sigma).^2)/2)/(sigma*sqrt(2*pi));
    end
    bp_kernel = fp_kernel;
    
    % Create initial high resultion image guess
    hr_image = lr_image;

    if plot_errors
        figure
        plot(0,norm(hr_image-ground_truth),'x')
        hold on
    end
    last_error = Inf;

    for i = 1:max_iter
        % Forward project
        lr_image_guess = conv(hr_image,fp_kernel,'same');
        % Find error
        lr_image_error = lr_image_guess-lr_image;
        % Back project error
        output_image_error = conv(lr_image_error,bp_kernel,'same');
        hr_image_temp = hr_image-output_image_error;
        % Calculate error metric and bail out if it is getting bigger
        current_error = norm(hr_image_temp-ground_truth);
        if current_error > last_error
            break
        else
            last_error = current_error;
            hr_image = hr_image_temp;
        end
        if plot_errors
            plot(i,norm(hr_image-ground_truth),'x')
        end
    end
end

