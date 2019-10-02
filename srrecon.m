function [hr_image] = srrecon(lr_image,fp_kernel_type,fp_width,bp_kernel_type,ground_truth)
%SRRECON Performs a super-resolution reconstruction
%   Takes a 1D image vector and, using a specified kernel_type, performed a
%   super-resolution reconstruction
%
%   input_image - 1D image
%   fp_kernel_type - FP kernel type, e.g. 'gaussian',<filename>
%   fp_width - the nominal width of the kernel, e.g. FWHM, in pixels
%                   (units of slice spacing)
%   bp_kernel - BP kernel type, e.g. 'gaussian',<filename>

    % Set number of iterations for iterative back projection
    max_iter = 10;
    
    % Switch to turn plotting errors on and off
    plot_errors = 0;
    
    % Switch for terminating optimisation if error increases
    diverge_stop = 0;

    % Create forward projection kernel
    m = 3; % Multiplier used to set how big vector representing kernel is
    n_kernel_pts = ceil(fp_width)*m+1;
    % Make sure kernel and image either both have an odd number or both
    % have an even number of pixels. This avoids conv function introducing 
    % a shift
    if mod(n_kernel_pts,2)~=mod(size(lr_image,2),2)
        n_kernel_pts = n_kernel_pts+1;
    end
    kernel_pts = linspace(-(n_kernel_pts-1)/2,(n_kernel_pts-1)/2,n_kernel_pts);
    switch fp_kernel_type
        case 'gaussian'
            sigma = fp_width/(2*sqrt(2*log(2)));
            fp_kernel = exp(-((kernel_pts/sigma).^2)/2)/(sigma*sqrt(2*pi));
        case 'generated'
            fp_kernel = create_fp_kernel(n_kernel_pts);
        otherwise
            % Load saved profile
            load(fp_kernel_type,'profile');
            spw = 6; % Conventional slice width for saved slice profile, mm
            spr = 0.001; % Conventional resolution for saved slice profile, mm
            fp_kernel = interp1((-0.24:1e-6:0.24)*fp_width/spw,profile,kernel_pts*spr,'linear',0);
            % Normalise
            fp_kernel = fp_kernel/sum(fp_kernel);
    end
    
    % Create backward projection kernel
    if strcmp(fp_kernel_type,bp_kernel_type) || strcmp(bp_kernel_type,'same')
        bp_kernel = fp_kernel;
    else
        if strcmp(bp_kernel_type,'gaussian')
            sigma = fp_width/(2*sqrt(2*log(2)));
            bp_kernel = exp(-((kernel_pts/sigma).^2)/2)/(sigma*sqrt(2*pi));
        end
    end
    
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
        if diverge_stop
            hr_image_temp = hr_image-output_image_error;
            % Calculate error metric and bail out if it is getting bigger
            current_error = norm(hr_image_temp-ground_truth);
            if current_error > last_error
                break
            else
                last_error = current_error;
                hr_image = hr_image_temp;
            end
        else
            hr_image = hr_image-output_image_error;
        end
        if plot_errors
            plot(i,norm(hr_image-ground_truth),'x')
        end
    end
end

