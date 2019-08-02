function fp_kernel = create_fp_kernel(fp_kernel_size)
%CREATE_FP_KERNEL Creates a forward projection kernel
%
%   Given a low-resolution slice profile and a high-resolution slice
%   profile, the function attempts to find the foward projection kernel
%   that, when convolved with the HR profile give the LR profile.
%
%   fp_kernel_size - size of kernel array (integer)

    % Set number of iterations for iterative back projection
    max_iter = 1000;
    
    % Switch to turn plotting errors on and off
    plot_errors = 0;
    
    res_ratio = 3; % Ratio of HR to LR
    
    if plot_errors
        figure
        plot(0,norm(fp_kernel-ground_truth),'x')
        hold on
    end
%     last_error = Inf;

    % Create forward projection kernel
%     m = 3; % Multiplier used to set how big vector representing kernel is
%     n_kernel_pts = ceil(fp_width)*m+1;
%     % Make sure kernel and image either both have an odd number or both
%     % have an even number of pixels. This avoids conv function introducing 
%     % a shift
%     if mod(n_kernel_pts,2)~=mod(size(lr_image,2),2)
%         n_kernel_pts = n_kernel_pts+1;
%     end
    kernel_pts = linspace(-(fp_kernel_size-1)/2,(fp_kernel_size-1)/2,fp_kernel_size);
    % Load saved slice profile and resample of kernel point spacing
    load('sg_150_100_167.mat','profile');
    spw = 6; % Conventional slice width for saved slice profile, mm
    spr = 0.001; % Conventional resolution for saved slice profile, mm
    slice_profile_lr = interp1((-0.24:1e-6:0.24)*res_ratio/spw,profile,kernel_pts*spr,'linear',0);
    slice_profile_hr = interp1((-0.24:1e-6:0.24)/spw,profile,kernel_pts*spr,'linear',0);
    % Normalise
    slice_profile_lr = slice_profile_lr/sum(slice_profile_lr);
    slice_profile_hr = slice_profile_hr/sum(slice_profile_hr);
    
    % Create initial high resultion image guess
    fp_kernel = slice_profile_lr;

%     figure
%     plot(slice_profile_lr,'+')
%     hold on
    
    for i = 1:max_iter
        % Forward project using HR slice profile
        slice_profile_lr_guess = conv(fp_kernel,slice_profile_hr,'same');
        % Find error
        slice_profile_lr_error = slice_profile_lr_guess-slice_profile_lr;
        % Back project error using HR slice profile
        slice_profile_hr_error = conv(slice_profile_lr_error,slice_profile_hr,'same');
        fp_kernel = fp_kernel-slice_profile_hr_error;
%         % Calculate error metric and bail out if it is getting bigger
%         current_error = norm(fp_kernel_temp-ground_truth);
%         if current_error > last_error
%             break
%         else
%             last_error = current_error;
%             fp_kernel = fp_kernel_temp;
%         end
%         if plot_errors
%             plot(i,norm(fp_kernel-ground_truth),'x')
%         end
%      plot(slice_profile_lr_guess)
    end
    % Normalise
    fp_kernel = fp_kernel/sum(fp_kernel);
end

