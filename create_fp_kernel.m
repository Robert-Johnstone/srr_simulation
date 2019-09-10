function fp_kernel = create_fp_kernel(fp_kernel_size)
%CREATE_FP_KERNEL Creates a forward projection kernel
%
%   Given a low-resolution slice profile and a high-resolution slice
%   profile, the function attempts to find the foward projection kernel
%   that, when convolved with the HR profile give the LR profile.
%
%   fp_kernel_size - size of kernel array (integer)

    % Set number of iterations for iterative back projection
    max_iter = 100;
    
    % Switch to turn plotting errors on and off
    plot_errors = 0;
    
    res_ratio = 3; % Ratio of HR to LR
    
    if plot_errors
        figure
        plot(0,norm(fp_kernel-ground_truth),'x')
        hold on
    end

    % Create forward projection kernel
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

    for i = 1:max_iter
        % Forward project using HR slice profile
        slice_profile_lr_guess = conv(fp_kernel,slice_profile_hr,'same');
        % Find error
        slice_profile_lr_error = slice_profile_lr_guess-slice_profile_lr;
        % Back project error using HR slice profile
        slice_profile_hr_error = conv(slice_profile_lr_error,slice_profile_hr,'same');
        fp_kernel = fp_kernel-slice_profile_hr_error;
    end
    % Normalise
    fp_kernel = fp_kernel/sum(fp_kernel);
end

