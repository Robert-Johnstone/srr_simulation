function [img] = mri_acq(phantom,fov,sim_resn,acq_resn,slice_thickness,slices,slice_profile,y)
%MRI_ACQ Acquires a 2D MR image

    % Set up parameters for kernel (slice profile)
    switch slice_profile
        case 'gaussian'
            sigma = slice_thickness/(2*sqrt(2*log(2)));
        case 'rect'
        case 'rect_adv'
        case 'sinc'
        otherwise
            % Assume that a filename has been specifed and load slice
            % profile
            load(slice_profile,'profile');
    end

    sim_x_pts = (fov/sim_resn)+1;

    % Iterate through the slices, exciting a slice and acquiring
    slices_y = linspace(-fov/2,+fov/2,slices);
    img = zeros(fov/acq_resn+1,slices);
    for slice = 1:slices
        % Excite phantom
        slice_pos = slices_y(slice);
        % slice profile is normalised in real space (mm)
        switch slice_profile
            case 'gaussian'
                kernel_shifted = exp(-(((y-slice_pos)/sigma).^2)/2)/(sigma*sqrt(2*pi));
            case 'rect'
                kernel_shifted = and((y-slice_pos)<slice_thickness/2, ...
                    (y-slice_pos)>=-slice_thickness/2)/(slice_thickness/sim_resn);
            case 'rect_adv'
                % Set positions fully inside slice to 1
                kernel_shifted = (abs(y-slice_pos)<=(slice_thickness-sim_resn)/2);
                % Set positions partially in slice to weighted value
                kernel_shifted = kernel_shifted+(abs(abs(y-slice_pos)-slice_thickness/2)<sim_resn/2) ...
                    .* (abs(abs(y-slice_pos)-(slice_thickness-sim_resn)/2)/sim_resn);
                % Normalise
                kernel_shifted = kernel_shifted/(slice_thickness/sim_resn);
            case 'sinc' % With FWHM = slice_thickness and truncated at first zero crossing
                kernel_shifted = sinc(2*(y-slice_pos)*1.895/(pi*slice_thickness));
                kernel_shifted = kernel_shifted.*(abs(2*(y-slice_pos)*1.895)<(pi*slice_thickness));
                % Normalise - incorrect calculation for slice partially in
                % volume
                kernel_shifted = kernel_shifted/(sim_resn*sum(kernel_shifted));
            otherwise
                kernel_shifted = interp1((-0.24:1e-6:0.24)*slice_thickness/6,profile,(y-slice_pos)/1000,'linear',0);
                % Normalise - incorrect calculation for slice partially in
                % volume
                kernel_shifted = kernel_shifted/(sum(profile)*slice_thickness/(1000*6));
       end
        excitation = repmat(kernel_shifted,(fov/sim_resn)+1,1);
        phant_excited = phantom.*excitation;
    %     imshow(phant_excited',[0 .1]);
    %     pause(.01)
        % Acquire echo
        echo = fft(sum(phant_excited,2)*sim_resn);
        % Truncate echo to number of acquired samples - fov/acq_resn+1
        echo_truncated = [echo(1:ceil(fov/(2*acq_resn))+1); echo(sim_x_pts-floor(fov/(2*acq_resn))+1:end)];
        % Reconstruct slice
        image_slice = (sim_resn/acq_resn)*abs(ifft(echo_truncated));
        % Store slice in 2D image
        img(:,slice) = image_slice;
    end
end

