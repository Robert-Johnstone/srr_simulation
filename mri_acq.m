function [img] = mri_acq(phantom,fov,sim_resn,acq_resn,slice_thickness,slices,slice_profile,y)
%MRI_ACQ Acquires a 2D MR image
%   Detailed explanation goes here

    % Set up parameters for kernel
    switch slice_profile
        case 'gaussian'
            sigma = slice_thickness/(2*sqrt(2*log(2)));
    end

    sim_x_pts = (fov/sim_resn)+1;

    % Iterate through the slices, exciting a slice and acquiring
    slices_y = linspace(-fov/2,+fov/2,slices);
    img = zeros(fov/acq_resn+1,slices);
    for slice = 1:slices
        % Excite phantom
        slice_pos = slices_y(slice);
        switch slice_profile
            case 'gaussian'
                kernel_shifted = exp(-(((y-slice_pos)/sigma).^2)/2)/(sigma*sqrt(2*pi));
        end
        excitation = repmat(kernel_shifted,(fov/sim_resn)+1,1);
        phant_excited = phantom.*excitation;
    %     imshow(phant_excited',[0 .1]);
    %     pause(.01)
        % Acquire echo
        echo = fft(sum(phant_excited,2));
        % Truncate echo to number of acquired samples - fov/acq_resn+1
        echo_truncated = [echo(1:fov/(2*acq_resn)+1); echo(sim_x_pts-fov/(2*acq_resn)+1:end)];
        % Reconstruct slice
        image_slice = abs(ifft(echo_truncated));
        % Store slice in 2D image
        img(:,slice) = image_slice;
    end
end

