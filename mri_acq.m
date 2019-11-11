function [img] = mri_acq(phantom,fov,sim_resn,acq_resn,slice_thickness,slices,slice_profile,y,snr)
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
            spw = 6; % Conventional slice width for saved slice profile, mm
            spr = 0.001; % Conventional resolution for saved slice profile, mm
    end

    sim_x_pts = (fov/sim_resn)+1;

    fprintf('Performing MRI scan: Slice ');
    cstr = ''; % Counter string

    % Iterate through the slices, exciting a slice and acquiring
    slices_y = linspace(-fov/2,+fov/2,slices);
    img = zeros(fov/acq_resn+1,slices);
    for slice = 1:slices
        % Update counter
        fprintf(repmat('\b',1,length(cstr))); % Perform carriage return
        cstr = [num2str(slice) ' of ' num2str(slices)];
        fprintf(cstr);
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
                kernel_shifted = kernel_shifted/sum(kernel_shifted);
            otherwise
                kernel_shifted = interp1((-0.24:1e-6:0.24)*slice_thickness/spw,profile,(y-slice_pos)*spr,'linear',0);
                % Normalise - incorrect calculation for slice partially in
                % volume
                kernel_shifted = kernel_shifted/sum(kernel_shifted);
        end
        excitation = repmat(kernel_shifted,(fov/sim_resn)+1,1);
        phant_excited = phantom.*excitation;
        % Acquire echo
        echo = fft(sum(phant_excited,2))/sqrt(size(phant_excited,1));
        % Truncate echo to number of acquired samples - fov/acq_resn+1
        echo_truncated = [echo(1:ceil(fov/(2*acq_resn))+1); echo(sim_x_pts-floor(fov/(2*acq_resn))+1:end)];
        % Scale to match the truncation
        echo_truncated = echo_truncated * sqrt((fov/acq_resn+1) / size(phant_excited,1));
        % Add some noise with amplitude such that background noise in image
        % has standard deviation equal to 1/SNR
        echo_truncated = echo_truncated + ...
            (randn(size(echo_truncated))+1i*randn(size(echo_truncated))) ...
            * (2-pi/2)^(-.5)/snr;
        % Reconstruct slice
        image_slice = abs(ifft(echo_truncated))*sqrt(size(echo_truncated,1));
        % Store slice in 2D image
        img(:,slice) = image_slice;
    end
    fprintf('\n');
end

