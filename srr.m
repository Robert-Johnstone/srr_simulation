%% Script to simulate MRI acquisition at super-resolution reconstruction.
%  The object is represented in 2D, with a set of 1D slices acquired
%  through it.
%  Doesn't currently simulate MR acquisition in-slice (in x-direction).

clear
close all

% Phantom parameters
phantom_radius = 100; % mm

% Acquisition parameters
fov = 304; % mm, governs number of slices, too
slice_thickness = 10; % mm
slice_spacing = 4; % mm - must give even number of pixels in slice
acq_resn = 4; % mm, in-slice resolution
slice_profile = 'gaussian';

% Simulation parameters
sim_resn = 0.1; % mm
disp_resn = 0.25; % mm

% Display options
interp = 'cubic'; % Can be a cell array representing a blurring kernel

% Derived parameters
sim_y_pts = (fov/sim_resn)+1;
sim_x_pts = (fov/sim_resn)+1;
y = linspace(-fov/2,+fov/2,sim_y_pts);
x = linspace(-fov/2,+fov/2,sim_x_pts);
acq_y_pts = (fov/acq_resn)+1;
acq_x_pts = (fov/acq_resn)+1;
slices = (fov/slice_spacing)+1;

% Set up parameters for kernel
switch slice_profile
    case 'gaussian'
        sigma = slice_thickness/(2*sqrt(2*log(2)));
end

% Define phantom
[X,Y] = meshgrid(x,y);
phantom=(X.^2+Y.^2)<phantom_radius^2;

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

% Display acquired image
img_disp = imresize(img,[(acq_resn/disp_resn)*size(img,1),(slice_spacing/disp_resn)*size(img,2)],interp);
imshow(img_disp',[])
title('Acquired image -- nearest-neighbour interpolation', 'Interpreter', 'latex')
xlabel('$x$ -- in-slice', 'Interpreter', 'latex');
ylabel('$y$ -- through-slice', 'Interpreter', 'latex');

% Perform SRR in through-slice (y) direction
srr_img = zeros(size(img));
kernel_width = sqrt(slice_thickness^2-acq_resn^2);
for column_x = 1:acq_x_pts
    srr_img(column_x,:) = srrecon(img(column_x,:),'gaussian',kernel_width);
end

% Display SR reconstructed image
img_disp = imresize(srr_img,[(acq_resn/disp_resn)*size(img,1),(slice_spacing/disp_resn)*size(img,2)],interp);
figure
imshow(img_disp',[])
title('SRR image -- nearest-neighbour interpolation', 'Interpreter', 'latex')
xlabel('$x$ -- in-slice', 'Interpreter', 'latex');
ylabel('$y$ -- through-slice', 'Interpreter', 'latex');

