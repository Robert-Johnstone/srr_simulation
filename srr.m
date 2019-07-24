%% Script to simulate MRI acquisition at super-resolution reconstruction.
%  The object is represented in 2D, with a set of 1D slices acquired
%  through it.
%  Doesn't currently simulate MR acquisition in-slice (in x-direction).

clear
close all

% Phantom parameters
phantom_radius = 100; % mm

% Acquisition parameters
fov = 300; % mm, governs number of slices, too
slice_thickness = 10; % mm
slice_spacing = 2; % mm
acq_resn = 2; % mm, in-slice resolution
slice_profile = 'gaussian';

% Simulation parameters
sim_resn = 0.1; % mm
disp_resn = 0.25; % mm
kernel_width = 40; % mm, full width of kernel, not slice width

% Display options
interp = 'nearest'; % Can be a cell array representing a blurring kernel

% Derived parameters
n_y_pts = (fov/sim_resn)+1;
n_x_pts = (fov/sim_resn)+1;
y = linspace(-fov/2,+fov/2,n_y_pts);
x = linspace(-fov/2,+fov/2,n_x_pts);

% Define blurring kernel
n_kernel_pts = (kernel_width/sim_resn)+1;
kernel_pts = linspace(-kernel_width/2,+kernel_width/2,n_kernel_pts);
switch slice_profile
    case 'gaussian'
        sigma = slice_thickness/(2*sqrt(2*log(2)));
        kernel = exp(-((kernel_pts/sigma).^2)/2)/(sigma*sqrt(2*pi));
end

% Define phantom
[X,Y] = meshgrid(x,y);
phantom=(X.^2+Y.^2)<phantom_radius^2;

% Iterate through the slices, exciting a slice and acquiring
slices = (fov/slice_spacing)+1;
slices_y = linspace(-fov/2,+fov/2,slices);
n_y_pts = (fov/sim_resn)+1;
y = linspace(-fov/2,+fov/2,n_y_pts);
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
    pause(.01)
    % Acquire echo
    echo = fft(sum(phant_excited,2));
    % Truncate echo to number of acquired samples - fov/acq_resn+1
    echo_truncated = [echo(1:fov/(2*acq_resn)+1); echo(n_x_pts-fov/(2*acq_resn)+1:end)];
    % Reconstruct slice
    image_slice = abs(ifft(echo_truncated));
%     plot(image_slice)
    % Store slice in 2D image
    img(:,slice) = image_slice;
end

% Display acquired image
img_disp = imresize(img,[(acq_resn/disp_resn)*size(img,1),(slice_spacing/disp_resn)*size(img,2)],'nearest');
imshow(img_disp',[])
title('Acquired image -- nearest-neighbour interpolation', 'Interpreter', 'latex')
xlabel('$x$ -- in-slice', 'Interpreter', 'latex');
ylabel('$y$ -- through-slice', 'Interpreter', 'latex');
