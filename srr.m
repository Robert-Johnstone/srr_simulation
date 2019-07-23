%% Script to simulate MRI acquisition at super-resolution reconstruction.
%  The object is represented in 2D, with a set of 1D slices acquired
%  through it.
%  Doesn't currently simulate MR acquisition in-slice (in x-direction).

clear
close all

% Phantom parameters
phantom_radius = 100; % mm

% Acquisition parameters
fov = 400; % mm, governs number of slices, too
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

% % Define sampling matrix
% % 1st dimension: in-slice co-ordinate
% % 2nd dimension: through slice co-ordinate
% % 3rd dimension: sample number
% n_samples = (fov/acq_resn) * (fov/slice_thickness);
% sampling = zeros(n_samples);

% Define blurring kernel
n_kernel_pts = (kernel_width/sim_resn)+1;
kernel_pts = linspace(-kernel_width/2,+kernel_width/2,n_kernel_pts);
switch slice_profile
    case 'gaussian'
        sigma = slice_thickness/(2*sqrt(2*log(2)));
        kernel = exp(-((kernel_pts/sigma).^2)/2)/(sigma*sqrt(2*pi));
        
end

% Define phantom vector
phant = zeros((fov/acq_resn)+1,(fov/slice_thickness)+1);
phant = phantom('Modified Shepp-Logan',200);

% Problem is separable, so iterate through pixels in slice (along x)
slices = (fov/slice_spacing)+1;
n_y_pts = (fov/sim_resn)+1;
n_x_pts = (fov/acq_resn)+1;
y = linspace(-fov/2,+fov/2,n_y_pts);
x = linspace(-fov/2,+fov/2,n_x_pts);
sampling_y = mod(y,slice_spacing)<sim_resn/2;
img = zeros(n_x_pts,(fov/slice_spacing)+1);
for slice_x = -fov/2:acq_resn:fov/2
    % Create phantom
    phant_y = (slice_x^2 + y.^2)<(phantom_radius^2);
    % Blur
    blurred_phant_y = conv(phant_y,kernel,'same');
    % Sample
    img_y = blurred_phant_y(sampling_y);
    % Store in image
    img(x==slice_x,:) = img_y;
end

% Simulate MR acquisition in x-direction
img_hybrid=fft(img,fov/sim_resn,1);
img_mri=ifft(img_hybrid,n_x_pts,1);

% Display acquired image
img_mri_disp = imresize(img,[(acq_resn/disp_resn)*size(img,1),(slice_spacing/disp_resn)*size(img,2)],'nearest');
imshow(img_disp')
title('Acquired image - nearest-neighbour interpolation', 'Interpreter', 'latex')
xlabel('$x - in-slice$', 'Interpreter', 'latex');
ylabel('$y - through-slice$', 'Interpreter', 'latex');
