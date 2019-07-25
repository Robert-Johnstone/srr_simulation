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
slice_spacing = 2; % mm - must give even number of pixels in slice
acq_resn = 2; % mm, in-slice resolution
slice_profile = 'gaussian';

% Simulation parameters
sim_resn = 0.1; % mm

% Display options
interp = 'nearest'; % Can be a cell array representing a blurring kernel
disp_resn = 0.5; % mm

% Derived parameters
sim_y_pts = (fov/sim_resn)+1;
sim_x_pts = (fov/sim_resn)+1;
y = linspace(-fov/2,+fov/2,sim_y_pts);
x = linspace(-fov/2,+fov/2,sim_x_pts);
acq_y_pts = (fov/acq_resn)+1;
acq_x_pts = (fov/acq_resn)+1;
slices = (fov/slice_spacing)+1;

% Define phantom
[X,Y] = meshgrid(x,y);
phantom=double((X.^2+Y.^2)<phantom_radius^2);

% Acquire MR image
img = mri_acq(phantom,fov,sim_resn,acq_resn,slice_thickness,slices,slice_profile,y);

% Create ground truth based on a slice thickness that corresponds to the
% slice spacing
ground_truth = mri_acq(phantom,fov,sim_resn,acq_resn,slice_spacing,slices,slice_profile,y);

% Perform SRR in through-slice (y) direction
srr_img = zeros(size(img));
kernel_width = sqrt(slice_thickness^2-acq_resn^2)/acq_resn; % In pixels
for column_x = 1:acq_x_pts
    srr_img(column_x,:) = srrecon(img(column_x,:),'gaussian',kernel_width,ground_truth(column_x,:));
end

% Display images
disp_size = [(acq_resn/disp_resn)*size(img,1),(slice_spacing/disp_resn)*size(img,2)];
show_image(img,disp_size,interp,'Acquired image',0)
show_image(ground_truth,disp_size,interp,'Ground truth image',0)
show_image(srr_img,disp_size,interp,'SRR image',0)
show_image((img-ground_truth),disp_size,interp,'Absolute error image for acquired image',1)
show_image((srr_img-ground_truth),disp_size,interp,'Absolute error image for SRR',1)

% Compare central lines profiles
figure
plot(img(ceil(acq_x_pts/2),:))
hold on
plot(srr_img(ceil(acq_x_pts/2),:))
plot(ground_truth(ceil(acq_x_pts/2),:))
title('Comparison of central line profiles', 'Interpreter', 'latex')
xlabel('y','Interpreter','latex')
