%% Script to simulate MRI acquisition at super-resolution reconstruction.
%  The object is represented in 2D, with a set of 1D slices acquired
%  through it.
%  Doesn't currently simulate MR acquisition in-slice (in x-direction).

clear
close all

% Phantom parameters
phantom_radius = 100; % mm

% Acquisition parameters
fov = 300; % mm - must be even multiple of slice spacing
slice_thickness = 6; % mm
slice_spacing = 2; % mm - must divide fov to give even number
acq_resn = 2; % mm, in-slice resolution
slice_profile = 'gaussian'; % gaussian, rect, rect_adv, sinc

% Simulation parameters
sim_resn = 0.2; % mm

% SRR parameters
% Project kernel width in y pixels (units of slice spacing)
kernel_width = sqrt(slice_thickness^2-slice_spacing^2)/slice_spacing; % The 'right' width
% kernel_width = slice_thickness/slice_spacing; % The 'wrong' width

% Display options
interp = 'cubic'; % Can be a cell array representing a blurring kernel
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
phantom = double((X.^2+Y.^2)<phantom_radius^2);
% Add rectangular insert
rect_centre = [-30 40];
rect_height = 40; % mm
rect_width = 40; % mm
rect_angle = 45; % mm
phantom(and(abs((X*cosd(rect_angle)-Y*sind(rect_angle))-rect_centre(1))<(rect_width/2), ...
    abs((Y*cosd(rect_angle)+X*sind(rect_angle))-rect_centre(2))<(rect_height/2))) = 0;
% Add ellipsoid
elip_centre = [40 -40];
elip_height = 25; % mm
elip_width = 60; % mm
phantom((((X-elip_centre(1)).^2)./(elip_width/2)^2) + ...
    (((Y-elip_centre(2)).^2)./(elip_height/2)^2) < 1) = 0.2;
% Add circle
circ_centre = [-40 -20];
circ_diam = 40; % mm
phantom(((X-circ_centre(1)).^2) + ...
    ((Y-circ_centre(2)).^2) < ((circ_diam/2)^2)) = 0.4;
show_image(phantom,[500 500],'cubic','Phantom',0)

% Acquire MR image
img = mri_acq(phantom,fov,sim_resn,acq_resn,slice_thickness,slices,slice_profile,y);

% Create ground truth based on a slice thickness that corresponds to the
% slice spacing
ground_truth = mri_acq(phantom,fov,sim_resn,acq_resn,slice_spacing,slices,slice_profile,y);

% Perform SRR in through-slice (y) direction
srr_img = zeros(size(img));
for column_x = 1:acq_x_pts
    srr_img(column_x,:) = srrecon(img(column_x,:),'gaussian',kernel_width,ground_truth(column_x,:));
end

% Display images
disp_size = [(acq_resn/disp_resn)*size(img,1),(slice_spacing/disp_resn)*size(img,2)];
show_image(img,disp_size,interp,'Acquired image',0)
show_image(ground_truth,disp_size,interp,'Ground truth image',0)
show_image(srr_img,disp_size,interp,'SRR image (magnitude)',0)
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
