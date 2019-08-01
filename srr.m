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
slice_profile = 'sg_150_100_167.mat'; % gaussian, rect, rect_adv, sinc, <filename>
acq_snr = inf; % Signal to noise ratio for acquisition

% Simulation parameters
sim_resn = 0.2; % mm

% SRR parameters
% Project kernel width in y pixels (units of slice spacing)
fp_kernel_type = 'generated';
bp_kernel_type = 'same';
kernel_width = sqrt(slice_thickness^2-slice_spacing^2)/slice_spacing; % The 'right' width
% kernel_width = slice_thickness/slice_spacing; % The 'wrong' width

% Derived parameters
sim_y_pts = (fov/sim_resn)+1; % Number of simulation points in y-direction
sim_x_pts = (fov/sim_resn)+1; % Number of simulation points in x-direction
y = linspace(-fov/2,+fov/2,sim_y_pts); % Simulated y points
x = linspace(-fov/2,+fov/2,sim_x_pts); % Simulated x points
acq_y_pts = (fov/acq_resn)+1; % Number of acquired points in y-direction
acq_x_pts = (fov/acq_resn)+1; % Number of acquired points in x-direction
slices = (fov/slice_spacing)+1; % Number of slices

% Display options
interp = 'cubic'; % Can be a cell array representing a blurring kernel
disp_resn = 0.5; % mm
disp_size = [(acq_resn/disp_resn)*(fov/acq_resn+1),(slice_spacing*slices/disp_resn)];
save_images = 1;

% Generate phantom
phantom = make_phantom(phantom_radius,fov,sim_resn);
show_image(phantom,disp_size,'cubic','Phantom',0)
if save_images
    imwrite(rot90(phantom), 'phantom.png')
end

% Acquire MR image
img = mri_acq(phantom,fov,sim_resn,acq_resn,slice_thickness,slices,slice_profile,y,acq_snr);

% Create ground truth based on a slice thickness that corresponds to the
% slice spacing
ground_truth = mri_acq(phantom,fov,sim_resn,acq_resn,slice_spacing,slices,slice_profile,y,inf);

% Perform SRR in through-slice (y) direction
srr_img = zeros(size(img));
fprintf('Performing SR recon: Column ');
cstr = ''; % Counter string
for column_x = 1:acq_x_pts
    % Update counter
    fprintf(repmat('\b',1,length(cstr))); % Perform carriage return
    cstr = [num2str(column_x) ' of ' num2str(acq_x_pts)];
    fprintf(cstr);
    % Do SRR
    srr_img(column_x,:) = srrecon(img(column_x,:),fp_kernel_type,kernel_width,bp_kernel_type,ground_truth(column_x,:));
end
fprintf('\n');

% Display images fov/acq_resn+1,slices
show_image(img,disp_size,interp,'Acquired image',0)
show_image(ground_truth,disp_size,interp,'Ground truth image',0)
show_image(srr_img,disp_size,interp,'SRR image (magnitude)',0)
show_image((img-ground_truth),disp_size,interp,'Absolute error image for acquired image',1)
show_image((srr_img-ground_truth),disp_size,interp,'Absolute error image for SRR',1)

% Save results as images
if save_images
    fn_root = [num2str(slice_thickness) 'mm_thick_' num2str(slice_spacing) 'mm_spacing_'];
    imwrite(rot90(img), [fn_root 'mri_acq.png'])
    imwrite(rot90(ground_truth), [fn_root 'mri_gt.png'])
    imwrite(rot90(srr_img), [fn_root 'srr.png'])
    imwrite(0.5+(rot90(img-ground_truth))*1, [fn_root 'acq_error.png'])
    imwrite(0.5+(rot90(srr_img-ground_truth))*1, [fn_root 'srr_error.png'])
end

% Compare central lines profiles
figure
plot(img(ceil(acq_x_pts/2),:))
hold on
plot(srr_img(ceil(acq_x_pts/2),:))
plot(ground_truth(ceil(acq_x_pts/2),:))
title('Comparison of central line profiles', 'Interpreter', 'latex')
xlabel('y','Interpreter','latex')

% Calculate errors
l1error_srr = sum(sum(abs(srr_img-ground_truth)));
fprintf('L1-error of SRR image: %d\n', l1error_srr)
l1error_lr = sum(sum(abs(img-ground_truth)));
fprintf('L1-error of LR image: %d\n', l1error_lr)
