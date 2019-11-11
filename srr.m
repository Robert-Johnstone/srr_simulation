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
fp_kernel_type = 'gaussian'; % guassian, <filename>, generated
bp_kernel_type = 'same'; % guassian, <filename>, generated, same [as FP kernel]
kernel_width = sqrt(slice_thickness^2-slice_spacing^2)/slice_spacing; % The 'right' width
% kernel_width = slice_thickness/slice_spacing; % The 'wrong' width
% Note - kernel_width not used for 'generated' FP kernel

% Derived parameters
sim_y_pts = (fov/sim_resn)+1; % Number of simulation points in y-direction
sim_x_pts = (fov/sim_resn)+1; % Number of simulation points in x-direction
y = linspace(-fov/2,+fov/2,sim_y_pts); % Simulated y points
x = linspace(-fov/2,+fov/2,sim_x_pts); % Simulated x points
acq_y_pts = (fov/acq_resn)+1; % Number of acquired points in y-direction
acq_x_pts = (fov/acq_resn)+1; % Number of acquired points in x-direction
y_acq = linspace(-fov/2,+fov/2,acq_y_pts); % Acquired y points
x_acq = linspace(-fov/2,+fov/2,acq_x_pts); % Acquired x points
slices = (fov/slice_spacing)+1; % Number of slices

% Display options
display_images = 0; % Whether to display images
interp = 'cubic'; % Can be a cell array representing a blurring kernel
disp_resn = 0.5; % mm
disp_size = [(acq_resn/disp_resn)*(fov/acq_resn+1),(slice_spacing*slices/disp_resn)];
save_images = 1;
bw = 1; % Black and white plots

% Generate phantom
phantom = make_phantom(phantom_radius,fov,sim_resn);
if display_images
    show_image(phantom,disp_size,'cubic','Phantom',0)
end
if save_images
    save_image(phantom,disp_size,'cubic','phantom.png')
end

% Acquire noisy LR MR image
lr_img = mri_acq(phantom,fov,sim_resn,acq_resn,slice_thickness,slices,slice_profile,y,acq_snr);

% Acquire noisy HR MR image with reduced SNR
hr_img = mri_acq(phantom,fov,sim_resn,acq_resn,slice_spacing,slices,...
    slice_profile,y,acq_snr*slice_spacing/slice_thickness);

% Create ground truth based on a slice thickness that corresponds to the
% slice spacing
ground_truth = mri_acq(phantom,fov,sim_resn,acq_resn,slice_spacing,slices,slice_profile,y,inf);

% Perform SRR in through-slice (y) direction
srr_img = zeros(size(lr_img));
fprintf('Performing SR recon: Column ');
cstr = ''; % Counter string
for column_x = 1:acq_x_pts
    % Update counter
    fprintf(repmat('\b',1,length(cstr))); % Perform carriage return
    cstr = [num2str(column_x) ' of ' num2str(acq_x_pts)];
    fprintf(cstr);
    % Do SRR
    srr_img(column_x,:) = srrecon(lr_img(column_x,:),fp_kernel_type,kernel_width,bp_kernel_type,ground_truth(column_x,:));
end
fprintf('\n');

if display_images
    % Display images fov/acq_resn+1,slices
    show_image(lr_img,disp_size,interp,'Acquired LR image',0)
    show_image(hr_img,disp_size,interp,'Acquired HR image',0)
    show_image(ground_truth,disp_size,interp,'Ground truth image',0)
    show_image(srr_img,disp_size,interp,'SRR image (magnitude)',0)
    show_image((lr_img-ground_truth),disp_size,interp,'Absolute error image for acquired image',1)
    show_image((srr_img-ground_truth),disp_size,interp,'Absolute error image for SRR',1)
end

if save_images
    % Save results as images
    fn_root = [num2str(slice_thickness) 'mm_at_' num2str(slice_spacing) 'mm_'];
    fn_root = [fn_root fp_kernel_type '_'];
    fn_root = regexprep(fn_root,'.mat',''); % Remove .mat from filename
    save_image(lr_img,disp_size,interp,[fn_root 'mri_acq_lr.png'])
    save_image(hr_img,disp_size,interp,[fn_root 'mri_acq_hr.png'])
    save_image(ground_truth,disp_size,interp,[fn_root 'mri_gt.png'])
    save_image(srr_img,disp_size,interp, [fn_root 'srr.png'])
    save_image(0.5+(lr_img-ground_truth),disp_size,interp, [fn_root 'acq_error.png'])
    save_image(0.5+(srr_img-ground_truth),disp_size,interp,[fn_root 'srr_error.png'])
end

% Compare central lines profiles
fig = figure;
if bw
    plot(x_acq,lr_img(ceil(acq_x_pts/2),:),'k--')
    hold on
    plot(x_acq,srr_img(ceil(acq_x_pts/2),:),'k')
    plot(x_acq,ground_truth(ceil(acq_x_pts/2),:),'k:')
else
    plot(x_acq,lr_img(ceil(acq_x_pts/2),:),'b')
    hold on
    plot(x_acq,srr_img(ceil(acq_x_pts/2),:),'r')
    plot(x_acq,ground_truth(ceil(acq_x_pts/2),:),'g')
end
% title('Comparison of central line profiles', 'Interpreter', 'latex')
xlabel('y / mm','Interpreter','latex')
ax = gca;
ax.YLim = [-0.05 1.05];
fig.Position(3:4) = [560 210];
if save_images
    saveas(gcf,[fn_root 'profiles'],'epsc')
end
% Plot Fourier transform
fig = figure;
if bw
    plotSpectrum(lr_img(ceil(acq_x_pts/2),:),fov/1000,'FT',...
        [0 500/slice_spacing],[-100 0],1,'k','--')
    hold on
    plotSpectrum(srr_img(ceil(acq_x_pts/2),:),fov/1000,'FT',...
        [0 500/slice_spacing],[-100 0],1,'k','-')
    plotSpectrum(ground_truth(ceil(acq_x_pts/2),:),fov/1000,'FT',...
        [0 500/slice_spacing],[-100 0],1,'k',':')
else
    plotSpectrum(lr_img(ceil(acq_x_pts/2),:),fov/1000,'FT',...
        [0 500/slice_spacing],[-100 0],1,'b','-')
    hold on
    plotSpectrum(srr_img(ceil(acq_x_pts/2),:),fov/1000,'FT',...
        [0 500/slice_spacing],[-100 0],1,'r','-')
    plotSpectrum(ground_truth(ceil(acq_x_pts/2),:),fov/1000,'FT',...
        [0 500/slice_spacing],[-100 0],1,'g','-')
end
fig.Position(3:4) = [560 210];
if save_images
    saveas(gcf,[fn_root 'profiles_ft'],'epsc')
end

% Calculate errors
l1error_srr = norm(srr_img(:)-ground_truth(:),1);
fprintf('L1 error of SRR image: %d\n', l1error_srr)
l1error_lr = norm(lr_img(:)-ground_truth(:),1);
fprintf('L1 error of LR image: %d\n', l1error_lr)
l1error_hr = norm(hr_img(:)-ground_truth(:),1);
fprintf('L1 error of HR image: %d\n', l1error_hr)
