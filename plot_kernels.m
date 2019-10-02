%   fp_kernel_type - FP kernel type, e.g. 'gaussian',<filename>
%   fp_width - the nominal width of the kernel, e.g. FWHM, in pixels
%                   (units of slice spacing)
slice_thickness = 6; % mm
slice_spacing = 2; % mm
fp_kernel_type = 'generated';
fp_width = sqrt(slice_thickness^2-slice_spacing^2)/slice_spacing; % The 'right' width

% Create forward projection kernel
m = 3; % Multiplier used to set how big vector representing kernel is
n_kernel_pts = ceil(fp_width)*m;
% Make sure kernel size is odd
if mod(n_kernel_pts,2)==0
    n_kernel_pts = n_kernel_pts + 1;
end

kernel_pts = linspace(-(n_kernel_pts-1)/2,(n_kernel_pts-1)/2,n_kernel_pts);
switch fp_kernel_type
    case 'gaussian'
        sigma = fp_width/(2*sqrt(2*log(2)));
        fp_kernel = exp(-((kernel_pts/sigma).^2)/2)/(sigma*sqrt(2*pi));
    case 'generated'
        fp_kernel = create_fp_kernel(n_kernel_pts);
    otherwise
        % Load saved profile
        load(fp_kernel_type,'profile');
        spw = 6; % Conventional slice width for saved slice profile, mm
        spr = 0.001; % Conventional resolution for saved slice profile, mm
        fp_kernel = interp1((-0.24:1e-6:0.24)*fp_width/spw,profile,kernel_pts*spr,'linear',0);
        % Normalise
        fp_kernel = fp_kernel/sum(fp_kernel);
end

% Compare central lines profiles
fig = figure;
plot(kernel_pts,fp_kernel,'k-')
xlabel('y / pixels','Interpreter','latex')
% ax = gca;
% ax.YLim = [-0.05 1.05];
% fig.Position(3:4) = [560 210];
% if save_images
%     saveas(gcf,[fn_root 'profiles'],'epsc')
% end
% Plot Fourier transform
fig = figure;
plotSpectrum(fp_kernel,n_kernel_pts*slice_spacing/1000,'FT',...
    [0 500/slice_spacing],[-100 0],1,'k','--')
% fig.Position(3:4) = [560 210];
% if save_images
%     saveas(gcf,[fn_root 'profiles_ft'],'epsc')
% end
