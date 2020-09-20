function plotSpectrum(img,image_width,label,xlimits,ylimits,dB,linecol,markers)
%PLOTSPECTRUM Plots a frequency spectrum of a signal (an image)
%   
%   img - the sampled signal - assumed uniform sampling
%   image_width - the physical width of the image (in metres)
%   label - the 'DisplayName property of the line
%   xlimits - the limits on the x-axis (in 1/metres)
%   dB - 0: linear y-axis, 1: log y-axis
    
    % Find the next power of 2 in size, for computational efficiency
    n = 2^nextpow2(size(img,2));
    % Fourier transform
    Y = fft(img,n);
    % Calculate the corresponding spatial frequencies
    f = ((size(img,2)-1)/image_width)*(0:(n/2))/n; % Stops at Nyquist limit
    %f = ((size(img,2)-1)/image_width)*(0:n-1)/n; % Full spectrum
    % Normalise and take magnitude
    P = abs(Y/n);
    % Plot
    if dB
        plot(f,10*log(P(1:n/2+1)/max(P)),'DisplayName',label,'Color',linecol,'LineStyle',markers) % Stops at Nyquist limit  
        ylabel(['$\textrm{fft}(signal)$ (dB)'], 'Interpreter', 'latex')
    else
        plot(f,P(1:n/2+1)/max(P),'DisplayName',label,'Color',linecol,'LineStyle',markers) % Stops at Nyquist limit  
        %plot(f,P/max(P),'DisplayName',label) % Full spectrum
        ylabel(['$\textrm{fft}(signal)$'], 'Interpreter', 'latex')
    end
    hold on 
    xlabel('$f\textrm{(m}^{-1})$', 'Interpreter', 'latex')
    xlim(xlimits)
    ylim(ylimits)
end

