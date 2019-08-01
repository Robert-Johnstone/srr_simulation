function show_image(img,disp_size,interp,img_title,scale_bar)
%SHOW_IMAGE Displays a 2D image

    % Switch for showing labels
    show_labels = 1;
    
    img_disp = imresize(img,disp_size,interp);
    figure
    imshow(rot90(img_disp),[]);
    if show_labels
        title(img_title, 'Interpreter', 'latex')
        xlabel('In-slice', 'Interpreter', 'latex');
        ylabel('Through-slice', 'Interpreter', 'latex');
    end
    if scale_bar
        colorbar
    end
end

