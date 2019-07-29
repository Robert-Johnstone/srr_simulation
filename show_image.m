function show_image(img,disp_size,interp,img_title,scale_bar)
%SHOW_IMAGE Displays a 2D image

    img_disp = imresize(img,disp_size,interp);
    figure
    imshow(rot90(img_disp),[]);
    title(img_title, 'Interpreter', 'latex')
    xlabel('$x$ -- in-slice', 'Interpreter', 'latex');
    ylabel('$y$ -- through-slice', 'Interpreter', 'latex');
    if scale_bar
        colorbar
    end
end

