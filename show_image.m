function [outputArg1,outputArg2] = show_image(img,disp_size,interp,img_title)
%SHOW_IMAGE Summary of this function goes here
%   Detailed explanation goes here
    img_disp = imresize(img,disp_size,interp);
    figure
    imshow(img_disp',[])
    title(img_title, 'Interpreter', 'latex')
    xlabel('$x$ -- in-slice', 'Interpreter', 'latex');
    ylabel('$y$ -- through-slice', 'Interpreter', 'latex');
end

