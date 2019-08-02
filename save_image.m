function save_image(img,disp_size,interp,fn)
%SAVE_IMAGE Saves a 2D image

    img_disp = imresize(img,disp_size,interp);
    imwrite(rot90(img_disp),fn);
end

