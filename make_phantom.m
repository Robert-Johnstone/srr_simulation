function phantom = make_phantom(phantom_radius,fov,sim_resn)
%MAKE_PHANTOM Makes a phantom

    sim_y_pts = (fov/sim_resn)+1;
    sim_x_pts = (fov/sim_resn)+1;
    y = linspace(-fov/2,+fov/2,sim_y_pts);
    x = linspace(-fov/2,+fov/2,sim_x_pts);
    % Define phantom
    [X,Y] = meshgrid(x,y);
    phantom = double((X.^2+Y.^2)<phantom_radius^2);
    % Add rectangular insert
    rect_centre = [-30 40]*(phantom_radius/100);
    rect_height = 40*(phantom_radius/100);
    rect_width = 40*(phantom_radius/100);
    rect_angle = 45*(phantom_radius/100);
    phantom(and(abs((X*cosd(rect_angle)-Y*sind(rect_angle))-rect_centre(1))<(rect_width/2), ...
        abs((Y*cosd(rect_angle)+X*sind(rect_angle))-rect_centre(2))<(rect_height/2))) = 0;
    % Add ellipsoid
    elip_centre = [40 -40]*(phantom_radius/100);
    elip_height = 25*(phantom_radius/100);
    elip_width = 60*(phantom_radius/100);
    phantom((((X-elip_centre(1)).^2)./(elip_width/2)^2) + ...
        (((Y-elip_centre(2)).^2)./(elip_height/2)^2) < 1) = 0.2;
    % Add circle
    circ_centre = [-40 -20]*(phantom_radius/100);
    circ_diam = 40*(phantom_radius/100);
    phantom(((X-circ_centre(1)).^2) + ...
        ((Y-circ_centre(2)).^2) < ((circ_diam/2)^2)) = 0.4;
end

