
% Variables
function [sph_n, com ,car_scene_v, bbox] = radar_car_model(model, range_translation, angle_threshold)
   
    [ sph_n, com ,car_scene_v, bbox] = car_model(model, range_translation, angle_threshold );
end

function [sph_n,com, car_scene_v, bbox] = car_model(model,range_translation , angle_threshold )
  
    path_to_pointCloud = 'CAD/CAD_model_%d.mat';
    CAD_idx = 1; % CAD_ID of the model car (1 -38)
    % 2 does not work
    % 5 works
    % 6 does not work
    % 7 does not work
    % 8 does not work
    disp_pc = 1; % Display images
    k = 6; % nearest neighbours to consider for change in normal
    car_v = loadCarPointCloud(path_to_pointCloud, CAD_idx , disp_pc);

    [car_scene_v,fig_1]  = rotateAndTranslate(car_v,range_translation, disp_pc);

    %Find the center of mass of the point cloud (Mean)
    com = mean(car_scene_v.cart_v);

    [visible_cart_v] = remove_occlusion(car_scene_v); % remove occluded body of the car
    % Theshold to detect strong variations in the surface normal...
    % Increase this to detect sharper edges
    threshold = 2;
    [strong_reflectors, normal_reflectors ,normals] =...
        model_reflectors(visible_cart_v, threshold, k, disp_pc, angle_threshold);

     % Hawkeye model
     blb_cart_v = model_point_reflector(visible_cart_v,car_scene_v.bbox,fig_1, disp_pc);

    if model == "edge"
        [sph_n(:,1),sph_n(:,2),sph_n(:,3)] = cart2sph(strong_reflectors(:,1),strong_reflectors(:,2),strong_reflectors(:,3));
    elseif model == "normal"
        [sph_n(:,1),sph_n(:,2),sph_n(:,3)] = cart2sph(normal_reflectors(:,1),normal_reflectors(:,2),normal_reflectors(:,3));
    else
        [sph_n(:,1),sph_n(:,2),sph_n(:,3)] = cart2sph(blb_cart_v(:,1),blb_cart_v(:,2),blb_cart_v(:,3));
    end
       
    bbox = car_scene_v.bbox;
end
%%

function car_v= loadCarPointCloud(path_to_pointCloud, CAD_idx , disp_pc)
    % CAD models are loaded as point clouds of size N_pt by 3, where N_pt
    % is the number of points and 3 values are the cartesian coordinates
    % unit is mm
    load(sprintf(path_to_pointCloud,CAD_idx));
    disp_pc = 0;
    if(disp_pc)
        % Visulize the original point cloud
        figure; 
        cart_v_plot = cart_v;
        cart_v_plot = datasample(cart_v, 1000); % downsampling when plotting
        scatter3(cart_v_plot(:,1),cart_v_plot(:,2),cart_v_plot(:,3),10,'filled','k'); hold on;
        xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); axis equal;
        set(gca,'FontSize',10) % Creates an axes and sets its FontSize to 18
        view([-12,14]);
        title('Original point cloud');
        clear cart_v_plot;
    end
    variable_library;
    % store point cloud in pc (point cloud) structure
    car_v = car_v_struct;
    car_v.CAD_idx = CAD_idx;
    car_v.N_pt = length(cart_v);
    car_v.cart_v = cart_v;
    car_v.lim = [min(cart_v);max(cart_v)]; % find the limits in all three dimensions 
    % 8 vertices of the bounding box of the point cloud
    [bbox_x, bbox_y, bbox_z] = meshgrid(car_v.lim(:,1),car_v.lim(:,2),car_v.lim(:,3)); 
    car_v.bbox = [bbox_x(:), bbox_y(:), bbox_z(:)]; 
    clear cart_v bbox N_pt car_idx;
end

function [car_scene_v,fig_1] = rotateAndTranslate(car_v,range_translation, disp_pc)
    car_scene_v = car_v;
    variable_library;
    %% Rotate     
    car_scene_v.rotate = rotate_ang(randi(length(rotate_ang))); % randomly select a rotation angle and store it in the pc structure
    car_scene_v.rotate = mod(car_scene_v.rotate,180);
    % inline function for 2D rotation
    rotate2d =  @(x, M) (x(:, 1:2) * M);
    rotate_angle_rad = car_scene_v.rotate/180*pi;
     % create rotation matrix
    rotation_matrix = [cos(rotate_angle_rad), -sin(rotate_angle_rad); sin(rotate_angle_rad), cos(rotate_angle_rad)];
    % rotate the point cloud 
    car_scene_v.cart_v(:,1:2) = rotate2d(car_scene_v.cart_v, rotation_matrix); 
    % rotate the bounding box
    car_scene_v.bbox(:,1:2) = rotate2d(car_scene_v.bbox, rotation_matrix); 
    % update the limits in all three dimensions
    car_scene_v.lim = [min(car_scene_v.cart_v);max(car_scene_v.cart_v)]; 

    %% Translation
    % range of translation along x and y axis
    %translate_x_rng = (translate_lim(1,1) - car_scene_v.lim(1,1)):translate_x_res:(translate_lim(1,2) - car_scene_v.lim(2,1)); 
    %translate_y_rng = (translate_lim(2,1) - car_scene_v.lim(1,2)):translate_y_res:(translate_lim(2,2) - car_scene_v.lim(2,2));
    % randomly select a translation distances
    translate_x = floor(1000 + rand(1)*range_translation*1000); % in mm
    translate_y = floor((rand(1))*range_translation*1000); % in mm
    translate_z = 0;

    % translate
    car_scene_v.translate = [translate_x, translate_y, translate_z]; % store translation information in the pc structure
    car_scene_v.cart_v = car_scene_v.cart_v + car_scene_v.translate; % translate the point cloud
    car_scene_v.bbox = car_scene_v.bbox + car_scene_v.translate; % translate the bounding box
    car_scene_v.lim = [min(car_scene_v.cart_v);max(car_scene_v.cart_v)]; % update the limits in all three dimensions
    
    % convert unit from mm to m
    car_scene_v.cart_v = car_scene_v.cart_v/1000; 
    car_scene_v.bbox = car_scene_v.bbox/1000; 
    
    if(disp_pc)
        % Visulize the rotated and translated point cloud
        fig_1 = figure(); 
        cart_v_plot = car_scene_v.cart_v; % downsampling when 
        scatter3(cart_v_plot(:,1),cart_v_plot(:,2),cart_v_plot(:,3),10,'Color',[128, 128, 128]); hold on;
        scatter3(car_scene_v.bbox(:,1), car_scene_v.bbox(:,2),car_scene_v.bbox(:,3),'red');
        hold on;
        xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal;
      
        xlim([min(cart_v_plot(:,1))-4 max(cart_v_plot(:,1))+4]);
        ylim([min(cart_v_plot(:,2))-4 max(cart_v_plot(:,2))+2]);
        zlim([0 max(cart_v_plot(:,3))+2]);
        set(gca,'FontSize',10) % Creates axes
        view([-12,14]);
        title('Rotated and translated point cloud');
        hold on;
    end
end

function [strong_reflectors,normal_reflectors, normals] = ...
    model_reflectors(visible_PC, threshold,k, disp_pc, angle_threshold)
    % This function calculates the change in the surface normal:
    % Normals which point directly towards the radar are strong reflectors
    % we assume points close to large changes in the surface normal also
    % strongly reflect towards the Radar
    
    
    %Point Cloud data structure
    ptCloud = pointCloud(visible_PC);
    disp_pc = 0;
    % Normal Esitmation 
    normals = pcnormals(ptCloud);
    ptCloud.Normal = normals; 
    % Flip normal towards the centre
    ptCloud = flip_Normal(ptCloud);
    % Look for strong variations in the change in normal 
    % Search for k nearest neighbors for each normal
    nearest_neighbors = zeros(ptCloud.Count, k);
    for i = 1 : (ptCloud.Count)
        [indices,dists] = findNearestNeighbors(ptCloud,ptCloud.Normal(i,:),k);
        nearest_neighbors(i,:) = indices';
    end
    strong_reflectors = [];
    normal_reflectors = [];
    normals =[];
    u = ptCloud.Normal(1:end,1);
    v = ptCloud.Normal(1:end,2);
    w = ptCloud.Normal(1:end,3);
    for j = 1 : numel(u)
        relfected_angle = acosd(dot (ptCloud.Location(j,:), ptCloud.Normal(j,:))...
            /norm(ptCloud.Location(j,:))/norm(ptCloud.Normal(j,:)));
        % normal_reflectors = [normal_reflectors; relfected_angle];
        normals = [normals ; relfected_angle];
        if relfected_angle > angle_threshold
            normal_reflectors = [normal_reflectors; ptCloud.Location(j,:)];
        end
        vec_diff = [(u(j) - u(nearest_neighbors(j,:))), ...
            v(j) - v(nearest_neighbors(j,:)), ...
            w(j) - w(nearest_neighbors(j,:))];
        vec_diff(isnan(vec_diff))=0;
        if norm(vec_diff) > threshold
            strong_reflectors = [strong_reflectors; ptCloud.Location(j,:)];
        end
    end

    if(disp_pc)
        step_size = 5; 
        x = ptCloud.Location(1:step_size:end,1);
        y = ptCloud.Location(1:step_size:end,2);
        z = ptCloud.Location(1:step_size:end,3);
        u = ptCloud.Normal(1:step_size:end,1);
        v = ptCloud.Normal(1:step_size:end,2);
        w = ptCloud.Normal(1:step_size:end,3);
        figure;
        pcshow(ptCloud);
        title('Estimated Normals of Point Cloud')
        hold on;
        quiver3(x,y,z,u,v,w);
        hold off;

        figure;
        scatter3(normal_reflectors(:,1),normal_reflectors(:,2),normal_reflectors(:,3),20,'filled','o');
        hold on;
        scatter3(strong_reflectors(:,1),strong_reflectors(:,2),strong_reflectors(:,3),20,'filled','k');
        title('Strong reflectors of Point Cloud');
        legend('Normal reflector', 'Edge reflector');
        xlim([min(normal_reflectors(:,1))-4 max(normal_reflectors(:,1))+4]);
        ylim([0 max(normal_reflectors(:,1))+2]);
        
        xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal;
        hold off;
    end
end

function ptCloud =  flip_Normal(ptCloud)
% Flip the normals to point towards the sensor location
    sensorCenter = [0,0,0];
    x = ptCloud.Location(1:end,1);
    y = ptCloud.Location(1:end,2);
    z = ptCloud.Location(1:end,3);
    u = ptCloud.Normal(1:end,1);
    v = ptCloud.Normal(1:end,2);
    w = ptCloud.Normal(1:end,3);
    for k = 1 : numel(x)
        p1 = sensorCenter - [x(k),y(k),z(k)];
        p2 = [u(k),v(k),w(k)];
        % Flip the normal vector if it is not pointing towards the sensor.
        angle = atan2(norm(cross(p1,p2)),p1*p2');
        if angle > pi/2 || angle < -pi/2
            u(k) = -u(k);
            v(k) = -v(k);
            w(k) = -w(k);
        end
    end
    ptCloud.Normal(1:end,1) = u;
    ptCloud.Normal(1:end,2) = v;
    ptCloud.Normal(1:end,3) = w;
end





