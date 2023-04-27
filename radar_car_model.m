
% Variables
function [sph_n, com ,normals] = radar_car_model(model, range_translation, angle_threshold)
   
    [ sph_n, com ,normals] = car_model(model, range_translation, angle_threshold );
end

function [sph_n,com, normals] = car_model(model,range_translation , angle_threshold )
  
    path_to_pointCloud = 'CAD/CAD_model_%d.mat';
    CAD_idx = 5; % CAD_ID of the model car (1 -38)
    disp_pc = 0; % Display images
    k = 6; % nearest neighbours to consider for change in normal
    car_v = loadCarPointCloud(path_to_pointCloud, CAD_idx , disp_pc);

    car_scene_v = rotateAndTranslate(car_v,range_translation, disp_pc);

    %Find the center of mass of the point cloud (Mean)
    com = mean(car_scene_v.cart_v);

    [visible_cart_v] = remove_occlusion(car_scene_v); % remove occluded body of the car
    % Theshold to detect strong variations in the surface normal...
    % Increase this to detect sharper edges
    threshold = 2;
    [strong_reflectors, normal_reflectors ,normals] =...
        model_reflectors(visible_cart_v, threshold, k, disp_pc, angle_threshold);

     % Hawkeye model
     blb_cart_v = model_point_reflector(visible_cart_v,car_scene_v.bbox);

    if model == "edge"
        [sph_n(:,1),sph_n(:,2),sph_n(:,3)] = cart2sph(strong_reflectors(:,1),strong_reflectors(:,2),strong_reflectors(:,3));
    elseif model == "normal"
        [sph_n(:,1),sph_n(:,2),sph_n(:,3)] = cart2sph(normal_reflectors(:,1),normal_reflectors(:,2),normal_reflectors(:,3));
    else
        [sph_n(:,1),sph_n(:,2),sph_n(:,3)] = cart2sph(blb_cart_v(:,1),blb_cart_v(:,2),blb_cart_v(:,3));
    end
       

end
%%

function car_v= loadCarPointCloud(path_to_pointCloud, CAD_idx , disp_pc)
    % CAD models are loaded as point clouds of size N_pt by 3, where N_pt
    % is the number of points and 3 values are the cartesian coordinates
    % unit is mm
    load(sprintf(path_to_pointCloud,CAD_idx));
    if(disp_pc)
        % Visulize the original point cloud
        figure; 
        cart_v_plot = cart_v;
        cart_v_plot = datasample(cart_v, 1000); % downsampling when plotting
        scatter3(cart_v_plot(:,1),cart_v_plot(:,2),cart_v_plot(:,3),10,'filled','k'); hold on;
        xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); axis equal;
        set(gca,'FontSize',30) % Creates an axes and sets its FontSize to 18
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

function car_scene_v = rotateAndTranslate(car_v,range_translation, disp_pc)
    car_scene_v = car_v;
    variable_library;
    %% Rotate     
    car_scene_v.rotate = rotate_ang(randi(length(rotate_ang))); % randomly select a rotation angle and store it in the pc structure
    car_scene_v.rotate = mod(car_scene_v.rotate*(randi(1)*2-1),180);
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
    translate_x = floor(-5000 + rand(1)*range_translation*1000); % in mm
    translate_y = floor(rand(1)*range_translation*3000); % in mm
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
        figure; 
        cart_v_plot = car_scene_v.cart_v; % downsampling when plotting
        scatter3(cart_v_plot(:,1),cart_v_plot(:,2),cart_v_plot(:,3),10,'filled','k'); hold on;
        scatter3(car_scene_v.bbox(:,1), car_scene_v.bbox(:,2),car_scene_v.bbox(:,3),'r'); hold on;
        xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal;
        xlim([0 12]);
        ylim([0 12]);
        set(gca,'FontSize',30) % Creates axes
        title('Rotated and translated point cloud');
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
        xlim([0 12]);
        ylim([0 12]);
        set(gca,'FontSize',30) % Creates axes
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hawkeye code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function blb_cart_v = model_point_reflector(visible_cart_v,bbox)
% Model radar point reflectors in the scene
% Input: 'visible_cart_v': point cloud of the visible body of the car (cartesian coordinates) 
%        'bbox': bounding box of the car
% Output: 'blb_cart_v': blobs of radar point reflectos

    % local function find the distance from a point pt to a line determined by 2 points v1 and v2
    function d = point_to_line(pt, v1, v2)
          a = [v1 - v2,0];
          b = [pt - v2,0];
          d = norm(cross(a,b)) / norm(a);
    end
    
    %% Blobs center
    blb_ctr_cart = []; % 
    N_pt_car = length(visible_cart_v);
    pt_specular = zeros(N_pt_car,1);

%     box_close = sqrt((bbox(2:4,1)-bbox(1,1)).^2+(bbox(2:4,2)-bbox(1,2)).^2);
%     [~,box_close] = sort(box_close);
%     bbox(2:4,1:2) = bbox(box_close+1,1:2);

    bbox_edge = zeros(4,1); % 4 edges of the bounding box represented as a complex number
    bbox_edge_cart = zeros(4,2,2); % coordinates of the vertices of the bounding box edges 

    bbox_edge_cart(1,:,:) = squeeze(bbox([1,2],1:2));
    bbox_edge_cart(2,:,:) = squeeze(bbox([1,3],1:2));
    bbox_edge_cart(3,:,:) = squeeze(bbox([2,4],1:2));
    bbox_edge_cart(4,:,:) = squeeze(bbox([3,4],1:2));
    for ke = 1:4
        bbox_edge(ke,:) = bbox_edge_cart(ke,2,1)-bbox_edge_cart(ke,1,1)+1j*(bbox_edge_cart(ke,2,2)-bbox_edge_cart(ke,1,2));
    end

    % find the specularity for each point
    for kp = 1:N_pt_car

        if min(sqrt((visible_cart_v(kp,1)-bbox(1:4,1)).^2+(visible_cart_v(kp,2)-bbox(1:4,2)).^2))<0.5
            % if the point is very close to a corner of the bounding box
            pt_specular(kp) = 0;
        else
            % find the distance between the point and all 4 edges to
            % find the nearest edge
            pt_edge_dist = zeros(4,1);
            for ke = 1:4
                pt_edge_dist(ke) = point_to_line([visible_cart_v(kp,1),visible_cart_v(kp,2)], squeeze(bbox_edge_cart(ke,1,:)).', squeeze(bbox_edge_cart(ke,2,:)).');
            end
            [~,edge_nearest] = min(pt_edge_dist);
            % find the angle between the insertion angle and angle of the nearest edge  
            pt_specular(kp) = abs(angle(bbox_edge(edge_nearest)) - angle(visible_cart_v(kp,1)+1j*visible_cart_v(kp,2)));
            pt_specular(kp) = abs(mod(pt_specular(kp),pi)/pi*180-90);
            % find the elevation angle theta
            pt_theta = angle(sqrt(visible_cart_v(kp,1)^2+visible_cart_v(kp,2)^2)+1j*abs(visible_cart_v(kp,3)));
            pt_theta = abs(mod(pt_theta,pi)/pi*180);
            % the larger one is the specularity of the point
            pt_specular(kp) = max(pt_specular(kp),pt_theta);
        end
    end

    % Blob center
    % Find the centers of the blobs
    blb_ctr = find(pt_specular==0); % all corners are selected
    blb_ctr_cart = [];
    if ~isempty(blb_ctr)
        blb_ctr_idx = datasample(blb_ctr,min(10,length(blb_ctr)));
        blb_ctr_cart = [blb_ctr_cart;visible_cart_v(blb_ctr_idx,:)];
    end

    blb_ctr = find((pt_specular>0)&(pt_specular<15));
    if ~isempty(blb_ctr)
        blb_ctr_idx = datasample(blb_ctr,min(20,length(blb_ctr)));
        blb_ctr_cart = [blb_ctr_cart;visible_cart_v(blb_ctr_idx,:)];
    end

    blb_ctr = find((pt_specular>15)&(pt_specular<25));
    if ~isempty(blb_ctr)
        blb_ctr_idx = datasample(blb_ctr,min(5,length(blb_ctr)));
        blb_ctr_cart = [blb_ctr_cart;visible_cart_v(blb_ctr_idx,:)];
    end

    %%
    if ~isempty(blb_ctr_cart)

%         figure;
%         scatter3(blb_ctr_cart(:,1),blb_ctr_cart(:,2),blb_ctr_cart(:,3)); hold on;
%         title('Center of blobs');

        blb_size = 0.3; % blob size around the center
        blb_cart_v = [];
        for kb = 1:size(blb_ctr_cart,1)
            dis_pt2blb = (visible_cart_v(:,1) - blb_ctr_cart(kb,1)).^2 + (visible_cart_v(:,2) - blb_ctr_cart(kb,2)).^2 + (visible_cart_v(:,3) - blb_ctr_cart(kb,3)).^2;
            ptInBlb = find(dis_pt2blb < blb_size^3);
            blb_cart_v = [blb_cart_v;[visible_cart_v(ptInBlb,1),visible_cart_v(ptInBlb,2),visible_cart_v(ptInBlb,3)]];
        end
        blb_cart_v = unique(blb_cart_v,'row');
    end
end


