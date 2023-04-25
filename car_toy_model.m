


function sph_v = car_toy_model(range_translation , car_size)
% Toy model of car
  % Model the car as 4 corners :
    %  0,car_size --------- car_size,car_size
    %  ---                  ---
    %  ---       CAR        ---
    %  ---                  ---
    %  0,0 ---------------- car_size,0
% Scale and rotate the car
    % 2 point model
%     N_target = 2;
%     input = [0, 0 ,0 ; 0, car_size, 0];

    N_target = 4;
    input = [0, 0 ,0 ; 0, car_size, 0 ; car_size, 0, 0 ; car_size, car_size,0];
   
    theta = -180+360*rand([1,3]);
    % Assuming car is 5 times the car length's distance away from radar
    % Translate 5*rand() distance in x, y, z with variances 1,1,0
    translation = range_translation.*([rand(),rand(),0]);
    tform = rigidtform3d(theta,translation);
    ptCloud = pointCloud(input);
    new_inputs = pctransform(ptCloud, tform);
    new_inputs = remove_occlusion(new_inputs);
    fig_1 = figure();
    plot3(new_inputs(:,1), new_inputs(:,2), new_inputs(:,3),'o')
    %axis([-5 5 -5 5 -5 5]);
    title('Rotated and translated point cloud');
    
    %Convert coordinates to amplitude vector and angle vectors

    sph_v = zeros([N_target, 3]);
    [sph_v(:,1),sph_v(:,2),sph_v(:,3)] = cart2sph(new_inputs(:,1),new_inputs(:,2),new_inputs(:,3));

end

function [visible_cart_v] = remove_occlusion(car_scene_v)
% Remove occluded body of the car
% Input: 'car_scene_v' point cloud of the car in the scene
% Output: 'visible_cart_v' point cloud of the visible body car to the radar (catesian coordinates)
    
    % convert cartesian coordinates of the point cloud to spherical coordinates
    sph_v = zeros(car_scene_v.Count,3); 
    [sph_v(:,1),sph_v(:,2),sph_v(:,3)] = cart2sph(car_scene_v.Location(:,1),car_scene_v.Location(:,2),car_scene_v.Location(:,3));
    sphlim = [min(sph_v);max(sph_v)]; % limits in all three dimensions in the spherical coordinates            
  
    % define voxel size of the grid
    phi_res = 0.2/180*pi;
    theta_res = 0.5/180*pi;
    rho_res = 0.02;
    
    sph_m_phi = sphlim(1,1):phi_res:sphlim(2,1)+phi_res;
    sph_m_theta = sphlim(1,2):theta_res:sphlim(2,2)+theta_res;
    sph_m_rho = sphlim(1,3):rho_res:sphlim(2,3)+rho_res;
    sph_m_size = [length(sph_m_phi),length(sph_m_theta),length(sph_m_rho)];
    sph_m = zeros(sph_m_size);

    phi_m_idx = round((sph_v(:,1) - sphlim(1,1))/phi_res)+1;
    theta_m_idx = round((sph_v(:,2) - sphlim(1,2))/theta_res)+1;
    rho_m_idx = round((sph_v(:,3) - sphlim(1,3))/rho_res)+1;
    for k_pt = 1:car_scene_v.Count
        sph_m(phi_m_idx(k_pt),theta_m_idx(k_pt),rho_m_idx(k_pt)) = 1;
    end

    %% Find first voxel with car body in every angle
    visible_sph_m = zeros(size(sph_m));
    for kphi = 1:sph_m_size(1)
        for ktheta = 1:sph_m_size(2)
            krho = find(sph_m(kphi,ktheta,:)>0,1);
            visible_sph_m(kphi,ktheta,krho) = sph_m(kphi,ktheta,krho);
        end
    end

    visible_sph_m_idx = find(visible_sph_m);
    sph_v_idx = [];
    [sph_v_idx(:,1),sph_v_idx(:,2),sph_v_idx(:,3)] = ind2sub(sph_m_size,visible_sph_m_idx);
    visible_sph_v = [sph_m_phi(sph_v_idx(:,1));sph_m_theta(sph_v_idx(:,2));sph_m_rho(sph_v_idx(:,3))].';

    visible_cart_v = zeros(size(visible_sph_v));
    [visible_cart_v(:,1),visible_cart_v(:,2),visible_cart_v(:,3)] = sph2cart(visible_sph_v(:,1),visible_sph_v(:,2),visible_sph_v(:,3));
end