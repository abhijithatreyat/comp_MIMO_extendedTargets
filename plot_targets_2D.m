function plot_targets_2D(true_targets,visible_cart_v, omega_vect, N_tx,...
     N_chirp, detected_targets,detected_angles, g_list, stringVal)
% Plots true and estimated targets in 2D space
% Input arguments:
% true_targets -  true range, doppler and azimuth
% omega_vect - [N_targets x 1] vector containing target angular frequencies
% N_tx - number of transmit antennas
% doppler_vect - [N_targets x 1] vector containing target Doppler frequencies
% N_chirp - number of chirps in the transmitted signal
% detected_range - [1 x N_detected] vector containing range estimates of the detected targets
% detected_targets - [N_detected x 3] matrix containing Cartesian coordinates of the detected targets
% Output arguments:
% fig_see - handle to the figure

    fig_see = figure();
    true_x = true_targets(:,1) .* cosd(true_targets(:,3));  
    true_y =  true_targets(:,1) .* sind(true_targets(:,3));

%     true_speed = mean(true_targets(:,2));
%     det_speed = mean(detected_targets(:,2));

    plot(visible_cart_v(:,1), visible_cart_v(:,2), 'o', 'Color',[.7 .7 .7]);
    hold on;

    plot(true_x, true_y, 'o','Color','b' );
    hold on;

    for i = 1: size(detected_angles,2)
        det_x = detected_targets(:,1) .* cosd(detected_angles(:,i));
        det_y = detected_targets(:,1) .* sind(detected_angles(:,i));
       if mean(g_list) ~= 0
        marker_size = int16(g_list./max(g_list)*50+100);
       else
           marker_size = 30;
       end
        scatter(det_x, det_y, marker_size, 'ro','LineWidth', 2.5);
        hold on;
    end
    
    plot(0,0,'ko');
    xlabel('x (m)'); ylabel('y (m)');
    xlim([(min(true_x)-5), (max(true_x)+5)]);
    ylim([(min(true_y)-5), (max(true_y)+5)]);
    legend('True', 'Estimated','Radar location');
%     txt = ['True speed: ' num2str(true_speed) 'm/s'];
%     txt2 =['Estimated speed: ' num2str(det_speed) 'm/s'];
%     text(1,1,txt);
%     text(1,3,txt2);
    title('Estimated and true targets in 2D space for' + stringVal);
end
