function plot_targets_2D(sph_vector, omega_vect, N_tx,...
    doppler_vect, N_chirp, detected_range, detected_targets)
% Plots true and estimated targets in 2D space
% Input arguments:
% sph_vector - [N_targets x 3] matrix containing spherical coordinates of the true targets
% omega_vect - [N_targets x 1] vector containing target angular frequencies
% N_tx - number of transmit antennas
% doppler_vect - [N_targets x 1] vector containing target Doppler frequencies
% N_chirp - number of chirps in the transmitted signal
% detected_range - [1 x N_detected] vector containing range estimates of the detected targets
% detected_targets - [N_detected x 3] matrix containing Cartesian coordinates of the detected targets
% Output arguments:
% fig_see - handle to the figure

    fig_see = figure();
    true_x = sph_vector(:,3) .* cosd(omega_vect(:)*N_tx/2/pi);
    true_y = sph_vector(:,3) .* sind(omega_vect(:)*N_tx/2/pi);
    true_speed = mean(doppler_vect(:)*N_chirp/2/pi);

    det_x = detected_range' .* cosd(detected_targets(:,3));
    det_y = detected_range' .* sind(detected_targets(:,3));
    det_speed = mean(detected_targets(:,2));

    plot(true_x, true_y, 'o');
    hold on;
    plot(det_x, det_y, 'x');
    hold on;
    plot(0,0,'ko');
    xlabel('x (m)'); ylabel('y (m)');
    xlim([-1, max(max(det_x),max(true_x))+2]);
    ylim([-2, max(max(det_y),max(true_y))+2]);
    legend('True', 'Estimated','Radar location');
    txt = ['True speed:' num2str(true_speed) 'm/s'];
    txt2 =['Estimated speed:' num2str(det_speed) 'm/s'];
    text(1,1,txt);
    text(1,1.5,txt2);
    title('Estimated and true targets in 2D space');
end
