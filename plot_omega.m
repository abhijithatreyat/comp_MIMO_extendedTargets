function plot_omega(true_targets, detected_targets,detected_angles , stringVal)
% Plot omega
fig_w = figure();
% Polar plot 
true_x= true_targets(:,1);
true_y = true_targets(:,2);
true_z = true_targets(:,3);
det_x = detected_targets(:,1);
det_y = detected_targets(:,2);
det_z = detected_angles;

plot3(true_x, true_y,true_z,'o');
hold on;
for i = 1: size(detected_angles,2)
    plot3(det_x, det_y, det_z(:,i),'x');
end
xlabel('Range'); ylabel('Doppler'); zlabel('Azimuth');
xlim([min(min(det_x),min(true_x))-2, max(max(det_x),max(true_x))+2]);
ylim([min(min(det_y),min(true_y))-2, max(max(det_y),max(true_y))+2]);
zlim([-180, 180]);

legend('True', 'Estimated');
title('True and Estimated for '+ stringVal);
end
