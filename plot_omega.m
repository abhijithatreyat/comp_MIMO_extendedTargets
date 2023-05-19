function plot_omega(true_targets, detected_targets, stringVal)
% Plot omega
fig_w = figure();
% Polar plot 
true_x= true_targets(:,1);
true_y = true_targets(:,2);
true_z = true_targets(:,3);
det_x = detected_targets(:,1);
det_y = detected_targets(:,2);
det_z = detected_targets(:,3);

plot3(true_x, true_y,true_z,'o');
hold on;
plot3(det_x, det_y, det_z,'x');
xlabel('Range'); ylabel('Doppler'); zlabel('Azimuth');
xlim([min(min(det_x),min(true_x))-2, max(max(det_x),max(true_x))+2]);
ylim([min(min(det_y),min(true_y))-2, max(max(det_y),max(true_y))+2]);
zlim([min(min(det_z),min(true_z))-2, max(max(det_z),max(true_z))+2]);

legend('True', 'Estimated');
title('True and Estimated for '+ stringVal);
end
