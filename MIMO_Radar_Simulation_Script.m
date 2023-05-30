close all;
clear;
clc;
rng(3);
single_run = false;
monte_carlo_run = false;

SimParameters = struct(); % if empty will use some default parameters
SimParameters.N_tx = 128; % number of elements in TX phased array
SimParameters.N_rx = 1; % number of elements in digital receiver
SimParameters.N_beacons = 15; % number of compressive beacons (subframes)
SimParameters.N_chirp = 32; % number of chirps in each subframe
SimParameters.N_symb = 256;  % number of samples in a single chirp
SimParameters.perSymb_SNR_dB = -5; % mean per sample SNR when one transmitter is active
SimParameters.T_gap = 512; % duration between consecutive chirps (as multiples of symbol period)
SimParameters.DR = 10; % dynamic range of target signal amplitudes

%%

%Car model parameters
% Choose between different models 
model = "normal" ; % Options : edge , normal , Hawkeye, Toy  (Default : Hawkeye)
range_translation = 20; % Expected amount of translation in the range axis (m)
                        % Keep below 70m
car_size = 10; % size of the car ( change this to make the car a truck :) ) 
angle_threshold = 160; % Parameter for normal model only

variables;
com_velocity = [];



if  model == "Toy"
    sph_v = car_toy_model(range_translation, car_size);
    sph_vector = sph_v;
else
    [ sph_n, com , car_scene_v, bbox] = radar_car_model(model, range_translation, angle_threshold);
    visible_cart_v = car_scene_v.cart_v;
    sph_vector = sph_n;
end

% single run for visual representation
if single_run
    PlotResultsFlag = true; % plot the pretty pictures
    N_target = 10;
    
    SimulateMIMORadarFrame(sph_vector,visible_cart_v,...
        N_target,model,SimParameters,PlotResultsFlag, single_run,com_velocity);
else
    PlotResultsFlag = true; % plot the pretty pictures
    N_target = 10;
    estimate = [];
    for runs = 0 :0.1:.5
        translate_x = 0.25* runs; 
        translate_y = (-6.4)* runs; 
        translate_z = 0;
        sph_v = [];
        
        translate = [translate_x, translate_y, translate_z];
        car_scene_v.cart_v = car_scene_v.cart_v + translate;
        bbox = car_scene_v.bbox + translate;
        visible_v = remove_occlusion(car_scene_v);

        blb_cart_v = model_point_reflector(visible_v,bbox, false);
        [sph_v(:,1),sph_v(:,2),sph_v(:,3)] = cart2sph(blb_cart_v(:,1),blb_cart_v(:,2),blb_cart_v(:,3));

        estimate = [estimate; SimulateMIMORadarFrame(sph_v,visible_cart_v,...
         N_target,model,SimParameters,PlotResultsFlag, single_run,com_velocity)];
    end
    save('multi-run');
end


%% run monte carlo to get error cdf
if monte_carlo_run
    N_tot_targets = 20; % total number of targets to be simulated
    min_targets = 5; % minimum number of targets in each realization
    max_targets = 8; % maximum number of targets in each realization
    PlotResultsFlag = false; % don't plot the pretty pictures
    error_mat_master = zeros(N_tot_targets,3); % contains absolute errors in range, doppler, spatial frequency
    % (in units of corresponding FFT grid size)
    number_of_simulated_targets = 0;
    while number_of_simulated_targets<N_tot_targets
        if N_tot_targets - number_of_simulated_targets <= max_targets
            N_target = N_tot_targets - number_of_simulated_targets;
        else
            N_target = min_targets + round(rand()*(max_targets-min_targets));
        end
        error_mat = SimulateMIMORadarFrame(N_target,SimParameters,PlotResultsFlag);
        error_mat_master(number_of_simulated_targets+1:number_of_simulated_targets+N_target,:) = error_mat;
        number_of_simulated_targets = number_of_simulated_targets + N_target;
        clc;  disp([num2str(number_of_simulated_targets),' out of ',num2str(N_tot_targets),' done.'])
    end
    figure()
    subplot(3,1,1)
    cdfplot(error_mat_master(:,1)); title('range error (normalized to FFT grid size)')
    subplot(3,1,2)
    cdfplot(error_mat_master(:,2)); title('doppler error (normalized to FFT grid size)')
    subplot(3,1,3)
    cdfplot(error_mat_master(:,3)); title('direction error (normalized to FFT grid size)')
end


% --------------------------------------------------------------------------------------------------------------
%% main function that runs sim:
% --------------------------------------------------------------------------------------------------------------

function detected_targets = SimulateMIMORadarFrame(sph_vector,visible_cart_v,...
    N_target,model,SimParameters,PlotResultsFlag, single_run, com_velocity)
% Simulates one frame of target acquisition for settings provided in
% SimParameters object. If none or incorrectly provided, default settings
% used.
% N_target: number of targets present in the scene
% PlotResultsFlag: if true, will plot range-doppler heatmaps after each
% target extraction

variables;


% sph_vector = datasample(sph_vector, 10); % Randomly choose 10 reflectors
N_points = numel(sph_vector(:,1));

true_range = sph_vector(:,3)'; % Range coming from (R,theta,phi) 
true_range = true_range(true_range>Rmin & true_range<= Rmax);

speed_to_doppler = 2*pi*(ts*(N_symb+T_gap))/lambda;
min_speed = min_doppler *speed_to_doppler; % m/s
max_speed = max_doppler *speed_to_doppler; % m/s
com_velocity = min_speed + (max_speed - min_speed)...
    .*mean(true_range)/max(true_range).*rand([1,2]);

omega_vect = sph_vector(:,1)'; % Azimuth cmg from (R,theta,phi)

% Doppler vector generation
% COM moves in a random direction in X-Y plane
speed_to_doppler = 2*pi*(ts*(N_symb+T_gap))/lambda;

% Calculate radial velocity of each point in the point cloud
% Dot product between the velocity vector and each point
[cart_x,cart_y,cart_z] = sph2cart(sph_vector(:,1), sph_vector(:,2), sph_vector(:,3));

% com_velocity = [0.2476, -6.3747] *speed_to_doppler ./ [mean(cart_x),mean(cart_y)] ;

doppler_vect = (com_velocity*[cart_x';cart_y']);

% Range to freq
ffreq_vect = chirp_slope /c *(2* true_range); 
% ffreq_vect = min_ffreq + (max_ffreq-min_ffreq)*ffreq_vect;

% amplitude vector for one chirp
amp_vect = calculate_amp_vect(TX_pos, RX_pos, cart_x, cart_y, cart_z,...
    c, fc, N_points, t_ax, array_x_m);

%%% normalize and ensure dynamic range does not exceed limit (DR)
amp_vect = amp_vect./sqrt(mean(amp_vect.^2));
while max(amp_vect)/min(amp_vect) > DR
    amp_vect = amp_vect + 0.1;
    amp_vect = amp_vect./sqrt(mean(amp_vect.^2));
end

%y_nom_chirp = amp_vect.*exp(1i*2*pi*0.05*rand(1,N_points)); % phase error on the rx signal
% with a standard deviation of 0.1

%%% generate compressive beacon weights
beacon_pool = 1/sqrt(N_tx)*exp(1i*2*pi*rand(N_tx,N_beacons));
A_CS = beacon_pool.';

%%% compute compressive beacon gain for each target direction
beacon_gains = exp(1i*(omega_vect)'*(0:N_tx-1))*beacon_pool; % Beacon gains: (N_points x N-beacons)

%%% generate raw received signal over one frame
y_rx = zeros(N_beacons,N_chirp,N_symb);
y_rx = generate_received_signal(N_beacons, N_chirp,...
    N_symb, ts, ffreq_vect, amp_vect, doppler_vect, beacon_gains);

% rx_sig = y_rx; % store noiseless version
%%% add noise
y_rx = y_rx + sqrt(noise_power/2)*(randn(size(y_rx)) + 1i*randn(size(y_rx)));
% Ok, we have our simulated measurements!
% now, let's detect

%% estimate targets from simulated measurements


%%% Window the signal 
% if (window)
% for beacon_i = 1:N_beacons
%     for chirp_i = 1:N_chirp
%         y_rx(beacon_i,chirp_i,:)  = reshape(y_rx(beacon_i,chirp_i,:),N_symb,1).* ...
%             kaiser(length(y_rx),30) ;
%     end
% end
% %%% Window the signal in doppler domain
% for beacon_i = 1:N_beacons
%     for syms_i = 1:N_symb
%         y_rx(beacon_i,:,syms_i)= reshape(y_rx(beacon_i,:,syms_i),N_chirp,1).*...
%             kaiser(size(y_rx,2),40);
%     end
% end
% end
%%% evaluate range-doppler domain signal ( 2D FFT)
y_rd = zeros(N_beacons, N_chirp*oversampling_chirp, N_symb*oversampling_symb);
for i_beacon = 1:N_beacons
    y_rd(i_beacon,:,:) = reshape(fft2(reshape(y_rx(i_beacon,:,:),N_chirp,N_symb), ...
                         N_chirp*oversampling_chirp,N_symb*oversampling_symb),1 ...
                        ,N_chirp*oversampling_chirp,N_symb*oversampling_symb);
end

% aggregate power in each range-doppler bin from all beacons
power_bins = reshape(sum(abs(y_rd).^2,1),N_chirp*oversampling_chirp, ...    
                     N_symb*oversampling_symb);
peak_power = max(max(power_bins));
ffreq_to_range = 1/chirp_slope*c/2; % scaling that translates 
                                    % delay in symbol periods to range in meters

% oversampled range-doppler axes:
range_resolution = 2*pi/N_symb/oversampling_symb*ffreq_to_range/ts;
range_axis = (0:N_symb*oversampling_symb-1)*range_resolution;

doppler_to_speed = lambda/2/pi/(ts*(N_symb+T_gap)); % scaling that translates 
    % doppler in radians per symbol period to radial speed in meters per second
doppler_resolution = doppler_to_speed/oversampling_chirp*2*pi/N_chirp;
speed_axis = angle(exp(-1i*(0:N_chirp*oversampling_chirp-1)/oversampling_chirp*2*pi/N_chirp))...
    *doppler_to_speed;
min_val = max(max(power_bins))/20000;

if PlotResultsFlag
    figure();
    h = surf(range_axis,speed_axis,20*log10(power_bins));
    set(h,'edgecolor','none'); view(2);
    xlim([Rmin-1,Rmax+1]); ylim([-1,1]*pi*doppler_to_speed)
    xlabel('range (m)'); ylabel('radial speed (m/s)');
    title('aggregate bin power');
end

% save the session
%save("preEstimation");
% 
% figure;
% h = surf(range_axis,speed_axis,20*log10(power_bins.*CFAR_mat));
% set(h,'edgecolor','none'); view(2);
% % xlim([Rmin-1,Rmax+1]); ylim([-1,1]*pi*doppler_to_speed)
% xlabel('range (m)'); ylabel('radial speed (m/s)');
% title('filtered bin power')

%% Estimation 
%%% find the (N_points + extra_bins) strongest components in range-doppler
%%% signal. This part is a crude version of NOMP crude and will be cleaned
%%% up in next version.
%load("preEstimation");
window = false;
extra_bins = floor(log(N_target));
target_bins = zeros(N_target+extra_bins,2);

%%% NOMP for range and doppler estimation
%%% Input is y_rd , A_CS = S = ((1024x2048)1D x (1024x2048)1D)
%%% Output must be target_bins(:,1) = Doppler bins and  target_bins(:,2) = Range bins
RangeList = []; % Each column represents 1 beacon ... row-wise across beacons
DopplerList = []; % Each column represents 1 beacon ... row-wise across beacons
gainListRD = []; % Each column represents 1 beacon ... row-wise across beacons
residueListRD = []; % Each column represents 1 beacon ... row-wise across beacons

% Subframe level NOMP for all targets
if model == "Toy"
    tau = sqrt(noise_power)*N_symb*N_chirp/4; 
else
    tau = sqrt(noise_power)*N_symb*N_chirp/2; 
end
[omegaRange, omegaDoppler, gainRD, residueRD] = extractRD(y_rx,y_rd, tau ,4,window, ...
    N_target, Rmin, Rmax  ,doppler_to_speed, range_axis, speed_axis,...
    oversampling_symb, oversampling_chirp, min_val, PlotResultsFlag);
RangeList = [RangeList,omegaRange];
DopplerList = [DopplerList,omegaDoppler ];
gainListRD = [gainListRD, gainRD];
residueListRD = [residueListRD, residueRD];


%%%% Post processiing info from 2D-NOMP
RangeListNOMP = RangeList .* ffreq_to_range ./ts;
DopplerListNOMP = DopplerList .* doppler_to_speed;

 true_targets = [ffreq_vect(:)*ffreq_to_range, doppler_vect(:)*doppler_to_speed,...
    omega_vect(:)*180/pi];

% Filter the data
% if model ~= "Toy"
%     % Filter out the range values more than size of the car
%     z_score = abs((RangeListNOMP - mean(RangeListNOMP)) ./ std(RangeListNOMP));
%     threshold = car_size/std(RangeListNOMP);
%     RangeListNOMP = RangeListNOMP(z_score <= threshold);
%      DopplerListNOMP = DopplerListNOMP(1:size(RangeListNOMP,1));
%      residueListRD = residueListRD(1:size(RangeListNOMP,1));
%     
%     % Filter out doppler values more than 2 m/s deviation from the mean
%     % This is done because in a practical scenario, a car would not have its
%     % points relatively move more than 2m/s in opposite directions
%     z_score = abs((DopplerListNOMP - mean(DopplerListNOMP)) ./ std(DopplerListNOMP));
%     threshold = 2/std(DopplerListNOMP);
%     DopplerListNOMP = DopplerListNOMP(z_score <= threshold);
% 
%  RangeListNOMP = RangeListNOMP(1:size(DopplerListNOMP,1));
%  residueListRD = residueListRD(1:size(DopplerListNOMP,1));
% end


% plot_NOMP_results(PlotResultsFlag, ffreq_vect, doppler_vect, ffreq_to_range, doppler_to_speed, ...
%     RangeList, DopplerList, gainListRD, Rmin, Rmax);

%%%%
if(single_run)
    power_residue = power_bins; % updated after each peak is identified and subtracted
    y_residue = y_rx; % updated after each peak is identified and subtracted
    y_rd_residue = y_rd; % updated after each peak is identified and subtracted
    detected_range = [];
    fig_2 = figure();
    if PlotResultsFlag
        subplot(num_subplot_rows,num_subplot_cols,1);
    end

    [target_bins, detected_range] = extract_targets(N_target,power_bins, ...
        extra_bins, power_residue, N_symb, oversampling_symb,...
        N_chirp, oversampling_chirp, y_residue, y_rd_residue,...
        N_beacons, ffreq_to_range, ts, PlotResultsFlag, fig_2,  ...
        num_subplot_rows, num_subplot_cols, Rmin, Rmax, doppler_to_speed,...
        range_axis, speed_axis, min_val, peak_power,window);

    %Plot the RD results
    plotRD_grid(ffreq_vect, ffreq_to_range, doppler_vect,...
        doppler_to_speed, range_resolution, doppler_resolution, ...
        target_bins, N_chirp, oversampling_chirp, power_bins,...
        peak_power, RangeListNOMP, DopplerListNOMP, residueListRD);

    omega_est = angle_OSFFT(N_target, target_bins,...
        y_rd, N_chirp, oversampling_chirp, N_beacons, A_CS, window, ...
        2*pi/N_symb/oversampling_symb*ffreq_to_range/ts);

    detected_targets = zeros(N_target,3); % delay, doppler (coarse), spatial frequency

    % Detected targets in oversampled FFT
    detected_targets(:,1) = detected_range; % delay
    detected_targets(:,2) = N_chirp/2/pi*...
        angle(exp(-1i*2*pi*(target_bins(:,1)-1)/N_chirp/oversampling_chirp));
    detected_targets(:,3) = omega_est*180/pi;
    
    
    %Plot omega
    %plot_omega(true_targets, detected_targets, "Oversampled FFT");
    % Plot the car image
    %plot_targets_2D(true_targets, omega_vect, N_tx,...
        % N_chirp, detected_targets,g_list, "Oversampled FFT")
else
    target_bins = [];
    plotRD_grid(ffreq_vect, ffreq_to_range, doppler_vect,...
        doppler_to_speed, range_resolution, doppler_resolution, ...
        target_bins, N_chirp, oversampling_chirp, power_bins,...
        peak_power, RangeListNOMP, DopplerListNOMP, residueListRD);
end
%%% now find the direction of each target (or targets inside each
%%% identified bin, for now just assuming one target in each bin)


%save('rdEstimation');
%%
% Estimate angle from NOMP-RD
%load('rdEstimation');
window = false;
[omega_estNOMP, g_list] = angle_NOMP(RangeListNOMP, ...
    range_axis, DopplerListNOMP, speed_axis, ...
    speed_to_doppler, N_chirp, N_symb, N_beacons, ...
    y_rd, oversampling_chirp, A_CS, window);
%

% Doppler Refinement 
% beacons_estimate = exp(1i*omega_estNOMP(1)*(0:N_tx-1))*beacon_pool; % (N_points x N-beacons)
% for beacon = 1:N_beacons
%         y_refine(beacon,:,:) = y_rx(beacon,:,:) * conj(beacons_estimate(:,beacon));
%     
% end
% 
% for i_beacon = 1:N_beacons
%     y_rd(i_beacon,:,:) = reshape(fft2(reshape(y_refine(i_beacon,:,:),N_chirp,N_symb), ...
%     N_chirp*oversampling_chirp,N_symb*oversampling_symb),1 ...
%     ,N_chirp*oversampling_chirp,N_symb*oversampling_symb);
% end
% power_bins = reshape(sum(abs(y_rd).^2,1),N_chirp*oversampling_chirp, ...
% N_symb*oversampling_symb);
% tau = sqrt(noise_power)*N_symb*N_chirp; 
% [omegaRange, omegaDoppler, gainRD, residueRD] = extractRD(y_refine,y_rd, tau ,4,window, ...
%     N_target, Rmin, Rmax, doppler_to_speed, range_axis, speed_axis,...
%     oversampling_symb, oversampling_chirp, min_val);
%  plotRD_grid(ffreq_vect, ffreq_to_range, doppler_vect,...
%         doppler_to_speed, range_resolution, doppler_resolution, ...
%         target_bins, N_chirp, oversampling_chirp, power_bins,...
%         peak_power, RangeListNOMP, DopplerListNOMP, gainNOMP);
% %%%

%% find absolute error in range, doppler and spatial frequency estimates (in grid size)

%Plot omega from NOMP
detected_targets = zeros(size(RangeListNOMP,1),3);

detected_targets(:,1) = RangeListNOMP;
detected_targets(:,2) = DopplerListNOMP;
detected_targets(:,3) = omega_estNOMP*180/pi;
if(PlotResultsFlag)
    plot_omega(true_targets, detected_targets, "NOMP");
    
    plot_targets_2D(true_targets,visible_cart_v, omega_vect, N_tx,...
         N_chirp, detected_targets,g_list, "NOMP");
end
% error_mat = zeros(N_target,3);
% % for each true target, look for closest match in estimated targets and
% % compute error
% for i_target = 1:N_target
%     error_vects = abs((detected_targets - true_targets(i_target,:)));
%     [~,ind_match] = min(sum(error_vects,2));
%     error_mat(i_target,:) = error_vects(ind_match,:);
% end

end


% -------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------

function plotRD_grid(ffreq_vect, ffreq_to_range, doppler_vect,...
    doppler_to_speed, range_resolution, doppler_resolution, ...
    target_bins, N_chirp, oversampling_chirp, power_bins,...
    peak_power, RangeListNOMP, DopplerListNOMP, residueRD)
fig_0 = figure();
plot(ffreq_vect*ffreq_to_range,doppler_vect*doppler_to_speed,'bo','LineWidth',.75);
hold on;
range_min = min(ffreq_vect*ffreq_to_range)-5;
range_max = max(ffreq_vect*ffreq_to_range)+5;
xlim([range_min,range_max]); 
doppler_min = min(doppler_vect*doppler_to_speed);
doppler_max = max(doppler_vect*doppler_to_speed);
ylim([doppler_min-2 , doppler_max+2]); 

if(~isempty(target_bins))
    ind_row = target_bins(:,1);
    ind_col = target_bins(:,2);
    w_doppler = angle(exp(-1i*2*pi*(ind_row-1)/N_chirp/oversampling_chirp));
    marker_size = 30 + 50 * power_bins(sub2ind(size(power_bins), ind_row, ind_col)) / peak_power;
    x = (ind_col-1)*range_resolution;
    y = w_doppler * doppler_to_speed;
    scatter(x', y', int16(marker_size'), 'rx','LineWidth', 2);
end

% w_doppler = angle(exp(-1i*2*pi*(target_bins(:,1)'-1)/N_chirp/oversampling_chirp)); 
% plot((target_bins(:,2)'-1)*range_resolution,w_doppler ...
%     *doppler_to_speed,'rx','LineWidth',2, 'MarkerSize', ...
%             5+5*power_bins(ind_row,ind_col)/peak_power);

marker_size = 30 + 50 * residueRD / max(residueRD);
scatter(RangeListNOMP, DopplerListNOMP, marker_size, 'm+', 'LineWidth', 2 );


xlabel('range (m)'); ylabel('radial speed (m/s)');
title('Range doppler Estimates');
grid on;
xticks(range_min:10*range_resolution:range_max); 
yticks(doppler_min:10*doppler_resolution:doppler_max);
hold off;

legend('True targets', 'Oversampled FFT', 'NOMP');

end