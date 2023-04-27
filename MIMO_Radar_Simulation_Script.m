close all;
clear;
clc;
rng(4);
single_run = true;
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

%% single run for visual representation
if single_run
    PlotResultsFlag = true; % plot the pretty pictures
    N_target = 10;
    SimulateMIMORadarFrame(N_target,SimParameters,PlotResultsFlag);
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

function error_mat = SimulateMIMORadarFrame(N_target,SimParameters,PlotResultsFlag)
% Simulates one frame of target acquisition for settings provided in
% SimParameters object. If none or incorrectly provided, default settings
% used.
% N_target: number of targets present in the scene
% PlotResultsFlag: if true, will plot range-doppler heatmaps after each
% target extraction

%% initialize and simulate frame response
%%% Fields required in SimParameters object:
%%% (initialized with default values)
N_tx = 128; % number of elements in TX phased array
N_rx = 1; % number of elements in digital receiver
N_beacons = 15; % number of compressive beacons (subframes)
N_chirp = 32; % number of chirps in each subframe(Sampling rate of the ADC)
N_symb = 256;  % number of samples in a single chirp
perSymb_SNR_dB = -5; % mean per sample SNR when one transmitter is active
T_gap = 512; % duration between consecutive chirps (as multiples of symbol period)
DR = 10; % dynamic range of target signal amplitudes
c = 3e8; % speed of light
fs = 60e6; % Sampling rate of the Rx ADC ( This is B in the paper)
ts = 1/fs; % symbol period 
t_ax = 0:ts:((ts*N_symb)-ts); t_ax = t_ax.';% time axis
Tc = (ts*N_symb); % Chirp duration
BW = 500e6; % BW of the chirp 
chirp_slope = BW/Tc; % slope of chirp frequency ramp in Hz/s
fc = 60e9; lambda = c/fc; % carrier frequency, wavelength

%Car model parameters
% Choose between different models 
model = "Hawkeye" ; % Options : edge , normal , Hawkeye  (Default : Hawkeye)
range_translation = 10; % Expected amount of translation in the range axis (m)
car_size = 5; % size of the car ( change this to make the car a truck :) ) 
angle_threshold = 160; % Parameter for normal model only

% Antenna specific variables
TX_pos = [0,0,0]; % Tx antenna position (x,y,z) 
RX_pos = [1,0,0]; % Rx antenna position (x,y,z)
array_spacing = lambda/2; % antenna element spacing
array_idx = (1:N_tx)-ceil(N_tx/2);
array_size = [1,1];

% Detection variables
N_target = 8;
window = true;

%Plot variables
PlotResultsFlag = true; % plot the pretty pictures
num_subplot_rows = 2; % floor(sqrt(N_target+extra_bins));
num_subplot_cols = 2; % ceil(N_target+extra_bins/num_subplot_rows);

noise_power = 10^(-perSymb_SNR_dB/10);
%%% Range limits
Rmin = 0;
Rmax = fs*c / 2 / chirp_slope ;

%%% Fast frequency guard bins are removed
ffreq_guard_bins = min(2,floor(N_symb*0.01));
min_ffreq = ffreq_guard_bins*2*pi/N_symb;
max_ffreq = 2*pi - ffreq_guard_bins*2*pi/N_symb;

%%% doppler offsets, unit: radians per chirp duration 
% % (total time between start of two consecutive chirps)
doppler_guard_bins = min(2,floor(N_chirp*0.01));
min_doppler = -pi + doppler_guard_bins*2*pi/N_chirp;
max_doppler =  pi - doppler_guard_bins*2*pi/N_chirp;

%%% angles, represented by SPATIAL FREQUENCY on the (linear) TX array
omega_guard_bins = min(2,floor(N_tx*0.01));
min_omega = -pi + omega_guard_bins*2*pi/N_tx;
max_omega =  pi - omega_guard_bins*2*pi/N_tx;

%%% draw target delay, doppler, angles (equivalent spatial frequencies)
%%% from car model
[ sph_n, com , normals] = radar_car_model(model, range_translation, angle_threshold);
sph_vector = sph_n;

% sph_v = car_toy_model(range_translation, car_size);
% sph_vector = sph_v;

% sph_vector = datasample(sph_vector, 10); % Randomly choose 10 reflectors
N_points = numel(sph_vector(:,1));

true_range = sph_vector(:,3)'; % Range coming from (R,theta,phi) 
true_range = true_range(true_range>Rmin & true_range<= Rmax);

omega_vect = sph_vector(:,1)'; % Azimuth cmg from (R,theta,phi)

% Doppler vector generation
% COM moves in a random direction in X-Y plane
speed_to_doppler = 2*pi*(ts*(N_symb+T_gap))/lambda;
min_speed = min_doppler *speed_to_doppler; % m/s
max_speed = max_doppler *speed_to_doppler; % m/s
% Velocity is generated in a random direction (Refinement: Find the
% orientation of the car and move it in that direction
com_velocity = min_speed + (max_speed - min_speed).*rand([1,2]);
clear min_speed max_speed;
% Calculate radial velocity of each point in the point cloud
% Dot product between the velocity vector and each point
[cart_x,cart_y,cart_z] = sph2cart(sph_vector(:,1), sph_vector(:,2), sph_vector(:,3));

doppler_vect = (com_velocity*[cart_x';cart_y']);

% Range to freq
ffreq_vect = chirp_slope /c *(2* true_range); 
% ffreq_vect = min_ffreq + (max_ffreq-min_ffreq)*ffreq_vect;

% amplitude vecotor for one chirp
amp_vect = zeros(1,N_points);
amp_vect = calculate_amp_vect(TX_pos, RX_pos, cart_x, cart_y, cart_z, c, fc, N_points, t_ax, chirp_slope);

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
beacon_gains = exp(1i*omega_vect'*(0:N_tx-1))*beacon_pool; % Beacon gains: (N_points x N-beacons)

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
oversampling_symb = 8; % range FFT oversampling rate
oversampling_chirp = 2^ceil(log(N_beacons)/log(2)+1); % doppler FFT oversampling rate

%%% Window the signal 
if (window)
for beacon_i = 1:N_beacons
    for chirp_i = 1:N_chirp
        y_rx(beacon_i,chirp_i,:)  = reshape(y_rx(beacon_i,chirp_i,:),N_symb,1).* ...
            hann(length(y_rx)) ;
    end
end
%%% Window the signal in doppler domain
for beacon_i = 1:N_beacons
    for syms_i = 1:N_symb
        y_rx(beacon_i,:,syms_i)= reshape(y_rx(beacon_i,:,syms_i),N_chirp,1).*...
            hann(size(y_rx,2));
    end
end
end
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
range_axis = (0:N_symb*oversampling_symb-1)*2*pi/N_symb/oversampling_symb*ffreq_to_range/ts;

doppler_to_speed = lambda/2/pi/(ts*(N_symb+T_gap)); % scaling that translates 
    % doppler in radians per symbol period to radial speed in meters per second
speed_axis = angle(exp(-1i*(0:N_chirp*oversampling_chirp-1)/oversampling_chirp*2*pi/N_chirp))...
    *doppler_to_speed;

if PlotResultsFlag
    min_val = max(max(power_bins))/20000;
    fig_0 = figure();
    plot(ffreq_vect*ffreq_to_range,doppler_vect*doppler_to_speed,'bo','LineWidth',.75);
    hold on;
    xlim([Rmin-1,Rmax+1]); ylim([-1,1]*pi*doppler_to_speed) % Xlim 
    xlabel('range (m)'); ylabel('radial speed (m/s)');
    title('true (o) and estimated ( SR:x / NOMP:+ ) targets')
    figure();
    h = surf(range_axis,speed_axis,20*log10(power_bins));
    set(h,'edgecolor','none'); view(2);
    xlim([Rmin-1,Rmax+1]); ylim([-1,1]*pi*doppler_to_speed)
    xlabel('range (m)'); ylabel('radial speed (m/s)');
    title('aggregate bin power')
end

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

extra_bins = floor(log(N_target));
target_bins = zeros(N_target+extra_bins,2);

%%% NOMP for range and doppler estimation
%%% Input is y_rd , A_CS = S = ((1024x2048)1D x (1024x2048)1D)
%%% Output must be target_bins(:,1) = Doppler bins and  target_bins(:,2) = Range bins
RangeList = []; % Each column represents 1 beacon ... row-wise across beacons
DopplerList = []; % Each column represents 1 beacon ... row-wise across beacons
gainListRD = []; % Each column represents 1 beacon ... row-wise across beacons
residueListRD = []; % Each column represents 1 beacon ... row-wise across beacons


% Subframe level NOMP
tau = sqrt(noise_power) * N_symb * N_chirp/2; % TODO : Verify
[omegaRange, omegaDoppler, gainRD, residueRD] = extractRD(y_rx,y_rd, tau ,4,window, ...
    N_target, Rmin, Rmax, doppler_to_speed, range_axis, speed_axis,...
    oversampling_symb, oversampling_chirp);
RangeList = [RangeList,omegaRange];
DopplerList = [DopplerList,omegaDoppler ];
gainListRD = [gainListRD, gainRD];
residueListRD = [residueListRD, residueRD];


%%%% Post processiing info from 2D-NOMP
RangeList = RangeList .* ffreq_to_range ./ts;
DopplerList = DopplerList .* doppler_to_speed;

% Take the mean of 80% data variance
RangeListNOMP = trimmean(RangeList, 10, 2);
DopplerListNOMP = trimmean(DopplerList, 10, 2);
gainNOMP = trimmean(abs(gainListRD).^2,10,2);

% plot_NOMP_results(PlotResultsFlag, ffreq_vect, doppler_vect, ffreq_to_range, doppler_to_speed, ...
%     RangeList, DopplerList, gainListRD, Rmin, Rmax);

%%%%

power_residue = power_bins; % updated after each peak is identified and subtracted
y_residue = y_rx; % updated after each peak is identified and subtracted
y_rd_residue = y_rd; % updated after each peak is identified and subtracted
detected_range = [];
if PlotResultsFlag
    fig_2 = figure();
    subplot(num_subplot_rows,num_subplot_cols,1);
end 

[target_bins, detected_range] = extract_targets(N_target,power_bins, ...
    extra_bins, power_residue, N_symb, oversampling_symb,...
    N_chirp, oversampling_chirp, y_residue, y_rd_residue,...
    N_beacons, ffreq_to_range, ts, PlotResultsFlag, fig_0, fig_2,  ...
    num_subplot_rows, num_subplot_cols, Rmin, Rmax, doppler_to_speed,...
    range_axis, speed_axis, min_val, peak_power,window);

    %Plot the RD results
    figure(fig_0); hold on;
    for i_target = 1:size(target_bins,1) 
        ind_row = target_bins(i_target,1);
        ind_col = target_bins(i_target,2);
        w_doppler = angle(exp(-1i*2*pi*(ind_row-1)/N_chirp/oversampling_chirp)); 
        plot((ind_col-1)*2*pi/N_symb/oversampling_symb*ffreq_to_range/ts,w_doppler ...
            *doppler_to_speed,'rx','LineWidth',2,'MarkerSize', ...
            5+5*power_bins(ind_row,ind_col)/peak_power)
        hold on;
    end
    figure(fig_0); hold on;
    for i_target = 1:size(RangeListNOMP,1) 
        plot(RangeListNOMP(i_target),DopplerListNOMP(i_target),"m+","LineWidth",2,...
            "MarkerSize",5+5*gainNOMP(i_target)/max(gainNOMP));
    end
    hold off;
%%% now find the direction of each target (or targets inside each
%%% identified bin, for now just assuming one target in each bin)
omega_est = zeros(N_target,1);
for i_bin = 1:length(target_bins(:,1))
    %  target_bins(:,1) = Doppler bin
    %  target_bins(:,2) = Range bin
    y_b = y_rd(:,target_bins(i_bin,1),target_bins(i_bin,2));
    y_b = y_b(:);
    % Estimate doppler phase drift across all subframes/beacons
    dop_i = (target_bins(i_bin,1)-1)*2*pi/N_chirp/oversampling_chirp;
    y_b = y_b.*exp(-1i*dop_i*N_chirp*(0:N_beacons-1)');
    
    tau = norm(y_b)^2/7;
    [om_list, g_list, r_list] = extractSpectrum_MR(y_b, A_CS, tau, 8, 6);
    omega_est(i_bin) = angle(exp(1i*om_list(1)));
end
%%%%
% Estimate angle from NOMP-RD
mean_speed_estimate = 0;
NOMP_targts = size(RangeList,1);

omega_estNOMP = zeros(NOMP_targts,1);
for i_bin = 1 : NOMP_targts
    % Find the closest NOMP est. range and doppler bins in the grid
    target_binsNOMP(i_bin,2) = ...
        find(abs(RangeListNOMP(i_bin) - range_axis) == ...
        min(abs(RangeListNOMP(i_bin) - range_axis)));
    target_binsNOMP(i_bin,1) = ...
        find(abs(DopplerListNOMP(i_bin) - speed_axis) == ...
        min(abs(DopplerListNOMP(i_bin) - speed_axis)));
    % TODO replace this with FFT2(y_rx(i_beacon,:,:) @ (Range_Est,
    % Doppler_est)
    % WIP
    x_estimate = exp(-1i*2*pi*DopplerListNOMP(i_bin)*speed_to_doppler*(0:N_chirp-1)')* ...
    exp(-1i*2*pi*RangeListNOMP(i_bin)/76.8*(0:N_symb-1)) ;
    y_est = zeros(N_beacons,1);
    for i_beacon = 1 : N_beacons
        y_est(i_beacon) = y_rx(i_beacon,:) * x_estimate(:);
    end
    % Commented as this is not working for now!
    %y_b = y_est;
    %y_b = y_b.*exp(-1i*2*DopplerListNOMP(i_bin)*(0:N_beacons-1)');

    y_b = y_rd(:,target_binsNOMP(i_bin,1),target_binsNOMP(i_bin,2));
    % Estimate doppler phase drift across all subframes/beacons
    dop_i = (target_binsNOMP(i_bin,1)-1)*2*pi/N_chirp/oversampling_chirp;
    y_b = y_b.*exp(-1i*dop_i*N_chirp*(0:N_beacons-1)');
    
    tau = norm(y_b)^2/7;
    [om_list, g_list, r_list] = extractSpectrum_MR(y_b, A_CS, tau, 8, 6);
    omega_estNOMP(i_bin) = angle(exp(1i*om_list(1)));

end
%
%%%

%% find absolute error in range, doppler and spatial frequency estimates (in grid size)
% (fine tuned doppler estimation not coded yet)
detected_targets = zeros(N_target,3); % delay, doppler (coarse), spatial frequency
true_targets = [ffreq_vect(:)*ffreq_to_range, doppler_vect(:)*doppler_to_speed,...
    omega_vect(:)*N_tx/2/pi];

% Detected targets in oversampled FFT
detected_targets(:,1) = detected_range; % delay
detected_targets(:,2) = N_chirp/2/pi*...
    angle(exp(-1i*2*pi*(target_bins(:,1)-1)/N_chirp/oversampling_chirp));
detected_targets(:,3) = omega_est*N_tx/2/pi;


%Plot omega
plot_omega(true_targets, detected_targets, "Oversampled FFT");
% Plot the car image
plot_targets_2D(true_targets, omega_vect, N_tx,...
     N_chirp, detected_targets, "Oversampled FFT")

%Plot omega from NOMP
detected_targets = zeros(size(RangeListNOMP,1),3);
detected_targets(:,1) = RangeListNOMP;
detected_targets(:,2) = DopplerListNOMP;
detected_targets(:,3) = omega_estNOMP*N_tx/2/pi;
plot_omega(true_targets, detected_targets, "NOMP");
plot_targets_2D(true_targets, omega_vect, N_tx,...
     N_chirp, detected_targets, "NOMP");

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

