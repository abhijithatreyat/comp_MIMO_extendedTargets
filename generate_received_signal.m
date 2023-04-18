function y_rx = generate_received_signal(N_beacons, N_chirp, N_symb, ...
    ts, ffreq_vect, amp_vect, doppler_vect, beacon_gains)
% Initialize the received signal
y_rx = zeros(N_beacons, N_chirp, N_symb);

% Get the number of points
N_points = length(ffreq_vect);

for target = 1:N_points
    % Calculate the frequency delay for this target
    f_delay = ts*ffreq_vect(target);

    % Calculate the nominal chirp response for this target
    y_nom_chirp = amp_vect(target)*exp(1i*f_delay*(0:N_symb-1));

    % Calculate the phase change over chirps due to Doppler (in one subframe)
    dop_phase = exp(-1i*doppler_vect(target)*(0:N_chirp-1)');

    % Calculate the nominal subframe response for this target
    y_nom = dop_phase*y_nom_chirp;

    % Calculate the phase change over subframes due to Doppler (in one frame)
    beac_dop_phase = exp(-1i*doppler_vect(target)*N_chirp*(0:N_beacons-1));

    for beacon = 1:N_beacons
        % Apply beacon and Doppler modulation on subframes
        y_tot_target = beac_dop_phase(beacon)*beacon_gains(target,beacon)...
            *reshape(y_nom,1,N_chirp,N_symb);
        % Add this target's signal to the overall received signal
        y_rx(beacon,:,:) = y_rx(beacon,:,:) + y_tot_target;
    end
end
end
