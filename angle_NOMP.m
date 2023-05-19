


function omega_estNOMP = angle_NOMP( RangeListNOMP, range_axis, DopplerListNOMP, speed_axis, speed_to_doppler, N_chirp, N_symb, N_beacons, y_rx, y_rd, oversampling_chirp, A_CS)

NOMP_targts = size(RangeListNOMP,1);
omega_estNOMP = zeros(NOMP_targts,1);
% Estimate the doppler drift on average
%dop_i = mean((DopplerListNOMP)*2*pi/N_chirp/oversampling_chirp);

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
%         x_estimate = exp(-1i*2*pi*DopplerListNOMP(i_bin)*speed_to_doppler*(0:N_chirp-1)')* ...
%         exp(-1i*2*pi*RangeListNOMP(i_bin)/76.8*(0:N_symb-1)) ;
%         y_est = zeros(N_beacons,1);
%         for i_beacon = 1 : N_beacons
%             y_est(i_beacon) = y_rx(i_beacon,:) * x_estimate(:);
%         end
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
end