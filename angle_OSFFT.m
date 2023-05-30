function omega_est = angle_OSFFT(N_target, target_bins, y_rd, N_chirp, oversampling_chirp, N_beacons, A_CS, ...
    window,range_scaling)
omega_est = zeros(N_target,1);
    for i_bin = 1:length(target_bins(:,1))
        %  target_bins(:,1) = Doppler bin
        %  target_bins(:,2) = Range bin
        y_b = y_rd(:,target_bins(i_bin,1),target_bins(i_bin,2));
        y_b = y_b(:);
        % Estimate doppler phase drift across all subframes/beacons
        dop_i = (target_bins(i_bin,1)-1)*2*pi/N_chirp/oversampling_chirp;
        y_b = y_b.*exp(-1i*dop_i*N_chirp*(0:N_beacons-1)');
        r_est = target_bins(i_bin,2) * range_scaling;

        tau = norm(y_b)^2/7;
        [om_list, g_list, r_list] = extractSpectrum_MR(y_b, A_CS, tau, 8, 6, window, r_est);
        omega_est(i_bin) = angle(exp(1i*om_list(1)));
    end
end