

function amp_vect = calculate_amp_vect(TX_pos, RX_pos, cart_x, cart_y, cart_z, c, fc, N_points, t_ax , chirp_slope)
% Input:
%   TX_pos: transmitter position (3x1 vector)
%   RX_pos: receiver position (3x1 vector)
%   cart_x: x-coordinates of points (1xN_points vector)
%   cart_y: y-coordinates of points (1xN_points vector)
%   cart_z: z-coordinates of points (1xN_points vector)
%   c: speed of light
%   fc: carrier frequency
%   N_points: number of points to consider
% Output:
%   amp_vect: received signal in the antenna array (length(t_ax) x N_points matrix)

amp_vect = zeros(length(t_ax), N_points);

    for target = 1:N_points
        % Calculate the distance from each of the transmitter to point target
        d_T2P = sqrt((TX_pos(1)-cart_x(target)).^2+(TX_pos(2)-cart_y(target)).^2+ (TX_pos(3)-cart_z(target)).^2);
        % distance between the Rx antenna and the point reflector
        d_P2R = sqrt((RX_pos(1)-cart_x(target)).^2+(RX_pos(2)-cart_y(target)).^2+ (RX_pos(3)-cart_z(target)).^2);
        
        tau = (d_T2P+d_P2R)/c;
        %path_loss = repmat(1/d_T2P/d_P2R,[length(t_ax),1,1]);
        path_loss = (1/d_T2P/d_P2R);
        pt_signal = (path_loss.*exp(-1j * 2 * pi * fc * tau));%.* exp(1j * 2*pi * chirp_slope * t_matrix .* tau); % beat signal from a single point reflector
        amp_vect(:,target) = pt_signal; % summing up signals from all point reflector to get the received signal in the antenna array   
    end
end



