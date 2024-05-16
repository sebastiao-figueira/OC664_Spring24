% Jacob Plumb
% Nearshore Sediment Transport, Final Project
% Code for "A. Intra-wave velocity time series"

% Version 3

close all
clear all
clc

% read-in field data
data = load("dataset_19971022_0700EST.mat");

% establish data for gauge 1
u = data.u(:, 1); % cross-shore velocity (m/s, 2Hz) 
p = data.p(:,1); % water surface elevation (m, 2Hz) relative to zp
x = data.x(:,1); % cross-shore position of sensor
y = data.y(:,1); % alongshore position of sensor
zu = data.zu(:,1); % vertical position for u-sensor
zp = data.zp(:,1); % vertical position for p-sensor
zbed = data.zbed(:,1); % vertical position of seabed (from sonar altimeter)

% create time array for 21504 points at 2Hz. 2Hz = .5s time step
t = (0:.5:21504 / 2 - .5)';

% concatenate "u" and "t" to associate time data with velocity data
u_with_t = [u, t];

% do zero up crossing analysis
% find zero-up crossings for the current gauge's velocity data
zero_upcrossings_indices = find(u_with_t(1:end - 1, 1) <= 0 & u_with_t(2:end, 1) > 0);

% extract wave data for each zero-up crossing
all_waves = cell(1, numel(zero_upcrossings_indices));
for i = 1:numel(zero_upcrossings_indices)
    start_idx = zero_upcrossings_indices(i);
    if i < numel(zero_upcrossings_indices)
        end_idx = zero_upcrossings_indices(i + 1);
    else
        end_idx = size(u_with_t, 1);
    end
    all_waves{i} = u_with_t(start_idx:end_idx, :);
end

% now, the variable, "waves" is a cell array with all the intra-waves
% extracted into their own double array

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% do zero down crossing analysis
% find zero-up crossings for gauge 1's velocity data
for i = 1:numel(all_waves)
    wave_data = all_waves{i};
    wave_velocity = wave_data(:, 1);
    wave_time = wave_data(:, 2);

    % find indices where velocity goes from positive to negative
    zero_downcrossing_indices = find(wave_velocity(1:end - 1) > 0 & wave_velocity(2:end) <= 0);

    % store the time values corresponding to velocity zero down crossings
    zero_downcrossings(1, i) = wave_time(zero_downcrossing_indices);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% identify time associated with the velocity positive maxima for each wave

% iterate for all the intra-waves
for i = 1:numel(all_waves)
    % extract the velocity and time data for each wave
    wave_data = all_waves{i};
    wave_velocity = wave_data(:, 1);
    wave_time = wave_data(:, 2);

    % find the maximum velocity value within the wave
    [max_velocity, max_index] = max(wave_velocity);

    % store the maximum velocity value as positive maxima
    max_pos_velocity(1, i) = max_velocity;
    
    % store the time associated with the maximum velocity
    max_pos_time(1, i) = wave_time(max_index);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% identify time associated with the velocity negative maxima for each wave

% iterate for all the intra-waves
for i = 1:numel(all_waves)
    % extract the velocity and time data for each wave
    wave_data = all_waves{i};
    wave_velocity = wave_data(:, 1);
    wave_time = wave_data(:, 2);
    
    % find the minimum velocity value within the wave
    [min_velocity, min_index] = min(wave_velocity);
    
    % store the minimum velocity value as negative maxima
    max_neg_velocity(1, i) = min_velocity;
    
    % store the time associated with the minimum velocity
    max_neg_time(1, i) = wave_time(min_index);
end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% add the zero_downcrossings, maxima_negative_time, and
% maxima_positive_time

for i = 1:numel(all_waves)
    % assign zero_downcrossings if data exists
    if ~isempty(zero_downcrossings(1, i))
        all_waves{1, i}(1, 3) = zero_downcrossings(1, i);
    end
    
    % assign maxima_positive_time if data exists
    if ~isempty(max_pos_time(1, i))
        all_waves{1, i}(1, 4) = max_pos_time(1, i);
    end
    
    % Assign maxima_negative_time if data exists
    if ~isempty(max_neg_time(1, i))
        all_waves{1, i}(1, 5) = max_neg_time(1, i);
    end
end

% all_waves 1st column is velocity data
% all_waves 2nd column is time data
% all_waves 3rd column is time at zero down crossing
% all_waves 4th column is time at positive maxima
% all_waves 5th column is time at negative maxima

% figure(1)
% plot(all_waves{1, 5}(:, 2), all_waves{1, 5}(:, 1))
% xlabel('t [s]');
% ylabel('u [m/s]')
% title('Example of Intra-Wave for Gauge 1')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% parameterize T, T_c, T_t, T_cu, T_tu for each intra-wave
for i = 1:numel(all_waves)
    T(1, i) = all_waves{1, i}(end, 2) - all_waves{1, i}(1, 2);
    T_c(1, i) = all_waves{1, i}(1, 3) - all_waves{1, i}(1, 2);
    T_t(1, i) = all_waves{1, i}(end, 2) - all_waves{1, i}(1, 3);
    T_cu(1, i) = all_waves{1, i}(1, 4) - all_waves{1, i}(1, 2);
    T_tu(1, i) = abs(all_waves{1, i}(1, 5) - all_waves{1, i}(1, 3));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (4): u_x = u_w + abs(u_delta) * cos(phi);

phi = 0; % obliquely incident current making angle with wave direction
% for this project, we are only focusing on x-direction, so phi = 0.

% iterate to get time varying free-stream orbital velocity subscript, "w"
% for wave

for i = 1:numel(all_waves)
    u_w{1, i} = all_waves{1, i}(:, 1);
end

% determine steady state current velocity, "u_delta"

for i = 1:numel(all_waves)
    u_delta{1, i} = mean(u_w{1, i}) ./ cos(phi);
end

% calculate u_x, time varying orbital velocity in x-direction
for i = 1:numel(all_waves)
    u_x{1, i} = u_w{1, i} + abs(u_delta{1, i}) .* cos(phi);
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (5): u_y = abs(u_delta) * sin(phi); 

% calculate u_y, orbital velocity in y-direction

for i = 1:numel(all_waves)
    u_y{1, i} = abs(u_delta{1, i}) .* sin(phi);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (8): u_hat = sqrt((2 / T) * trapz((u_w) ^ 2));

% iterate to calculate u_hat for each wave

for i = 1:numel(all_waves)
    u_hat(1, i) = sqrt((2 ./ T(1, i) .* trapz(linspace(0, T(1, i), numel(u_w{1, i})), u_w{1, i} .^ 2)));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (8a) u_hat_c = sqrt((2 / T_c) * trapz((u_w) ^ 2));

% iterate to calculate u_hat_c for each wave

for i = 1:numel(all_waves)
    u_hat_c(1, i) = sqrt((2 ./ T_c(1, i) .* trapz(linspace(0, T_c(1, i), numel(u_w{1, i})), u_w{1, i} .^ 2)));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (8b) u_hat_t = sqrt((2 / T_t) * trapz((u_w) ^ 2));

% iterate to calculate u_hat_t for each wave

for i = 1:numel(all_waves)
    u_hat_t(1, i) = sqrt((2 ./ T_t(1, i) .* trapz(linspace(0, T_t(1, i), numel(u_w{1, i})), u_w{1, i} .^ 2)));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (9): a_hat = (u_hat * T) / (2 * pi);

% iterate to calculate a_hat for each wave
for i = 1:numel(all_waves)
    a_hat(1, i) = (u_hat(1, i) .* T(1, i)) ./ (2 * pi);
end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (10): u_tilda_cr = (1 / 2) * sqrt(2) * u_hat_c;

% iterate to calculate u_tilda_cr for each wave

for j = 1:numel(all_waves)
    u_tilda_cr(1, i) = (1 / 2) .* sqrt(2) .* u_hat_c(1, i);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (11): u_tilda_tr = (1 / 2) * sqrt(2) * u_hat_t;

% iterate to calculate u_tilda_tr

for i = 1:numel(all_waves)
    u_tilda_tr(1, i) = (1 / 2) .* sqrt(2) .* u_hat_t(1, i);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (12)
% u_crx = u_tilda_cr + abs(u_delta * cos(phi));
% u_cry = abs(u_delta * sin(phi));
% u_vec_cr = [u_crx u_cry];

% iterate to calculate u_crx and u_cry
for i = 1:numel(all_waves)
    u_crx(1, i) = u_tilda_cr(1, i) + abs(u_delta{1, i} * cos(phi));
    u_cry(1, i) = abs(u_delta{1, i} * cos(phi));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (13)
% u_trx = -u_tilda_tr + abs(u_delta * cos(phi));
% u_try = abs(u_delta * sin(phi));
% u_vec_tr = [u_trx u_try];

% iterate to calculate u_trx and u_try

for i = 1:numel(all_waves)
    u_trx(1, i) = -u_tilda_tr(1, i) + abs(u_delta{1, i} * cos(phi));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (celerity)

% iterate over all gauges to calculate mean eta (via pressure), then water depth

p_avg = mean(p);
h = p_avg - zp;

% call "wavelength" function at bottom of script to calculate the
% wavelength for each intra-wave

for i = 1:numel(all_waves)
    L_new(1, i) = wavelength(T(1, i), h);
end

% iterate to calculate c_w, celerity

for i = 1:numel(all_waves)
    c_w(1, i) = L_new(1, i) ./ T(1, i);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% wave dispersion equation and iteration for wavelength 
% input conditions: T in m/s, h in m

function [ L ] = wavelength(T, h)
% Wavelength Calculator - based on period (T) and and depth (h)
% where L = wavelength(T, h)

g = 9.81; % acceleration due to gravity in m/s^2
L0 = (g * T ^ 2) / (2 * pi);
guess = L0;
L = (g * T ^ 2) / (2 * pi) * tanh((2 * pi) * (h / guess));
diff = abs(L - guess);

    while diff > 0.01
        diff = abs(L - guess);
        guess = L + (0.5 * diff);
        L = (g * T ^ 2) / (2 * pi) * tanh((2 * pi) * (h / guess));
    end

end
