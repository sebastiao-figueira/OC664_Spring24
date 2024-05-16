% Jacob Plumb
% Nearshore Sediment Transport, Final Project
% Code for "A. Intra-wave velocity time series"

close all
clear all
clc

% read-in field data
data = load("dataset_19971022_0700EST.mat");

% iterate over all gauges & assign arrays for different data types
for i = 1:8
    u{i} = data.u(:,i); % cross-shore velocity (m/s, 2Hz) 
    p{i} = data.p(:,i); % water surface elevation (m, 2Hz) relative to zp
    x(1,i) = data.x(:,i); % cross-shore position of sensors
    y(1,i) = data.y(:,i); % alongshore position of sensors
    zu(1,i) = data.zu(:,i); % vertical position for u-sensor
    zp(1,i) = data.zp(:,i); % vertical position for p-sensor
    zbed(1,i) = data.zbed(:,i); % vertical position of seabed (from sonar altimeter)
end 

% create time array for 21504 points at 2Hz. 2Hz = .5s time step
t = [0:.5:21504 / 2 - .5]';

% concatentate "t" and "u" to associate time data with velocity data
for i = 1:8
    u{1, i} = [u{1, i}, t];
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% do zero up crossing analysis
% create cell arrays to store results
all_waves = cell(1, 8);

% iterate over each gauge's velocity data
for i = 1:8
    % find zero-up crossings for the current gauge's velocity data
    zero_upcrossings_indices = find(u{i}(1:end - 1, 1) <= 0 & u{i}(2:end, 1) > 0);
    
    % extract wave data for each zero-up crossing
    waves = cell(1, numel(zero_upcrossings_indices));
    for j = 1:numel(zero_upcrossings_indices)
        start_idx = zero_upcrossings_indices(j);
        if j < numel(zero_upcrossings_indices)
            end_idx = zero_upcrossings_indices(j + 1);
        else
            end_idx = size(u{i}, 1);
        end
        waves{j} = u{i}(start_idx:end_idx, :);
    end
    
    % store extracted waves into "all_waves" cell array
    all_waves{i} = waves;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% do zero down crossing analysis 
% find time associated with the velocity zero down crossing for each wave
% establish cell array to store results
zero_downcrossings = cell(1, 8);

% iterate over each gauge's waves
for i = 1:8
    % start zero down crossing for each gauge
    zero_downcrossings{i} = cell(1, numel(all_waves{i}));
    
    % iterate over each wave for each gauge
    for j = 1:numel(all_waves{i})
        % extract the velocity and time data for each wave
        wave_data = all_waves{i}{j};
        wave_velocity = wave_data(:, 1);
        wave_time = wave_data(:, 2);
        
        % find indices where velocity goes from positive to negative
        zero_downcrossing_indices = find(wave_velocity(1:end - 1) > 0 & wave_velocity(2:end) <= 0);
        
        % store the time values corresponding to velocity zero down crossings
        zero_downcrossings{i}{j} = wave_time(zero_downcrossing_indices);
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% identify time associated with the velocity positive maxima for each wave
% establish cell arrays to store values
maxima_positive_velocity = cell(1, 8);
maxima_positive_time = cell(1, 8);

% iterate over each gauge's waves
for i = 1:8
    % start positive maxima and associated times for each gauge
    maxima_positive_velocity{i} = zeros(1, numel(all_waves{i}));
    maxima_positive_time{i} = cell(1, numel(all_waves{i}));
    
    % iterate over each wave for each gauge
    for j = 1:numel(all_waves{i})
        % extract the velocity and time data for each wave
        wave_data = all_waves{i}{j};
        wave_velocity = wave_data(:, 1);
        wave_time = wave_data(:, 2);
        
        % find the maximum velocity value within the wave
        [max_velocity, max_index] = max(wave_velocity);
        
        % store the maximum velocity value as positive maxima
        maxima_positive_velocity{i}(j) = max_velocity;
        
        % store the time associated with the maximum velocity
        maxima_positive_time{i}{j} = wave_time(max_index);
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% identify time associated with the velocity negative maxima for each wave
% establish cell arrays to store values
maxima_negative_velocity = cell(1, 8);
maxima_negative_time = cell(1, 8);

% iterate over each gauge's waves
for i = 1:8
    % start negative maxima and associated times for the current gauge
    maxima_negative_velocity{i} = zeros(1, numel(all_waves{i}));
    maxima_negative_time{i} = cell(1, numel(all_waves{i}));
    
    % iterate over each wave for each gauge
    for j = 1:numel(all_waves{i})
        % extract the velocity and time data for each wave
        wave_data = all_waves{i}{j};
        wave_velocity = wave_data(:, 1);
        wave_time = wave_data(:, 2);
        
        % find the minimum velocity value within the wave
        [min_velocity, min_index] = min(wave_velocity);
        
        % store the minimum velocity value as negative maxima
        maxima_negative_velocity{i}(j) = min_velocity;
        
        % store the time associated with the minimum velocity
        maxima_negative_time{i}{j} = wave_time(min_index);
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% add the zero_downcrossings, maxima_negative_time, and
% maxima_positive_time
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        % assign zero_downcrossings if data exists
        if ~isempty(zero_downcrossings{i}{j})
            all_waves{1, i}{1, j}(1, 3) = zero_downcrossings{1, i}{1, j};
        end
        
        % assign maxima_positive_time if data exists
        if ~isempty(maxima_positive_time{i}{j})
            all_waves{1, i}{1, j}(1, 4) = maxima_positive_time{1, i}{1, j};
        end
        
        % Assign maxima_negative_time if data exists
        if ~isempty(maxima_negative_time{i}{j})
            all_waves{1, i}{1, j}(1, 5) = maxima_negative_time{1, i}{1, j};
        end
    end
end

% all_waves 1st column is velocity data
% all_waves 2nd column is time data
% all_waves 3rd column is time at zero down crossing
% all_waves 4th column is time at positive maxima
% all_waves 5th column is time at negative maxima

% figure(1)
% plot(all_waves{1, 1}{1, 2}(:, 2), all_waves{1, 1}{1, 2}(:, 1))

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% parameterize T, T_c, T_t, T_cu, T_tu for each intra-wave
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        T{1, i}{1, j} = all_waves{1, i}{1, j}(end, 2) - all_waves{1, i}{1, j}(1, 2);
        T_c{1, i}{1, j} = all_waves{1, i}{1, j}(1, 3) - all_waves{1, i}{1, j}(1, 2);
        T_t{1, i}{1, j} = all_waves{1, i}{1, j}(end, 2) - all_waves{1, i}{1, j}(1, 3);
        T_cu{1, i}{1, j} = all_waves{1, i}{1, j}(1, 4) - all_waves{1, i}{1, j}(1, 2);
        T_tu{1, i}{1, j} = abs(all_waves{1, i}{1, j}(1, 5) - all_waves{1, i}{1, j}(1, 3));
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (4): u_x = u_w + abs(u_delta) * cos(phi);

phi = 20; % obliquely incident current making angle with wave direction

% iterate to get time varying free-stream orbital velocity subscript, "w"
% for wave
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        u_w{1, i}{1, j} = all_waves{1, i}{1, j}(:, 1);
    end
end

% determine steady state current velocity, "u_delta"
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        u_delta{1, i}{1, j} = mean(u_w{1, i}{1, j}) ./ cos(phi);
    end
end

% calculate u_x, time varying orbital velocity in x-direction
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        u_x{1, i}{1, j} = u_w{1, i}{1, j} + abs(u_delta{1, i}{1, j}) .* cos(phi);
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (5): u_y = abs(u_delta) * sin(phi); 

% calculate u_y, orbital velocity in y-direction
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        u_y{1, i}{1, j} = abs(u_delta{1, i}{1, j}) .* sin(phi);
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (8): u_hat = sqrt((2 / T) * trapz((u_w) ^ 2));

% iterate to calculate u_hat for each wave
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        u_hat{1, i}{1, j} = sqrt((2 ./ T{1, i}{1, j}) .* trapz(linspace(0, T{1, i}{1, j}, numel(u_w{1, i}{1, j})), u_w{1, i}{1, j} .^ 2));
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (8a) u_hat_c = sqrt((2 / T_c) * trapz((u_w) ^ 2));

% iterate to calculate u_hat_c for each wave
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        u_hat_c{1, i}{1, j} = sqrt((2 ./ T_c{1, i}{1, j}) .* trapz(linspace(0, T_c{1, i}{1, j}, numel(u_w{1, i}{1, j})), u_w{1, i}{1, j} .^ 2));
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (8b) u_hat_t = sqrt((2 / T_t) * trapz((u_w) ^ 2));

% iterate to calculate u_hat_t for each wave
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        u_hat_t{1, i}{1, j} = sqrt((2 ./ T_t{1, i}{1, j}) .* trapz(linspace(0, T_t{1, i}{1, j}, numel(u_w{1, i}{1, j})), u_w{1, i}{1, j} .^ 2));
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (9): a_hat = (u_hat * T) / (2 * pi);

% iterate to calculate a_hat for each wave
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        a_hat{1, i}{1, j} = (u_hat{1, i}{1, j} .* T{1, i}{1, j}) ./ (2 * pi);
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (10): u_tilda_cr = (1 / 2) * sqrt(2) * u_hat_c;

% iterate to calculate u_tilda_cr for each wave
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        u_tilda_cr{1, i}{1, j} = (1 / 2) .* sqrt(2) .* u_hat_c{1, i}{1, j};
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (11): u_tilda_tr = (1 / 2) * sqrt(2) * u_hat_t;

% iterate to calculate u_tilda_tr
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        u_tilda_tr{1, i}{1, j} = (1 / 2) .* sqrt(2) .* u_hat_t{1, i}{1, j};
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (12)
% u_crx = u_tilda_cr + abs(u_delta * cos(phi));
% u_cry = abs(u_delta * sin(phi));
% u_vec_cr = [u_crx u_cry];

% iterate to calculate u_crx and u_cry
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        u_crx{1, i}{1, j} = u_tilda_cr{1, i}{1, j} + abs(u_delta{1, i}{1, j} * cos(phi));
        u_cry{1, i}{1, j} = abs(u_delta{1, i}{1, j} * cos(phi));
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (13)
% u_trx = -u_tilda_tr + abs(u_delta * cos(phi));
% u_try = abs(u_delta * sin(phi));
% u_vec_tr = [u_trx u_try];

% iterate to calculate u_trx and u_try
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        u_trx{1, i}{1, j} = -u_tilda_tr{1, i}{1, j} + abs(u_delta{1, i}{1, j} * cos(phi));
        u_try{1, i}{1, j} = abs(u_delta{1, i}{1, j} * sin(phi));
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (celerity)

% iterate over all gauges to calculate mean eta (via pressure), then water depth
for i = 1:8
    p_avg(1, i) = mean(p{i});
    h(1, i) = p_avg(1, i) - zp(1, i);
end

% call "wavelength" function at bottom of script to calculate the
% wavelength for each intra-wave
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        L_new{1, i}{1, j} = wavelength(T{1, i}{1, j}, h(1, i));
    end
end

% iterate to calculate c_w, celerity
for i = 1:8
    for j = 1:numel(all_waves{1, i})
        c_w{1, i}{1, j} = L_new{1, i}{1, j} ./ T{1, i}{1, j};
    end
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
