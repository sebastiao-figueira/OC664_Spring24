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
u    = data.u(:, end); % cross-shore velocity (m/s, 2Hz) 
p    = data.p(:,end); % water surface elevation (m, 2Hz) relative to zp
x    = data.x(:,end); % cross-shore position of sensor
y    = data.y(:,end); % alongshore position of sensor
zu   = data.zu(:,end); % vertical position for u-sensor
zp   = data.zp(:,end); % vertical position for p-sensor
zbed = data.zbed(:,end); % vertical position of seabed (from sonar altimeter)

% create time array for 21504 points at 2Hz. 2Hz = .5s time step
t = (0:.5:21504 / 2 - .5)'; % time in s

% convert time array to a time_serial for purposes of plotting only
% dur = duration(2, 59, 12, 'Format', 'hh:mm');
dur = seconds(t);
time_serial = datenum(dur);

fig = figure(1);
plot(time_serial, u)
xlabel('t [hh:mm]');
ylabel('u [m/s]');
title('Measured Near-Bed Velocities at Gauge 1 (Closest to Shore)')
datetick('x', 'HH:MM', 'keeplimits');
xlim([0 .1250]);
hold on
reference_line = plot(get(gca, 'XLim'), [0 0], 'r-', 'LineWidth', 1.5);
legend('u', 'Reference Line (y = 0)', 'Location', 'northeast');
hold off;
% print(fig, 'PicName', '-dpng');
ax = gca;
ax.FontSize = 20;

fig = figure();
plot(u); hold on 
yline(0)
xlim([0,30])

% interpolating u and t
tq = t(1):0.05:t(end);
tq = tq';
uq = interp1(t,u,tq);

u = uq;
t = tq;

fig = figure();
plot(u); hold on 
yline(0)
xlim([0,300])
grid()

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

zero_downcrossings_indices = find(u_with_t(1:end - 1, 1) > 0 & u_with_t(2:end, 1) <= 0);
zero_downcrossings         = t(zero_downcrossings_indices);

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
    [min_velocity, min_index] = min(wave_velocity(2:end));
    
    % store the minimum velocity value as negative maxima
    max_neg_velocity(1, i) = min_velocity;
    
    % store the time associated with the minimum velocity
    max_neg_time(1, i) = wave_time(min_index+1);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% add the zero_downcrossings, maxima_negative_time, and
% maxima_positive_time

for i = 1:numel(all_waves)
    % assign zero_downcrossings if data exists
    if ~isempty(zero_downcrossings(i))
        all_waves{1, i}(1, 3) = zero_downcrossings(i);
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

figure();
plot(all_waves{1, 5}(:, 2), all_waves{1, 5}(:, 1), 'linewidth', 2)
xlabel('t [s]');
ylabel('u [m/s]');
title('Example of Intra-Wave for Gauge 1')
grid on
ax = gca;
ax.FontSize = 20;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% parameterize T, T_c, T_t, T_cu, T_tu for each intra-wave
for i = 1:numel(all_waves)
    T(1, i)    = all_waves{1, i}(end, 2) - all_waves{1, i}(1, 2);
    T_c(1, i)  = all_waves{1, i}(1, 3) - all_waves{1, i}(1, 2);
    T_t(1, i)  = all_waves{1, i}(end, 2) - all_waves{1, i}(1, 3);
    T_cu(1, i) = all_waves{1, i}(1, 4) - all_waves{1, i}(1, 2);
    T_tu(1, i) = all_waves{1, i}(1, 5) - all_waves{1, i}(1, 3);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (4): u_x = u_w + abs(u_delta) * cos(phi);

% phi = 0; % obliquely incident current making angle with wave direction
% for this project, we are only focusing on x-direction, so phi = 0.

% iterate to get time varying free-stream orbital velocity subscript, "w"
% for wave

for i = 1:numel(all_waves)
    u_w{1, i} = detrend(all_waves{1, i}(:, 1), 'constant');
end

% determine steady state current velocity, "u_delta"

for i = 1:numel(all_waves)
    u_delta(1, i) = mean(all_waves{1, i}(:, 1));
end

% calculate u_x, time varying orbital velocity in x-direction
for i = 1:numel(all_waves)
    u_x{1, i} = u_w{1, i} + u_delta(1, i);
%     u_x{1, i} = u_w{1, i} + abs(u_delta(1, i));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (5): u_y = abs(u_delta) * sin(phi); 

% calculate u_y, orbital velocity in y-direction

% for i = 1:numel(all_waves)
%     u_y{1, i} = abs(u_delta{1, i}) .* sin(phi);
% end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (8): u_hat = sqrt((2 / T) * trapz((u_w) ^ 2));

% iterate to calculate u_hat for each wave

for i = 1:numel(all_waves)
    u_hat(1, i) = sqrt((2 ./ T(1, i) .* trapz(linspace(0, T(1, i), numel(u_w{1, i})), u_w{1, i} .^ 2)));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% iterate to calculate u_hat_c for each wave

for i = 1:numel(all_waves)
    u_hat_c(1, i) = max(u_x{1,i}) - u_delta(1, i);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% iterate to calculate u_hat_t for each wave

for i = 1:numel(all_waves)
%     u_hat_t(1, i) = min(u_w{1,i});
    u_hat_t(1, i) = u_delta(1,i) - min(u_x{1,i});
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

for i = 1:numel(all_waves)
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
    u_crx(1, i) = u_tilda_cr(1, i) + u_delta(1, i);
%     u_crx(1, i) = u_tilda_cr(1, i) + abs(u_delta(1, i));
    % u_cry(1, i) = abs(u_delta(1, i) * cos(phi));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (13)
% u_trx = -u_tilda_tr + abs(u_delta * cos(phi));
% u_try = abs(u_delta * sin(phi));
% u_vec_tr = [u_trx u_try];

% iterate to calculate u_trx and u_try

for i = 1:numel(all_waves)
    u_trx(1, i) = -u_tilda_tr(1, i) + u_delta(1, i);
%     u_trx(1, i) = -u_tilda_tr(1, i) + abs(u_delta(1, i));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (celerity)

% iterate over all gauges to calculate mean eta (via pressure), then water depth

p_avg = mean(p);
h = zp - zbed + p_avg;

keyboard

% call "wavelength" function at bottom of script to calculate the
% wavelength for each intra-wave
g=9.81;
for i = 1:numel(all_waves)
    L(1, i) = 2*pi/wavenumber(T(1, i), h, g);
end

% iterate to calculate c_w, celerity

for i = 1:numel(all_waves)
    c_w(1, i) = L(1, i) ./ T(1, i);
end

% create an output file of parameters needed for other parts of the project
output_file = 'intrawave_outputs_offshore.mat'; % file name

% save the variables to a .mat file
save(output_file, 'T', 'T_c', 'T_t', 'T_cu', 'T_tu', 'u_w', 'u_delta', 'u_x', ...
    'u_hat', 'u_hat_c', 'u_hat_t', 'a_hat', 'u_tilda_cr', 'u_tilda_tr', ...
    'u_crx', 'u_trx', 'c_w');


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% This function solves the linear wave dispersion relation. This is based 
% on the two-step solution method given by Simarro & Orfila (2013). The 
% method uses an initial guess (from Beji, 2013) with a one-step correction 
% given by Newton’s method. 
%
% by: Luis Daniel Perez-Squeo
% 
% Inputs: T --- wave period [s] (scalar)
%         h --- water depth [m] (vector or scalar)
%         g --- gravitational acceleration [9.81m/s^2 or 32.17 ft/s^2]
%               (scalar)
% Output: k --- wavenumber [1/m] (same dimensions as h)

function k = wavenumber(T,h,g)

    % kh - Deep water approximation
    kh_0 = 4*(pi^2)*h./(g*(T.^2));
    
    % kh - Explicit approximation (Beji, 2013)
    kh_B = kh_0./sqrt(tanh(kh_0)) .* (1 + kh_0.^(1.09).*exp(-1.55 - 1.3*kh_0 - 0.216*(kh_0.^2)));
    
    % kh - 1-step correction with Netwon's method (Simarro & Orfila, 2013)
    kh_next = (kh_B.^2 + kh_0.*(cosh(kh_B).^2))./(kh_B + sinh(kh_B).*cosh(kh_B));
    
    % k
    k = kh_next./h;
end
