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
p = data.p(:,1); % water surface elevation (m, 2Hz) relative to zp
x = data.x(:,1); % cross-shore position of sensor
y = data.y(:,1); % alongshore position of sensor
zu = data.zu(:,1); % vertical position for u-sensor
zp = data.zp(:,1); % vertical position for p-sensor
zbed = data.zbed(:,1); % vertical position of seabed (from sonar altimeter)

u = [
-0.242259
-0.134112
-0.0181004
0.103613
0.228319
0.352905
0.474073
0.58859
0.69355
0.786594
0.866061
0.931035
0.981328
1.01737
1.04007
1.05068
1.05062
1.0414
1.02448
1.00124
0.972919
0.940625
0.905304
0.867756
0.828648
0.788526
0.747833
0.706921
0.666069
0.625491
0.585351
0.545772
0.50684
0.468617
0.431142
0.394434
0.358503
0.323346
0.288951
0.255302
0.222376
0.19015
0.158595
0.127684
0.0973852
0.0676702
0.0385086
0.00987062
-0.0182728
-0.0459501
-0.0731886
-0.100015
-0.126454
-0.15253
-0.178266
-0.203683
-0.228802
-0.25364
-0.278214
-0.302539
-0.326629
-0.350493
-0.374141
-0.397579
-0.42081
-0.443833
-0.466646
-0.489241
-0.511604
-0.533719
-0.555561
-0.577099
-0.598294
-0.619096
-0.639444
-0.659265
-0.67847
-0.696953
-0.714586
-0.731216
-0.746666
-0.760723
-0.773138
-0.78362
-0.791829
-0.797371
-0.799791
-0.79857
-0.793118
-0.782772
-0.766801
-0.744411
-0.71476
-0.676984
-0.630241
-0.573761
-0.506922
-0.429344
-0.340987
-0.242259];

t = [
0
0.130536
0.261072
0.391608
0.522145
0.652681
0.783217
0.913753
1.04429
1.17483
1.30536
1.4359
1.56643
1.69697
1.82751
1.95804
2.08858
2.21911
2.34965
2.48019
2.61072
2.74126
2.87179
3.00233
3.13287
3.2634
3.39394
3.52448
3.65501
3.78555
3.91608
4.04662
4.17716
4.30769
4.43823
4.56876
4.6993
4.82984
4.96037
5.09091
5.22145
5.35198
5.48252
5.61305
5.74359
5.87413
6.00466
6.1352
6.26573
6.39627
6.52681
6.65734
6.78788
6.91841
7.04895
7.17949
7.31002
7.44056
7.5711
7.70163
7.83217
7.9627
8.09324
8.22378
8.35431
8.48485
8.61538
8.74592
8.87646
9.00699
9.13753
9.26807
9.3986
9.52914
9.65967
9.79021
9.92075
10.0513
10.1818
10.3124
10.4429
10.5734
10.704
10.8345
10.965
11.0956
11.2261
11.3566
11.4872
11.6177
11.7483
11.8788
12.0093
12.1399
12.2704
12.4009
12.5315
12.662
12.7925
12.923];

fig = figure();
plot(t, u); hold on 
yline(0)

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
zero_downcrossings = t(zero_downcrossings_indices);

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
    [min_velocity, min_index] = min(wave_velocity(1:end));
    
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

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% parameterize T, T_c, T_t, T_cu, T_tu for each intra-wave
for i = 1:numel(all_waves)
    T(1, i) = all_waves{1, i}(end, 2) - all_waves{1, i}(1, 2);
    T_c(1, i) = all_waves{1, i}(1, 3) - all_waves{1, i}(1, 2);
    T_t(1, i) = all_waves{1, i}(end, 2) - all_waves{1, i}(1, 3);
    T_cu(1, i) = all_waves{1, i}(1, 4) - all_waves{1, i}(1, 2);
    T_tu(1, i) = all_waves{1, i}(1, 5) - all_waves{1, i}(1, 3);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (4): u_x = u_w + abs(u_delta) * cos(phi);

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
    u_x{1, i} = u_w{1, i} + abs(u_delta(1, i));
end

keyboard

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (5): u_y = abs(u_delta) * sin(phi); 

% calculate u_y, orbital velocity in y-direction

% for i = 1:numel(all_waves)
%     u_y{1, i} = abs(u_delta{1, i}) .* sin(phi);
% end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (8): u_hat = sqrt((2 / T) * trapz((u_w) ^ 2));

% iterate to calculate u_hat for each wave

keyboard

for i = 1:numel(all_waves)
    u_hat(1, i) = sqrt((2 ./ T(1, i) .* trapz(linspace(0, T(1, i), numel(u_w{1, i})), u_w{1, i} .^ 2)));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (8a)

% iterate to calculate u_hat_c for each wave

for i = 1:numel(all_waves)
    u_hat_c(1, i) = max(detrend(u_w{1,i}, 'constant'));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (8b)

% iterate to calculate u_hat_t for each wave

for i = 1:numel(all_waves)
    u_hat_t(1, i) = min(detrend(u_w{1,i}, 'constant'));
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
    u_crx(1, i) = u_tilda_cr(1, i) + abs(u_delta(1, i));
    % u_cry(1, i) = abs(u_delta(1, i) * cos(phi));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (13)
% u_trx = -u_tilda_tr + abs(u_delta * cos(phi));
% u_try = abs(u_delta * sin(phi));
% u_vec_tr = [u_trx u_try];

% iterate to calculate u_trx and u_try

for i = 1:numel(all_waves)
    u_trx(1, i) = -u_tilda_tr(1, i) + abs(u_delta(1, i));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% eqn (celerity)

% iterate over all gauges to calculate mean eta (via pressure), then water depth

p_avg = mean(p);
h = zp - zbed + p_avg;

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
output_file = 'A_intrawave_velocity_timeseries_output_file.mat'; % file name

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
