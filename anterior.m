function [data_matrix] =anterior(data_file, age, zone_vec)
close all; clc;

%% Read & Process Data

% Read data
data_path = strcat("data/", data_file);

% Reads all X/Y data from correct sheet and stores in matrix M. 
% First column = x, Second column = y

M = readmatrix(data_path, 'Sheet', 'Centered and Aligned', 'Range', 'A:B');
figure; scatter(M(:,2), M(:,1)); title("Original plot of data"); xlim([-5,5]); ylim([-5,5]); % Plot to confirm


% Bottom in original data, is indented from suture, so I
% replicate from top side to smoothen out data and make consistent

ant = M(M(:, 1) < 0, :); % Filter out anterior
figure; scatter(ant(:,1), ant(:,2)); title("Anterior Raw"); % Plot to check

ant_top = ant(ant(:,2) > 0, :); % Filter out top of anterior
ant_bot = [ant_top(:,1), -1*ant_top(:,2)]; % Flip across x-axis

ant_new = cat(1, ant_top, ant_bot); % Concat to form new anterior
figure; scatter(ant_new(:,1), ant_new(:,2)); title("Anterior w/ Fixed Suture") % Plot to check

X = ant_new(:,1); Y = ant_new(:,2); % Draw out X, Y from anterior data

% Process data for algorithm (places anterior on top, optic axis is x-axis)
X_data = Y; Y_data = -X;
figure; scatter(X_data, Y_data, 6); hold on;

%% Fourier 

%Find parameterized x, y of Fourier Curve
[x_fourier, y_fourier] = fourier_curve(age);

% Obtain x,y coords for anterior
fourier_bounds = [pi/2, 3*pi/2];
fp = fplot(x_fourier, y_fourier, fourier_bounds, 'LineWidth', 2); Y_fourierAnt = fp.YData; X_fourierAnt = fp.XData;
%plot(X_fourierAnt, -Y_fourierAnt, 'LineWidth', 2);

a_ant_fourier = max(X_fourierAnt);

% fourier_bounds = [-pi/2, pi/2]
% 
% [x_fourier, y_fourier] = fourier_fit(M(:,2), M(:,1));
% [X_fourier, Y_fourier] = fplot(x_fourier, y_fourier, fourier_bounds);
% plot(X_fourier, Y_fourier)

%% Chien - fit to raw
syms t

b0_ant = max(Y_data);
a_ant = max(X_data) + 0.0001; % add epsilon for numerical stability

[b1_ant, b3_ant] = findChienCoefficients(X_data', Y_data', a_ant, b0_ant);

x_chienAnt = a_ant*sin(t);
y_chienAnt = (b0_ant + b1_ant*t^2 + b3_ant*t^4)*cos(t);

chien_bounds = [-pi/2, pi/2];
fp = fplot(x_chienAnt, y_chienAnt, chien_bounds, 'LineWidth', 2); X_chienAnt = fp.XData; Y_chienAnt = fp.YData;
%plot(X_chienAnt, Y_chienAnt, 'LineWidth', 2);


%% Forbes
% format data to forbes specs
Y_forbes = -1*Y_data + max(Y_data); %figure; scatter(X_data, Y_forbes);

% fit forbes to data
syms rho;
[forbes_eq, Y_forbes_raw, A] = forbes(X_data', Y_forbes', 6);

forbes_reformat = -1*forbes_eq + eval(subs(forbes_eq, rho, a_ant));
forbes_eq = forbes_reformat;

fp = fplot(rho, forbes_eq, [min(X_data), max(X_data)], 'LineWidth', 2); X_forbes = fp.XData; Y_forbes = fp.YData;
%Y_forbes = -Y_forbes + max(Y_forbes); % Reconvert to standard format (matching other graphs)
%plot(X_forbes, Y_forbes, 'LineWidth', 2);

forbes_eq = subs(forbes_eq, rho, t);
% Note -- the t in the forbes equation is cartesian! Stands for x (not
% theta)
%% Ellipse
x_elipAnt = a_ant*cos(t); % in mm
y_elipAnt = b0_ant*sin(t); % in mm

elip_bounds = [0 pi];
fp = fplot(x_elipAnt, y_elipAnt, elip_bounds, 'LineWidth', 2); X_elipAnt = fp.XData; Y_elipAnt = fp.YData;
%plot(X_elipAnt, Y_elipAnt, 'LineWidth', 2);

legend("Raw", "Fourier", "Chien", "Forbes", "Ellipse")

data_matrix=zeros(10,4);
for zone_int=1: length(zone_vec)
zone=zone_vec(zone_int)/2; % optical zone [-zone, +zone]
offset_chien = abs(chien_bounds(1)) - asin(zone/a_ant); % polar - difference to come in from edges
offset_elip = -1* (abs(elip_bounds(1)) - acos(zone/a_ant));
offset_fourier = abs(fourier_bounds(1)) - asin(zone/a_ant_fourier); 


% Find curvature
k_chienAnt = findCurvature(x_chienAnt, y_chienAnt, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien);
k_elipAnt = findCurvature(x_elipAnt, y_elipAnt, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip);
k_fourierAnt = findCurvature(x_fourier, y_fourier, fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier);
k_forbes = findCurvature(t, forbes_eq, -zone, zone);

% Plot curvature
figure; hold on;
fplot(a_ant*sin(t), abs(1/k_chienAnt), chien_bounds);
fplot(a_ant*sin(t-pi/2), abs(1/k_elipAnt), elip_bounds);
fplot(a_ant_fourier*sin(t-pi), abs(1/k_fourierAnt), fourier_bounds);
fplot(t, abs(1/k_forbes));
legend("Chien", "Ellipse", "Fourier", "Forbes"); title("Curvature")

% Find smoothing energy (integral of derivative of curvature squared)
smth_chienAnt = eval(vpaintegral(diff(k_chienAnt, t, 1) ^ 2, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien));
smth_elipAnt = eval(vpaintegral(diff(k_elipAnt, t, 1) ^ 2, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip));
smth_fourierAnt = eval(vpaintegral(diff(k_fourierAnt, t, 1) ^ 2, fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier));
smth_forbes = eval(vpaintegral(diff(k_forbes, t, 1) ^ 2, -zone, zone));

% Find bending energy
[bendE_chienAnt, firstD_chienAnt, expr_chienAnt] = findBendingEnergy(x_chienAnt, y_chienAnt, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien);
[bendE_elipAnt, firstD_elipAnt, expr_elipAnt] = findBendingEnergy(x_elipAnt, y_elipAnt, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip);
[bendE_fourierAnt, firstD_fourierAnt, expr_fourierAnt] = findBendingEnergy(x_fourier, y_fourier, fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier);
[bendE_forbes, firstD_forbes, expr_forbes] = findBendingEnergy(t, forbes_eq, -zone, zone);

% Mean/Variance of RoC - variance found numerically
meanROC_chienAnt = abs(1/((chien_bounds(2)-offset_chien) - (chien_bounds(1)+offset_chien)) * eval(vpaintegral(1/k_chienAnt, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien)));
fp = fplot(1/k_chienAnt, [chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien], 'MeshDensity', 200);
yROC_chienAnt = fp.YData;
varROC_chienAnt = var(yROC_chienAnt);
valROC_chienAnt = abs(eval(subs(1/k_chienAnt, t, chien_bounds(1)+offset_chien)));


meanROC_elipAnt = abs(1/((elip_bounds(2)-offset_elip) - (elip_bounds(1)+offset_elip)) * eval(vpaintegral(1/k_elipAnt, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip)));
fp = fplot(1/k_elipAnt, [elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip], 'MeshDensity', 200);
yROC_elipAnt = fp.YData;
varROC_elipAnt = var(yROC_elipAnt);
valROC_elipAnt = abs(eval(subs(1/k_elipAnt, t, elip_bounds(1)+offset_elip)));

meanROC_fourierAnt = abs(1/((fourier_bounds(2)-offset_fourier) - (fourier_bounds(1)+offset_fourier)) * eval(vpaintegral(1/k_fourierAnt, fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier)));
fp = fplot(1/k_fourierAnt, [fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier], 'MeshDensity', 200);
yROC_fourierAnt = fp.YData;
varROC_fourierAnt = var(yROC_fourierAnt);
valROC_fourierAnt = abs(eval(subs(1/k_fourierAnt, t, fourier_bounds(1)+offset_fourier)));

% meanROC_forbesAnt = abs(1/(2*zone) * eval(vpaintegral(1/k_forbes, -zone, zone)));
%[unused, yROC_forbesAnt] = fplot(1/k_forbes, [-zone, zone], 'MeshDensity', 200);
%varROC_forbesAnt = var(yROC_forbesAnt);
valROC_forbesAnt = abs(eval(subs(1/k_forbes, t, -zone)));

% Fit

data = [X_data, Y_data];
data_fit = data(-3 < data(:,1), :);
data_fit = data_fit(data_fit(:,1) < 3, :);
X_fit = data_fit(:, 1); Y_fit = data_fit(:, 2);
figure; scatter(X_fit, Y_fit)

fit_forbes = getFit(X_fit, Y_fit, forbes_eq) * 10^6;
fit_elip = getFit(acos(X_fit ./ a_ant), Y_fit, y_elipAnt) * 10^6;
fit_chien = getFit(asin(X_fit ./ a_ant), Y_fit, y_chienAnt) * 10^6;
fit_fourier = getFit(asin(X_fit ./ a_ant_fourier), Y_fit, y_fourier) * 10^6;

data_matrix(4+zone_int,1)=valROC_chienAnt;
data_matrix(4+zone_int,2)=valROC_forbesAnt;
data_matrix(4+zone_int,3)=valROC_fourierAnt;
data_matrix(4+zone_int,4)=valROC_elipAnt;

if zone==3
    data_matrix(1:end-length(zone_vec),1)= [fit_chien, bendE_chienAnt, smth_chienAnt, varROC_chienAnt]';
    data_matrix(1:end-length(zone_vec),1)= [fit_forbes, bendE_forbes, smth_forbes, 0]'; % the 0 should be replaced when forbes is fixed
    data_matrix(1:end-length(zone_vec),1)= [fit_fourier, bendE_fourierAnt, smth_fourierAnt, varROC_fourierAnt]';    
    data_matrix(1:end-length(zone_vec),4)= [fit_elip, bendE_elipAnt, smth_elipAnt, varROC_elipAnt]';
end
        
end

end



