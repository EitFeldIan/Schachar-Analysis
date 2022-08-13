function [z, Y_forbes, A] = forbes(X_Ant, Y_Ant, M)
% Please ensure curve is concave up & centered at the origin

%% Forbes polynomial

% Curvature of best-fit circle
syms t;

c = (2*max(Y_Ant) / (max(X_Ant)^2 + max(Y_Ant)^2)); 
r = 1/c; % radius of curvature

%% OLS Formulation
% Y
Y = Y_Ant'; % (d x 1)

% S
rho = X_Ant;
u = rho ./ max(rho);

lambda = (c * rho.^2) ./ (1+sqrt(1-c^2 .* rho.^2));
beta = (u.^2 .* (1 - u.^2)) ./ (sqrt(1-c^2.*rho.^2));

mini_S = (u.^2 .* ones(M+1, size(u,2))) .^ repmat((0:M)', 1, size(u,2));
mini_S = mini_S .* beta;

S = [lambda ; mini_S];

% A
A = inv(S * S') * S * Y;

%% Assembling curve 

% Points
Y_forbes = A'*S;

% Equation (rho is equivalent to x, z is equivalent to y)
syms rho
u = rho / max(X_Ant);

z = (c*rho^2) / (1+sqrt(1-c^2*rho^2)) + (u^2 * (1-u^2)) / (sqrt(1-c^2*rho^2)) .* (A(2:end)' * (u.^2 .^ (0:M)'));
end

