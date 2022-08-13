function [variation] = findVariation(x, y, a, b)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Finds energy
syms t;
% pre-calculating first/second order derivatives of x(t), y(t)
x_p = diff(x, t, 1);
x_pp = diff(x, t, 2);
y_p = diff(y, t, 1);
y_pp = diff(y, t, 2);

% formula for a curvature of a curve parameterized by x(t), y(t) -
% wikipedia
k = (x_p*y_pp - y_p*x_pp) / (x_p^2 + y_p^2)^(3/2);

%% Finds variation
norm = (x_p^2 + y_p^2) ^ (1/2);
d_k = diff(k, t, 1);
integrand = (d_k * (1/norm))^2 * norm;
variation = eval(vpaintegral(integrand, a, b));

end
