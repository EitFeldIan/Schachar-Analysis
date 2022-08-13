function [E, k] = findEnergy(x, y, a, b)
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

% expression to be integrated -- from Horn pg. 5
norm = (x_p^2 + y_p^2) ^ (1/2);
expr = k^2 * norm;

eval(vpaintegral(norm, a, b))

E = eval(vpaintegral(expr, a, b)) * eval(vpaintegral(norm, a, b));

% figure; hold on; fplot(x); fplot(y); title("x & y as a function of t"); legend('x', 'y')
% figure; hold on; fplot(x_p); fplot(y_p); title("x' & y' as a function of t"); legend('x', 'y')
end

