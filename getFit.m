function [fit] = getFit(raw_x, raw_y, eq)
% Finds fit of equation to raw data
syms t
syms y(t)

y(t) = eq;

eq_y = double(y(raw_x)); % Value of equation evaluated at all points

% Calculate mean-sqaured fit

fit = mean((eq_y - raw_y).^2);
end

