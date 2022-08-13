function [bend_E, first_d, expr] = findBendingEnergy(x_eq, y_eq, a, b)
%% Finds bending energy
syms t;

% y_eq
% 
% % Take inverse of x to find t in terms of x
% x = finverse(x_eq)
% y_ofx = subs(y_eq, t, x) % Substitute t for x
% 
% y_pp = diff(y_ofx, t, 2)
% 
% expr = y_pp ^ 2 % Dropped absolute value because we are squaring
% 
% % Formula comes from Eberly
% eval(subs(x_eq,t,a))
% eval(subs(x_eq,t,b))
% bend_E = eval(vpaintegral(expr, subs(x_eq,t,a), subs(x_eq,t,b)));

%% Part 2
y_p = diff(y_eq, t, 1);
x_p = diff(x_eq, t, 1);

num = diff(y_p/x_p, t, 1);

y_pp_x = num / x_p;

expr = y_pp_x ^ 2;
bend_E = eval(vpaintegral(expr, a, b));
first_d = eval(vpaintegral((y_p/x_p) ^ 2, a, b));
end