function y = Reaction_2( V,U )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
alfa = 1.5;
b = 100;
ro = 13;
K = 0.05;
    y = alfa*(b-V) - (ro*U.*V)./(1+U+K*(U.^2));
end

