function [ y ] = reaction( U , V , t )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   for seeing monostability set a=0.135 , for bi stability set a=-0.1
    a = 0.135;
    o = size(U);
    I = zeros(o);
    if ( t < 6 )
    for i = 1 : 10
        I(i) = 1;
    end
    end
    
%     y = U.*(1-U).*(U-a) + U - 4*V + I;\
y = U.*(a-U).*(U-1) - V + I;
end

