function [ y ] = reaction1( V , U , epsilon )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here   
    b = 0.008;
    c = 0.008*2.54;
    y = 0;
  y = b*U-c*V;

end

