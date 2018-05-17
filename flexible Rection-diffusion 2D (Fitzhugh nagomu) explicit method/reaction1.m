function [ y ] = reaction1( V , U , epsilon )
%UNTITLED8 Summary of this function goes here
%   this is the reaction term of V
    b = 1.9;
    c = 0.1*2.54;
    y = 0;
  y = b*U-c*V;

end

