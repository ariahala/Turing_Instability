function y = Reaction_1( U,V )

    a = 150;
    ro = 13;
    K = 0.05;
    y = a - U - (ro*U.*V)./(1+U+K*(U.^2));

end

