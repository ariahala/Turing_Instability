% PDE general non-linear solver Godunov method
% Aria Halvati
%The initial condition can be set to anything and the transition function
%is also arbitrary
clc
clear
dT = 0.02;
dX = 0.1;
xLength = 16;
tLength = 10;
len = xLength/dX;
U = zeros(xLength/dX,tLength/dT);
for i = 1:xLength/dX
    if ( i*dX - 1 >= 0 && i*dX - 1 <= pi )
        U(i,1) = sin(i*dX - 1);
    end
end

for i = 1:(tLength/dT)-1
    for j = 1:len
        if ( j == len )
            q = 1;
        else
            q = j+1;
        end
        if ( U(j,i) >= U(q,i) )
            k = [U(q,i):0.001:U(j,i)];
            F(j) = max(f(k));
        else
            k = [U(j,i):0.001:U(q,i)];
            F(j) = min(f(k));
        end
    end
    for j = 1:len
        if ( j == 1 )
            U(j,i+1) = U(j,i) - (dT/dX)*(F(j) - F(len));
        else
            U(j,i+1) = U(j,i) - (dT/dX)*(F(j) - F(j-1));
        end
    end
      
    plot(U(:,i));
      xlim( [-1 150] )
    ylim ( [-1 1] )
    drawnow;

        
end