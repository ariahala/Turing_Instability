clc
clear

dx = 0.01;
N = 1 / dx;
U = zeros(N,N);

%set the boundary conditions below (anything)
for i = 1:N
    U(1,i) = 0.01*exp(i*dx*10);
    U(i,1) = 100*sin(i*dx*10);
    U(N,i) = 100*cos(sqrt(i^2 + N^2)*dx*20);
    U(i,N) = 100*sin(sqrt(i^2 + N^2)*dx*10);
end

s = 1;
for i = 1:N^2
    u(i) = i;
    v(i) = i;
    if( i < N+1 || i > (N^2 - N) || mod(i,N) == 1 ||  mod(i,N) == 0 )
        c(i) = 1;
    else
        c(i) = 4/(dx^2);
    end
    s = s+1;
end


for i = 1:(N^2)-1
    if ( i > N && i < (N^2) - N +1)
        if ( mod(i,N) ~= 0 && mod(i,N) ~= N-1)
            u(s) = i+1;
            v(s) = i;
            c(s) = -1/(dx^2);
            s = s+1;
        end
    end
end

for i = 2:(N^2)
    if ( i > N && i < (N^2) - N +1)
        if ( mod(i,N) ~= 2 && mod(i,N) ~= 1)
            u(s) = i-1;
            v(s) = i;
            c(s) = -1/(dx^2);
            s = s+1;
        end
    end
end

for i = N+1:((N^2) - N)
    if( mod(i,N) ~= 1 && mod(i,N) ~= 0 )
        u(s) = i;
        v(s) = i - N;
        c(s) = -1/(dx^2);
        s = s+1;
        u(s) = i;
        v(s) = i + N;
        c(s) = -1/(dx^2);
        s = s+1;
    end
end
b = zeros(N^2);
A = sparse(u,v,c);

%this section sets the function which you can play with in the functions
%definition
for i = 1:N
    for j = 1:N
        b((i-1)*N + j) = F(i*dx,j*dx);
    end
end
for i = 1:N
    b(i) = U(1,i);
    b((i-1)*N + 1) = U(i,1);
    b((N-1)*N + i) = U(N,i);
    b((i-1)*N + N) = U(i,N);
end


x = A\b;
for i = 1:N
    for j = 1:N
        U(i,j) = x((i-1)*N + j);
    end
end

surf(U);
