clear

dx = 0.05;
dt = 0.00009;
N = 2/dx;
T = 1000006;
no = 500;
d1 = 1; % the laplacian parameter for V after nondimensionalizing
d2 = 80; % the other laplace param

%creating a sparse matrice to solve imlplicitly
% Implicit Solver of Turing pattern

U = random('Normal',0,0.01,N,N);
V = random('Normal',0,0.01,N,N);
s = 1;

for i = 2:N-1
    for j= 2:N-1
        u(s) = (i-1)*N + j;
        v(s) = (i-1)*N + j;
        c(s) = (dt*4/(dx^2));
        s = s+1;
        u(s) = (i-1)*N + j;
        v(s) = (i-1)*N + j-1;
        c(s) = -dt/(dx^2);
        s = s+1;
        u(s) = (i-1)*N + j;
        v(s) = (i-1)*N + j+1;
        c(s) = -dt/(dx^2);
        s = s+1;
        u(s) = (i-1)*N + j;
        v(s) = (i)*N + j;
        c(s) = -dt/(dx^2);
        s = s+1;
        u(s) = (i-1)*N + j;
        v(s) = (i-2)*N + j;
        c(s) = -dt/(dx^2);
        s = s+1;
    end
end

% enforcing neumann boundary conditions
for i=2:N-1
    u(s) = i;
    v(s) = i;
    c(s) = -3;
    s = s+1;
    u(s) = i;
    v(s) = i+N;
    c(s) = 4;
    s = s+1;
    u(s) = i;
    v(s) = i+N+N;
    c(s) = -1;
    s = s+1;
    
    u(s) = (i-1)*N + 1;
    v(s) = (i-1)*N + 1;
    c(s) = -3;
    s = s+1;
    u(s) = (i-1)*N + 1;
    v(s) = (i-1)*N + 2;
    c(s) = 4;
    s = s+1;
    u(s) = (i-1)*N + 1;
    v(s) = (i-1)*N + 3;
    c(s) = -1;
    s = s+1;
    
    u(s) = (i-1)*N + N;
    v(s) = (i-1)*N + N;
    c(s) = -3;
    s = s+1;
    u(s) = (i-1)*N + N;
    v(s) = (i-1)*N + N-1;
    c(s) = 4;
    s = s+1;
    u(s) = (i-1)*N + N;
    v(s) = (i-1)*N + N-2;
    c(s) = -1;
    s = s+1;
    
    u(s) = (N-1)*N + i;
    v(s) = (N-1)*N + i;
    c(s) = -3;
    s = s+1;
    u(s) = (N-1)*N + i;
    v(s) = (N-2)*N + i;
    c(s) = 4;
    s = s+1;
    u(s) = (N-1)*N + i;
    v(s) = (N-3)*N + i;
    c(s) = -1;
    s = s+1;
end

u(s) = N^2;
v(s) = N^2;
c(s) = 1;
s = s+1;
u(s) = 1;
v(s) = 1;
c(s) = 1;
s = s+1;
u(s) = N;
v(s) = N;
c(s) = 1;
s = s+1;
u(s) = N^2 - N + 1;
v(s) = N^2 - N + 1;
c(s) = 1;
s = s+1;
A = sparse(u,v,c);
%A is now the negetive laplacian and neuann boundary condition of the system
A1 = d1*A;
A2 = d2*A;
for i = N+1:(N^2) - N
        if ( mod(i,N) ~= 1 && mod(i,N) ~= 0)
            A1(i,i) = A1(i,i) + 1;
            A2(i,i) = A2(i,i) + 1;
        end
end

    bU = zeros(N^2,1);
    bV = zeros(N^2,1);
for i = 2:N-1
    for j = 2:N-1
         bU((i-1)*N + j) = U(i,j);
         bV((i-1)*N + j) = V(i,j);
    end
end
    
for i = 1:T
    bU1 = bU;
    bV1 = bV;
    bU = bU1 + no*dt*Reaction_1(bU1,bV1);
    bV = bV1 + no*dt*Reaction_2(bV1,bU1);
    
    for i = 1:N
        bU(i) = 0;
        bU((i-1)*N + 1) = 0;
        bU((N-1)*N + i) = 0;
        bU((i-1)*N + N) = 0;
        bV(i) = 0;
        bV((i-1)*N + 1) = 0;
        bV((N-1)*N + i) = 0;
        bV((i-1)*N + N) = 0;
    end
    x1 = A1\bU;
    x2 = A2\bV;
    for i = 1:N
        for j = 1:N
            U(i,j) = x1((i-1)*N + j);
            V(i,j) = x2((i-1)*N + j);

        end
    end
    bU = x1;
    bV = x2;
    colormap('jet');
    imagesc(U);drawnow;
end
