%this simulation is the simulation of the fitzhugh nagomu reaction diffusion system in an axon which causes travelling wave

clear
N = 401;
T = 20000;
dt = 0.7;
dx = 20/(N-1);
x = [0:dx:20];


u(1) = 1;
v(1) = 1;
c(1) = -1/2;
for i = 2:N-1
    u(i) = i;
    v(i) = i;
    c(i) = 2;
end
u(N) = N;
v(N) = N;
c(N) = -1/2;

u(N+1) = 1;
v(N+1) = 2;
c(N+1) = 2;

for i = N+2:(2*N-1)
    u(i) = i-N;
    v(i) = i-N+1;
    c(i) = -1;
end

u(2*N) = N;
v(2*N) = N-1;
c(2*N) = 2;

for i = (2*N+1):(3*N-2)
    u(i) = i - 2*N + 1;
    v(i) = u(i) - 1;
    c(i) = -1;
end

u(3*N-1) = 1;
v(3*N-1) = 3;
c(3*N-1) = -3/2;

u(3*N) = N;
v(3*N) = N-2;
c(3*N) = -3/2;

A = 0.001*sparse(u,v,c)/(dx^2);

U = zeros(N,T);
V = zeros(N,T);
DU = 1;
DV = 1;

for i = 1:N
    U(i,1) = 0;
    V(i,1) = 0;
end

for i = 2:T
    Ua = U(:,i-1) + dt*(reaction(U(:,i-1),V(:,i-1),i)-A*U(:,i-1));
    Ub = ((3/4)*U(:,i-1) + (1/4)*Ua) + (dt/4)*(reaction(Ua,V(:,i-1),i)-A*Ua);
    U(:,i) = ((1/3)*U(:,i-1) + (2/3)*Ub ) + (2*dt/3)*(reaction(Ub,V(:,i-1),i) - A*Ub);
    Va = V(:,i-1) + dt*(reaction1(V(:,i-1),U(:,i-1)));
    Vb = ((3/4)*V(:,i-1) + (1/4)*Va) + (dt/4)*(reaction1(Va,U(:,i-1)));
    V(:,i) = ((1/3)*V(:,i-1) + (2/3)*Vb ) + (2*dt/3)*(reaction1(Vb,U(:,i-1)));
    plot(transpose(x),[V(:,i), U(:,i)]);
    ylim([-2 , 2]);
    drawnow;
end
colormap('jet')
imagesc(U)