clear
N = 301;
dx = 10/(N-1);
dt = 0.03;
T = 10000;
U = zeros(N,N);
V = zeros(N,N);
 U = random('Normal' , 0 , 1 , N , N );
 V = random('Normal' , 0 , 1 , N , N );

 for u = 1:10
     for i = 2:N-1
         for j = 2:N-1
%                            U(i,j) = sin(-((i-(N/2))^2 +(j-(N/2))^2)*dx) ;
                                U(i, j) = (U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1))/4;
               V(i,j) = 0;
         end
     end
end
for i = 1:N
    U(1,i) = 0;
    U(N,i) = 0;
    U(i,1) = 0;
    U(i,N) = 0;
end
%diffusion parameter for U
D = 0.001;
diffusion parameter for V
D1 = 0;
epsilon = 1;

for i = 2:T
    
    Ua = U + dt*(reaction(U,V,i)-D*laplacian2D(U,N,N,dx));
    for j = 2:N-1
         Ua(1,j) = ((-1/2)*Ua(3,j) + 2*Ua(2,j))*2/3;
        Ua(N,j) = ((-1/2)*Ua(N-2,j) + 2*Ua(N-1,j))*2/3;
        Ua(j,1) = ((-1/2)*Ua(j,3) + 2*Ua(j,2))*2/3;
         Ua(j,N) = ((-1/2)*Ua(j,N-2) + 2*Ua(j,N-1))*2/3;
    end
    
    Ub = ((3/4)*U + (1/4)*Ua) + (dt/4)*(reaction(Ua,V,i)-D*laplacian2D(Ua,N,N,dx));
    for j = 2:N-1
         Ub(1,j) = ((-1/2)*Ub(3,j) + 2*Ub(2,j))*2/3;
        Ub(N,j) = ((-1/2)*Ub(N-2,j) + 2*Ub(N-1,j))*2/3;
        Ub(j,1) = ((-1/2)*Ub(j,3) + 2*Ub(j,2))*2/3;
         Ub(j,N) = ((-1/2)*Ub(j,N-2) + 2*Ub(j,N-1))*2/3;
    end
    
    Uc = ((1/3)*U + (2/3)*Ub ) + (2*dt/3)*(reaction(Ub,V,i) - D*laplacian2D(Ub,N,N,dx));
    
     Va = V + dt*(reaction1(V,U,epsilon)- D1*laplacian2D(Ub,N,N,dx));
     
     for j = 2:N-1
         Va(1,j) = ((-1/2)*Va(3,j) + 2*Va(2,j) )*2/3;
        Va(N,j) = ((-1/2)*Va(N-2,j) + 2*Va(N-1,j))*2/3;
        Va(j,1) = ((-1/2)*Va(j,3) + 2*Va(j,2))*2/3;
         Va(j,N) = ((-1/2)*Va(j,N-2) + 2*Va(j,N-1))*2/3;
     end
     
     Vb = ((3/4)*V + (1/4)*Va) + (dt/4)*(reaction1(Va,U,epsilon)- D1*laplacian2D(Ub,N,N,dx));
     Vb(1,1) = 0;
     Vb(1,N) = 0;
     Vb(N,1) = 0;
     Vb(N,N) = 0;
     for j = 1:N
         Vb(1,j) = ((-1/2)*Vb(3,j) + 2*Vb(2,j) )*2/3;
        Vb(N,j) = ((-1/2)*Vb(N-2,j) + 2*Vb(N-1,j))*2/3;
        Vb(j,1) = ((-1/2)*Vb(j,3) + 2*Vb(j,2))*2/3;
         Vb(j,N) = ((-1/2)*Vb(j,N-2) + 2*Vb(j,N-1))*2/3;
     end
     Vc = ((1/3)*V + (2/3)*Vb) + (2*dt/3)*(reaction1(Vb,U,epsilon)- D1*laplacian2D(Ub,N,N,dx));
     U = Uc;
     V = Vc;
     for j = 2:N-1
         U(1,j) = ((-1/2)*U(3,j) + 2*U(2,j))*2/3;
        U(N,j) = ((-1/2)*U(N-2,j) + 2*U(N-1,j))*2/3;
        U(j,1) = ((-1/2)*U(j,3) + 2*U(j,2))*2/3;
         U(j,N) = ((-1/2)*U(j,N-2) + 2*U(j,N-1))*2/3;
         V(1,j) = ((-1/2)*V(3,j) + 2*V(2,j) )*2/3;
        V(N,j) = ((-1/2)*V(N-2,j) + 2*V(N-1,j))*2/3;
        V(j,1) = ((-1/2)*V(j,3) + 2*V(j,2))*2/3;
         V(j,N) = ((-1/2)*V(j,N-2) + 2*V(j,N-1))*2/3;
     end
     
     colormap('jet')
     imagesc(U);
     zlim([-10,10]);
     drawnow;
end
