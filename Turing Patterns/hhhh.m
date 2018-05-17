clear
[t,y] = meshgrid([-15:0.1:15],[0:0.1:30]);
s = 100;
q = 20;
g = t*0;
for i = 1:1000000
    z = ((sin((-i/5) + 2*sqrt((q-y).^2+(t+s).^2))./sqrt((q-y).^2+(t+s).^2)) + (sin((-i/5) + 2*sqrt((q-y).^2+(t-s).^2))./sqrt((q-y).^2+(t-s).^2)));
    
    g = max(abs(z),g);
        
imagesc(z);
zlim([-1 1]);
drawnow;

    
    
end