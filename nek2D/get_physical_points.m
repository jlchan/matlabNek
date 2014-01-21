function [X Y] = get_physical_points(N,kx,ky)

[Ah,Bh,Ch,Dh,z,w] = SEMhat(N); % to get z only...

Nq = N+1;
dx = 1/kx;
dy = 1/ky;
% define physical points
xx = 0; %starting x point
for i = 1:kx
    xk = dx*(z+1)/2 + dx*(i-1);
    xx = [xx xk(2:Nq)'];
end
yy = 0;
for j = 1:ky
    yk = dy*(z+1)/2 + dy*(j-1);
    yy = [yy yk(2:Nq)'];
end
[X Y] = meshgrid(xx,yy);
