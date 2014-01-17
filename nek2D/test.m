N = 4; % order
K = 8;
CFL = .1;
dt = CFL/(K*N^2);
eps = .01;
c = 1;

% clean this part up - consolidate
kx = K;ky = K;
dx = 1/kx;dy = 1/ky;
Nq = N+1;
Nqkx = kx*(Nq-1) + 1; Nqky = ky*(Nq-1) + 1; % num global dofs along line
[Ah,Bh,Ch,Dh,z,w] = SEMhat(N); % to get z only...

[Mg, Kg, Cg, bcInds, f, u0] = assemble(N,K); % get global SEM matrices

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

% initial condtion
x0 =0.6; y0=0.3; delta = 0.10; R=(X-x0).^2+(Y-y0).^2;
U0 = exp(-((R./(delta^2)).^1)).*X.*(1-X).*Y.*(1-Y); % pulse * bubble
u0 = reshape(U0,Nqkx*Nqky,1); 

fg_t =  f + (1/dt)*Mg*u0;
%fg_t = fg + (1/dt)*Mg*u0 - Kg*u0; %explicit
A = (1/dt)*Mg + eps*Kg + c*Cg;
% bcs
A(bcInds,:) = zeros(size(A(bcInds,:)));
A(:,bcInds) = zeros(size(A(:,bcInds)));
A(bcInds,bcInds) = eye(length(bcInds));

figure
pcolor(xx,yy,U0)
ax = axis;
cax = caxis;
pause
Nsteps = 1/dt;
[L U] = lu(A);
for i = 1:Nsteps
    % implicit
    ug = U\(L\fg_t);
    fg_t = f + (1/dt)*Mg*ug; % next timestep
    fg_t(bcInds) = 0;
    
    % explicit
    %ug = dt*(1./diag(Mg)).*fg_t;
    %fg_t = fg + (1/dt)*Mg*ug - Kg*ug;
    
    if (mod(i,Nsteps/10)==0)
        pcolor(xx,yy,reshape(ug,Nqkx,Nqky))
        caxis(cax)
        colorbar
        title(['Time = ',num2str(i*dt)])
        axis(ax)        
        drawnow
    end
end

% view(90,0)
% lam = eig(Kglob);
% plot(real(lam),imag(lam),'.')
% title(['eps = ' num2str(eps) ', dt = ', num2str(dt)])

%
% starts = cumsum(starts);
% for i = 1:length(galnums)
% gid = galnums(i)
% indices(starts(gid) + count(gid)) = i % global to local
% count(gid) = count(gid) + 1
% end