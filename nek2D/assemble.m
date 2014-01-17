% N = 10; % order
% K = 16; % num elems

function [Mg, Kg, Cg, bcInds, fg, u0] = assemble(N,K)

if (nargin<2)
   N = 4;
   K = 4;
end

Kx = K;
Ky = K;
Nq = N+1;
Nq2 = Nq*Nq;
dx = 1/Kx;
dy = 1/Ky;

galnums = zeros(Kx*Ky*Nq2,1);
numStarts = (Kx*N+1)*(Ky*N+1)+1;
starts = zeros(numStarts,1);
count = zeros(numStarts,1);
indices = zeros(Kx*Ky*Nq2,1);
yOff = Kx*N + 1; % offset going up in the y direction
for kx = 1:Kx
    for ky = 1:Ky
        off = (kx-1)*N + ((ky-1)*N)*yOff;
        for i = 1:Nq
            for j = 1:Nq
                gid = off + (i-1) + (j-1)*yOff + 1;
                elem_offset = ((ky-1)*Kx + (kx-1))*Nq2;
                local_offset = ((i-1) + (j-1)*Nq) + 1;
                galnums(elem_offset + local_offset) = gid;
                %starts(gid) = starts(gid) + 1; % should be global to local map
            end
        end
    end
end
numdofs = gid;

% make 2D element matrices
%dt = .5*min(dx,dy)/N^2; % cfl condition
dt = .01;
eps = 0.01;

% bb = [1 1];
% %cx = -cos(pi2*xx).*sin(pi2*yy); cy =  sin(pi2*xx).*cos(pi2*yy);
% a = @(x) ones(size(x))*bb(1); % define beta = (ab,cd)
% b = @(y) ones(size(y))*bb(1);
% c = @(x) ones(size(x))*bb(2);
% d = @(y) ones(size(y))*bb(2);

a = @(x) ones(size(x));
b = @(y) 5*(y-.5);
c = @(x) -5*(x-.5);
d = @(y) ones(size(y));

% some nonconst nonzero forcing
Fx = @(x) .5*exp(x);
Fy = @(y) y.^2;

% 1D matrices and GLL points
[Ah,Bh,Ch,Dh,z,w] = SEMhat(N);       

% preassemble element mass/stiffness matrices
Jx = dx/2; Jy = dy/2;
M = zeros(Nq2,Nq2);
K = zeros(Nq2,Nq2);
for i = 1:Nq
    for j = 1:Nq
        r = i + Nq*(j-1);                
        % diagonal mass
        M(r,r) = Jx*Jy*w(i)*w(j);
        
        % delta_jl -> j = l, loop over k
        for k = 1:Nq
            l = j;
            q = k + Nq*(l-1);
            K(r,q) = K(r,q) + ...
                Ah(i,k)*w(j);
            
        end
        % delta_ik -> i = k, loop over l
        for l = 1:Nq
            k = i;
            q = k + Nq*(l-1);
            K(r,q) = K(r,q) + ...
                Ah(j,l)*w(i);
        end
        
    end
end

% assemble global mass/stiffness/convection (also local)
Mg = sparse(numdofs,numdofs);
Kg = sparse(numdofs,numdofs);
Cg = sparse(numdofs,numdofs);
fg = zeros(numdofs,1);
for kx = 1:Kx
    for ky = 1:Ky        
        Jx = dx/2; Jy = dy/2;
        % assemble element matrices
        C = sparse(Nq2,Nq2);
        xp = dx*(z+1)/2 + dx*(kx-1);
        yp = dy*(z+1)/2 + dy*(ky-1);
        for i = 1:Nq
            for j = 1:Nq
                r = i + Nq*(j-1);                
                % delta_jl -> j = l, loop over k
                for k = 1:Nq
                    l = j;
                    q = k + Nq*(l-1);
                    C(r,q) = C(r,q) + ...
                        Jy*a(xp(i))*b(yp(j))*w(j)*w(i)*Dh(i,k);                    
                end
                % delta_ik -> i = k, loop over l
                for l = 1:Nq
                    k = i;
                    q = k + Nq*(l-1);
                    C(r,q) = C(r,q) + ...
                        Jx*c(xp(i))*d(yp(j))*w(i)*w(j)*Dh(j,l);
                end                
            end
        end                    
        elem_offset = ((ky-1)*Kx + (kx-1))*Nq2;
        local_inds = 1:Nq2;
        inds = galnums(elem_offset + local_inds);
        Kg(inds,inds) = Kg(inds,inds) + K;
        Mg(inds,inds) = Mg(inds,inds) + M;
        Cg(inds,inds) = Cg(inds,inds) + C;
        %fg(inds) = fg(inds) + f;
        disp(['kx = ',num2str(kx), ', ky = ', num2str(ky)])
    end
end
Nqkx = kx*(Nq-1) + 1; % num global dofs along x line
Nqky = ky*(Nq-1) + 1; % num global dofs along y line

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

% forcing 
F = Fx(X).*Fy(Y);
f = reshape(F,Nqkx*Nqky,1);
fg = Mg*f;

% BCs
bottom = 1:Nqkx;
left = (0:Nqky-1)*Nqkx + 1;
right = (1:Nqky)*Nqkx;
top = (numdofs-Nqkx):numdofs;
bcInds = unique([bottom, left, right, top]);
% bcInds = unique([bottom, left]);

% % Dirichlet BCs
Cg(bcInds,:) = zeros(size(Cg(bcInds,:)));
Cg(:,bcInds) = zeros(size(Cg(:,bcInds)));
Cg(bcInds,bcInds) = eye(length(bcInds));

Mg(bcInds,:) = zeros(size(Mg(bcInds,:)));
Mg(:,bcInds) = zeros(size(Mg(:,bcInds)));
Mg(bcInds,bcInds) = eye(length(bcInds));

Kg(bcInds,:) = zeros(size(Kg(bcInds,:)));
Kg(:,bcInds) = zeros(size(Kg(:,bcInds)));
Kg(bcInds,bcInds) = eye(length(bcInds));

fg(bcInds) = 0;
u0(bcInds) = 0;

% fg_t =  fg + (1/dt)*Mg*u0;

% %fg_t = fg + (1/dt)*Mg*u0 - Kg*u0; %explicit
% 
% figure
% pcolor(xx,yy,U0)
% ax = axis;
% cax = caxis;
% view(2)
% pause
% Nsteps = 1/dt;
% for i = 1:Nsteps
%     % implicit
%     ug = ((1/dt)*Mg+eps*Kg+Cg)\fg_t;
%     fg_t = fg + (1/dt)*Mg*ug; % next timestep
%     
%     % explicit
%     %ug = dt*(1./diag(Mg)).*fg_t;
%     %fg_t = fg + (1/dt)*Mg*ug - Kg*ug;
%     
%     pcolor(xx,yy,reshape(ug,Nqkx,Nqky))    
%     caxis(cax)
%     colorbar
%     title(['Time = ',num2str(i*dt)])
%     axis(ax)
%     view(2)
%     drawnow
% end

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