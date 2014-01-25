% N = 10; % order
% K = 16; % num elems

function [Mg, Kg, Cg, bcInds, fg, u0, galnums] = assemble(N,K,beta,ff)

if (nargin<2)
    N = 4;
    K = 4;
end
if (nargin<3)
    % circulatory flow
    beta = {@(x) ones(size(x)), @(y) 5*(y-.5), @(x) -5*(x-.5), @(y) ones(size(y))};    
end
if (nargin<4)   
   % some nonconst nonzero forcing
   ff = @(x,y) .5*exp(x).*y.^2;
end
a = beta{1};b = beta{2};c=beta{3};d=beta{4};

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
% a = @(x) ones(size(x));
% b = @(y) 5*(y-.5);
% c = @(x) -5*(x-.5);
% d = @(y) ones(size(y));


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
%Mg = sparse(numdofs,numdofs);
%Kg = sparse(numdofs,numdofs);
%Cg = sparse(numdofs,numdofs);
% fg = zeros(numdofs,1);
rowInds = [];colInds = [];
cvec = [];
for kx = 1:Kx
    for ky = 1:Ky        
        Jx = dx/2; Jy = dy/2;
        % assemble element matrices
        C = zeros(Nq2,Nq2);
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
                        Jy*c(xp(i))*d(yp(j))*w(j)*w(i)*Dh(i,k);                    
                end
                % delta_ik -> i = k, loop over l
                for l = 1:Nq
                    k = i;
                    q = k + Nq*(l-1);
                    C(r,q) = C(r,q) + ...
                        Jx*a(xp(i))*b(yp(j))*w(i)*w(j)*Dh(j,l);
                end                
            end
        end                    
        elem_offset = ((ky-1)*Kx + (kx-1))*Nq2;
        local_inds = 1:Nq2;
        inds = galnums(elem_offset + local_inds);
        [I J] = meshgrid(inds);
        rowInds = [rowInds; I(:)];
        colInds = [colInds; J(:)];
        cvec = [cvec; C(:)];
        %Kg(inds,inds) = Kg(inds,inds) + K;
        %Mg(inds,inds) = Mg(inds,inds) + M;
        %Cg(inds,inds) = Cg(inds,inds) + C;
        %fg(inds) = fg(inds) + f;        
    end
    disp(['kx = ',num2str(kx), ', ky = ', num2str(ky)])
end
Cg = sparse(colInds,rowInds,cvec,numdofs,numdofs,numdofs*Nq2*4);
mvec = repmat(M(:),kx*ky,1);
Mg = sparse(colInds,rowInds,mvec,numdofs,numdofs,numdofs*Nq2*4);
kvec = repmat(K(:),kx*ky,1);
Kg = sparse(colInds,rowInds,kvec,numdofs,numdofs,numdofs*Nq2*4);

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
%F = Fx(X).*Fy(Y); % interpolate and integrate
F = ff(X,Y);
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
