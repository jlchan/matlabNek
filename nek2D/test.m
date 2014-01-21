N = 3; % order
K = 8;
CFL = 10;
dt = CFL/(K*N^2);
eps = .01;
c = 1;
beta = {@(x) ones(size(x)),@(x) ones(size(x)),@(x) zeros(size(x)),@(x) zeros(size(x))};

% clean this part up - consolidate
kx = K;ky = K;
Nq = N+1;
Nqkx = kx*(Nq-1) + 1; Nqky = ky*(Nq-1) + 1; % num global dofs along line

[Mg, Kg, Cg, bcInds, f, u0] = assemble(N,K,beta); % get global SEM matrices

[X Y] = get_physical_points(N,kx,ky);

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

middle = round(Nqkx*Nqky/2);
e = zeros(Nqkx*Nqky,1);
e(middle) = .01;
G = A\e;
surf(X,Y,reshape(G,Nqkx,Nqky));

break

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
