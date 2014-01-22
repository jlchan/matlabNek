% constants
one = @(x) ones(size(x));
zero = @(x) zeros(size(x));

N = 4; % order
K = 16;
% CFL = 1; dt = CFL/(K*N^2);
dt = .01;
eps = .01;
c = 1;

%beta = {one,one,zero,zero};
beta = {one,one,one, one};
%f = {one, one};

% manufactured solns
ue = {@(x) sin(pi*x), @(y) sin(pi*y)};
f = @(x,y) zeros(size(x));
%ue = {@(x) x.*(1-x), @(y) y.*(1-y)};
%f = @(x,y) -2*x.*(x-1) - 2*y.*(y-1);

% clean this part up - consolidate
kx = K;ky = K;Nq = N+1;Nqkx = kx*(Nq-1) + 1; Nqky = ky*(Nq-1) + 1; % num global dofs along line
[Mg, Kg, Cg, bcInds, fvec, u0] = assemble(N,K,beta,f); % get global SEM matrices
[X Y] = get_physical_points(N,kx,ky);

% initial condtion
% initial condtion
x0 =0.6; y0=0.3; delta = 0.10; R=(X-x0).^2+(Y-y0).^2;
U0 = exp(-((R./(delta^2)).^1)).*X.*(1-X).*Y.*(1-Y); % pulse * bubble
u0 = reshape(U0,Nqkx*Nqky,1); 

fg_t =  fvec + (1/dt)*Mg*u0;
A = (1/dt)*Mg + eps*Kg + c*Cg;

% time-dependent
figure
pcolor(X,Y,U0)
ax = axis;
cax = caxis;
pause
Nsteps = 1/dt;
[L U] = lu(A);
for i = 1:Nsteps
    % implicit
    ug = U\(L\fg_t);
    fg_t = fvec + (1/dt)*Mg*ug; % next timestep
    fg_t(bcInds) = 0;
    
    % explicit
    %ug = dt*(1./diag(Mg)).*fg_t;
    %fg_t = fg + (1/dt)*Mg*ug - Kg*ug;
    
    if (mod(i,Nsteps/10)==0)
        surf(X,Y,reshape(ug,Nqkx,Nqky))
%         caxis(cax)
        colorbar
        title(['Time = ',num2str(i*dt)])
%         axis(ax)        
        drawnow
    end
end
