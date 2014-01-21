N = 4; % order
K = 16;
CFL = 10;
dt = Inf*CFL/(K*N^2);
eps = .00001;
c = 1;
beta = {@(x) ones(size(x)),@(x) ones(size(x)),@(x) zeros(size(x)),@(x) zeros(size(x))};

% clean this part up - consolidate
kx = K; ky = K; Nq = N+1; Nqkx = kx*(Nq-1) + 1; Nqky = ky*(Nq-1) + 1; % num global dofs along line

[Mg, Kg, Cg, bcInds, f, u0] = assemble(N,K,beta); % get global SEM matrices

[X Y] = get_physical_points(N,kx,ky);

A = (1/dt)*Mg + eps*Kg + c*Cg;
% bcs
A(bcInds,:) = zeros(size(A(bcInds,:)));
A(:,bcInds) = zeros(size(A(:,bcInds)));
A(bcInds,bcInds) = eye(length(bcInds));

middle = round(Nqkx*Nqky/2);
e = zeros(Nqkx*Nqky,1);
e(middle) = .01;
G = A\e;
figure
surf(X,Y,reshape(G,Nqkx,Nqky));
shading interp

% agmg inverse
levels = agmg_setup(A);
g_agmg = precond_agmg(levels, e);

figure
surf(X,Y,reshape(g_agmg,Nqkx,Nqky));
shading interp