% test_agmg

dt = .1;
eps = .01;
Nvec = 2:8;
Kvec = 2:8;
Iters = zeros(length(Kvec),length(Nvec));

for i = 1:length(Kvec)
    for j = 1:length(Nvec)
        K = Kvec(i);N = Nvec(j);       
        [Mg, Kg, Cg, bcInds, u0] = assemble(N,K); % get global SEM matrices
        % full matrix
        A = (1/dt)*Mg + eps*Kg + Cg;
        
        % time-rhs
        Rhs = (1/dt)*Mg*u0;
        
        % bcs
        A(bcInds,:) = zeros(size(A(bcInds,:)));
        A(:,bcInds) = zeros(size(A(:,bcInds)));
        A(bcInds,bcInds) = eye(length(bcInds));
        Rhs(bcInds) = 0;
        
        levels = agmg_setup(A);
        
        tol = 1e-7;
        maxiter = 100;
        
        % solve
        [x flag relres iter resvec] = agmg_solve(levels, Rhs, maxiter, tol);
        Iters(i,j) = iter;
    end    
end

