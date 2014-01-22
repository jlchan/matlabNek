% constants
one = @(x) ones(size(x));
zero = @(x) zeros(size(x));
for N = 1:4
    for Ki = 1:3
        % N = 1; % order
        
        K = 2^(Ki+1);
        
        dt = .01;
        eps = .01;
        c = 0;
        
        beta = {one,one,zero,zero};
        %f = {one, one};
        
        % manufactured solns
        ue = {@(x) sin(pi*x), @(y) sin(pi*y)};
        f = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);
        %ue = {@(x) x.*(1-x), @(y) y.*(1-y)};
        %f = @(x,y) -2*x.*(x-1) - 2*y.*(y-1);
        
        % clean this part up - consolidate
        kx = K;ky = K;Nq = N+1;
        Nqkx = kx*(Nq-1) + 1; Nqky = ky*(Nq-1) + 1; % num global dofs along line
        [Mg, Kg, Cg, bcInds, fvec, u0] = assemble(N,K,beta,f); % get global SEM matrices
        [X Y] = get_physical_points(N,kx,ky);
        
        % manufactured solution checking mass matrix
        fp = Mg\fvec;
        fex = reshape(f(X,Y),Nqkx*Nqky,1);
        fex(bcInds) = zeros(size(bcInds));
        norm(fp-fex)
        
        % manufactured solution checking diffusion matrix
        u = Kg\fvec;
        uex = reshape(ue{1}(X).*ue{2}(Y),Nqkx*Nqky,1);
        err(N,Ki) = (u-uex)'*Mg*(u-uex);
        hN(N,Ki) = Nqkx*Nqky;
    end
end