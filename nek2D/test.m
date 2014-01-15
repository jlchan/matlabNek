clear
N = 1; % order
Nq = N+1;
Nq2 = Nq*Nq;

Kx = 4;
Ky = 4;

xv = linspace(0,1,Kx+1);
yv = linspace(0,1,Ky+1);
[Xv Yv] = meshgrid(xv,yv);

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
%                 starts(gid) = starts(gid) + 1; % should be global to local map
            end
        end
    end
end
numdofs = gid;

% make 2D element matrices
dt = .1;
eps = 1.0;
bb = 1.0;
a = bb; % define beta = (ab,cd)
b = bb;
c = bb;
d = bb;
M = zeros(Nq2);
K = zeros(Nq2);
C = zeros(Nq2);
[Ah,Bh,Ch,Dh,z,w] = SEMhat(N); 
for i = 1:Nq
    for j = 1:Nq
        r = i + Nq*(j-1);
        M(r,r) = w(i)*w(j);        
        
        % delta_jl -> j = l, loop over k
        for k = 1:Nq
            l = j;
            q = k + Nq*(l-1);
            K(r,q) = K(r,q) + ...
                Ah(i,k)*w(j); 
            C(r,q) = C(r,q) + ...
                a*b*w(j)*w(k)*Dh(k,i);
                
        end        
        % delta_ik -> i = k, loop over l
        for l = 1:Nq
            k = i;
            q = k + Nq*(l-1);
            K(r,q) = K(r,q) + ...
                Ah(j,l)*w(i);             
            C(r,q) = C(r,q) + ...
                c*d*w(i)*w(l)*Dh(l,j);            
        end
        
    end
end
% for i = 1:Nq
%     for j = 1:Nq
%         r = i + Nq*(j-1);
%         for k = 1:Nq
%             for l = 1:Nq
%                 q = k + Nq*(l-1);
%                 K(r,q) = K(r,q) + ...
%                     Ah(i,k)*w(j)*(j==l) + ...
%                     Ah(j,l)*w(i)*(i==k);
%             end
%         end
%     end
% end

A = M + eps*K + (1/dt)*C;
Kglob = zeros(numdofs);
for kx = 1:Kx
    for ky = 1:Ky
        elem_offset = ((ky-1)*Kx + (kx-1))*Nq2;
        local_inds = 1:Nq2;
        inds = galnums(elem_offset + local_inds);
        Kglob(inds,inds) = Kglob(inds,inds) + A;
    end
end
spy(Kglob)
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