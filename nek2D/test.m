clear
N = 1; % order
Nq = N+1;
Nq2 = Nq*Nq;

Kx = 2;
Ky = 2;

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
                starts(gid) = starts(gid) + 1; % should be global to local map
            end
        end
    end
end
numdofs = gid;

A = ones(Nq2);
K = zeros(numdofs);
for kx = 1:Kx
    for ky = 1:Ky
        elem_offset = ((ky-1)*Kx + (kx-1))*Nq2;
        local_inds = 1:Nq2;
        inds = galnums(elem_offset + local_inds);
        K(inds,inds) = K(inds,inds) + A;
    end
end
spy(K)


starts = cumsum(starts);
for i = 1:length(galnums)
    gid = galnums(i)
    indices(starts(gid) + count(gid)) = i % global to local 
    count(gid) = count(gid) + 1
end