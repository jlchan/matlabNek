      function[Ah,Bh,Ch,Dh,z,w] =  SEMhat(N)
%
%                                                 ^
%     Compute the single element 1D SEM Stiffness Mass, and Convection
%     matrices, as well as the points and weights for a polynomial
%     of degree N
%

      [z,w] = zwgll(N);

      Bh    = diag(w); % diag of weights (mass matrix)
      Dh    = Dhat(z); % derivative matrix

      Ah    = Dh'*Bh*Dh; % stiffness matrix
      Ch    = Bh*Dh; % convection matrix

