function [PsiOm] = SystemMatComp3(numUn, nP, nO, M, N, alpha, dr, dth, dt)

global Re;
% does sparse storage

    i = zeros(1, 11*M*N);
    j = i;
    PsiOm = i;
    
%
%  Boundary conditions at origin 
    k = 1;
    for jrow = 1: N
        ij = nP(jrow, 1);
        ijp1 = nP(jrow, 2);
        ijO = nO(jrow, 1);
        i(k:k+2) = ij;
        j(k:k+2) = [ij ijp1 ijO];
        PsiOm(k) = (N-2);
        PsiOm(k+1) = -1;
        PsiOm(k+2) = -alpha*(dr/2)^2/dth;
        k = k+3;        
        ijOp1 = nO(jrow, 2);
        i(k:k+1) = ijO;
        j(k:k+1) = [ijO ijOp1];
        PsiOm(k) = 1/dr;
        PsiOm(k+1) = -1/dr;
        k = k+2;
    end
    for icol = 2: M-1
        r_i = (icol-1)*dr;
%
%  Boundary condition at top - Dirichlet 
        ij = nP(1, icol);
        i(k) = ij;
        j(k) = ij;
        PsiOm(k) = 1;
        k = k+1;
        ijO = nO(1, icol);
        ijm1 = nP(2, icol);
        ijm2 = nP(3, icol);
        i(k: k+2) = ijO;
        j(k: k+2) = [ijO ijm1 ijm2];
        PsiOm(k) = 1;
        PsiOm(k+1) = 4/(r_i*dth)^2;
        PsiOm(k+2) = -0.5/(r_i*dth)^2;
        k = k+3;
%
%  Interior points - 5 point stencil
        for jrow = 2: N-1
            ij = nP(jrow, icol);
            ip1j = nP(jrow, icol+1);
            im1j = nP(jrow, icol-1);
            ijp1 = nP(jrow-1, icol);
            ijm1 = nP(jrow+1, icol);
            ijO = nO(jrow, icol);
            i(k: k+5) = ij;
            j(k: k+5) = [ij im1j ip1j ijm1 ijp1 ijO];
            PsiOm(k+1) = -(r_i - 0.5*dr)/(r_i*dr^2);
            PsiOm(k+2) = -(r_i + 0.5*dr)/(r_i*dr^2);
            PsiOm(k+3) = -1/(r_i*dth)^2;
            PsiOm(k+4) = -1/(r_i*dth)^2;
            PsiOm(k) = 2/dr^2 + 2/(r_i*dth)^2;
            PsiOm(k+5) = -1;   % lapl(psi) = -Omega;
            k = k+6;
            
            alpha = 2*dt/(3*Re);
            ij = ijO;
            ip1j = nO(jrow, icol+1);
            im1j = nO(jrow, icol-1);
            ijp1 = nO(jrow-1, icol);
            ijm1 = nO(jrow+1, icol);
            i(k: k+4) = ij;
            j(k: k+4) = [ij im1j ip1j ijm1 ijp1];
            PsiOm(k+1) = 1-alpha*(r_i - 0.5*dr)/(r_i*dr^2);
            PsiOm(k+2) = 1-alpha*(r_i + 0.5*dr)/(r_i*dr^2);
            PsiOm(k+3) = 1-alpha/(r_i*dth)^2;
            PsiOm(k+4) = 1-alpha/(r_i*dth)^2;
            PsiOm(k) = 1+alpha*(2/dr^2 + 2/(r_i*dth)^2);
            k = k + 5;
        end
%
%  Boundary condition at bottom - Dirichlet
        ij = nP(N, icol);
        i(k) = ij;
        j(k) = ij;
        PsiOm(k) = 1;
        k = k+1;
        ijO = nO(N, icol);
        ijp1 = nP(N-1, icol);
        ijp2 = nP(N-2, icol);
        i(k:k+2) = ijO;
        j(k:k+2) = [ijO ijp1 ijp2];
        PsiOm(k) = 1;
        PsiOm(k+1) = 4/(r_i*dth)^2;
        PsiOm(k+2) = -0.5/(r_i*dth)^2;
        k = k+3;
    end
%
%  Boundary conditions at right - Dirichlet
    for jrow = 1: N
        ij = nP(jrow, M);
        i(k) = ij;
        j(k) = ij;
        PsiOm(k) = 1;
        k = k+1;
        ijO = nO(jrow, M);
        im1j = nP(jrow, M-1);
        im2j = nP(jrow, M-2);
        i(k:k+2) = ijO;
        j(k:k+2) = [ijO im1j im2j];
        PsiOm(k) = 1;
        PsiOm(k+1) = 4/dr^2;
        PsiOm(k+2) = -0.5/dr^2;
        k = k+3;
    end  
    k = k - 1;
    
    i = i(1:k);
    j = j(1:k);
    PsiOm = PsiOm(1:k);
    PsiOm = sparse(i, j, PsiOm, numUn, numUn);
end