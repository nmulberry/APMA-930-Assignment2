function r = ConstructRhs2(numUn, nO, M, N, R, dr, dth, U, psiOm, Re)
   

r = zeros(numUn, 1);


psi   = psiOm(1:numUn/2);
omega = psiOm(numUn/2+1:end);

rg = (0:dr:R)';
rg = 1./rg;
rg(1) = 0;

e = ones(M,1);
idr  = speye(N); idr(1,1) = 0; idr(N,N) = 0;
idth = speye(M); idth(1,1) = 0; idth(M,M) = 0;

% R-direction
Dr_0  = spdiags([-rg, rg], [-1,1] , N, N);

% Theta-direction 
Dth_0 = spdiags([-e, e], [-1,1] , M, M);


% Mesh

Upsi = kron(Dth_0, idr)*psi;
Vpsi = kron(idth,Dr_0)*psi;

Uom = kron(idth, Dr_0)*omega;
Vom = kron(Dth_0,idr)*omega;


J = Uom.*Upsi - Vom.*Vpsi;
J = (Re/(4*dr*dth))*J;

r(numUn/2+1:end) = full(J);

%  Boundary conditions at right: omega given
for jrow = 2: N-1
    ijO = nO(jrow, M);
    r(ijO) = U/R + 3*U/dr;
end

end