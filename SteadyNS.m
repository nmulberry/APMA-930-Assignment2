%% Steady Navier Stokes

% Solve the driven cavity problem for Steady Navier Stokes flow in a wedge
% using the streamfunction/vorticity formulation. Calls functions
% SystemMatComp2 and ConstructRhs2


close all;  clear; clc;


% Problem parameters:
    U = -1;
    Rmax = 1;
    alpha = pi/2;
    
    % Works up to ~24 ??
    Re = 24;

% Set up finite difference grid
    M = 50; dr = Rmax/(M-1);
    N = 50; dth = alpha/(N-1);
    [rg, thg] = meshgrid(0: dr :Rmax, ...
                         alpha: -dth: 0);                     

% Unknowns and numbering
    numUn = M*N;
    nP = reshape(1:numUn, size(rg));
    nO = reshape(numUn+1:2*numUn, size(rg));
    numUn = 2*numUn;

%% Iterations

itmax = 10000;
tol = 1.d-9;
   

% Initialize

iter = 0;
psivort   = zeros(numUn, 1);
tic;

% System Matrix 
PsiOmSys = SystemMatComp2(numUn, nP, nO, M, N, alpha, dr, dth);

for reTest = Re   % continuation goes here, doesn't seem to be helping?
    normPsi = 1;
    normOm = 1;
    
    while normPsi > tol || normOm > tol && iter < itmax

    % Build rhs

        rhs  = ConstructRhs2(numUn, nO, M, N, Rmax, dr, dth, U, psivort, reTest);

        psivort0 = psivort;

    % Solve
        psivort = PsiOmSys \ rhs ;  

        normPsi = norm(psivort(1:numUn/2) - psivort0(1:numUn/2));
        normOm  = norm(psivort(numUn/2+1:numUn) - psivort0(numUn/2+1:numUn));
        disp(['Iteration: ', num2str(iter), ...
            ' Residual Psi = ', num2str(normPsi), ...
            ' Residual Omega = ', num2str(normOm), ...
            ' Re = ', num2str(reTest)])
        
        iter = iter + 1;

    end  
end
%%
t = toc;

% Plot
psi   = reshape(psivort(1:numUn/2), size(rg));
omega = reshape(psivort(numUn/2+1:numUn), size(rg));
    
    figure()
    subplot(1, 2, 1)
        contour(rg.*cos(thg), rg.*sin(thg), psi, [0 0],'k','LineWidth',2); 
        hold on;
        contour(rg.*cos(thg), rg.*sin(thg), psi, 40);         
%         c = linspace(-5d-7, 0, 40);
%         contour(rg.*cos(thg), rg.*sin(thg), psi, c); 
        hold on
        colorbar;
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Streamfunction')
        axis([0 1 0 1])
        axis square
    subplot(1, 2, 2)
        c = linspace(-20, 20, 40);
        contour(rg.*cos(thg), rg.*sin(thg), omega, c);
        colorbar;
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Vorticity')
        axis([0 1 0 1])
        axis square