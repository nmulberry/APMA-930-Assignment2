%% Unsteady Navier Stokes

% Solve the driven cavity problem for unsteady Navier Stokes flow in a wedge
% using the streamfunction/vorticity formulation.  
% Uses SBDF time stepping scheme & central differences
% Calls functions SystemMatComp3 and ConstructRhs3


close all;  clear; clc;

global Re

% =======================================================================%
%                          Setup                                         %
% =======================================================================%

% Parameters
U     = -1;
Rmax  = 1;
alpha = pi/2;



% Reynold's number. 
Re = 20;

% Set up finite difference grid
    M = 50; dr = Rmax/(M-1);
    N = 50; dth = alpha/(N-1);
    [rg, thg] = meshgrid(0: dr :Rmax, ...
                        alpha: -dth: 0);  
                    
    T  = 3; dt = 0.1;
    
% Unknowns and numbering
    numUn = M*N;
    nP = reshape(1:numUn, size(rg));
    nO = reshape(numUn+1:2*numUn, size(rg));
    numUn = 2*numUn;

% =======================================================================%
%                          Solver                                        %
% =======================================================================%
 
% Initialize fluid at rest
psivort0   = zeros(numUn, 1);
psivort   = zeros(numUn, 1);
tic;

% System Matrix 
PsiOmSys = SystemMatComp3(numUn, nP, nO, M, N, alpha, dr, dth, dt);

% Evolve in time 

for time = 0:dt:T
    
    %-------------- Build rhs---------------------------------------%
    
    rhs  = ConstructRhs3(numUn, nO, M, N, Rmax, dr, dth, dt, U, psivort, psivort0);
    
    psivort0 = psivort;
     
        
    
    %--------------- Solve -----------------------------------------%
    psivort = PsiOmSys \ rhs ;  
       
           
    
    % -------------- Plot ------------------------------------------%
    psi   = reshape(psivort(1:numUn/2), size(rg));
    omega = reshape(psivort(numUn/2+1:numUn), size(rg));
    
    clf(); 
    subplot(1, 2, 1)
        contour(rg.*cos(thg), rg.*sin(thg), psi, [0 0],'k','LineWidth',2); 
        hold on;
        contour(rg.*cos(thg), rg.*sin(thg), psi, 40);         
        hold on;
        colorbar;
        shading flat;  colormap(jet);  
        title('Streamfunction')
        axis([0 1 0 1])
        axis square
        
    subplot(1, 2, 2)
        c = linspace(-20, 20, 40);
        k = contour(rg.*cos(thg), rg.*sin(thg), omega, c);
        colorbar;
        shading flat;  colormap(jet);  
        title('Vorticity')
        axis([0 1 0 1])
        axis square
        
     ha = axes('Position',[0 0  1 1],'Xlim',[0 1],'Ylim',[0 ...
            1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
     
     
     text(0.5, 1, ([' Re = ' num2str(Re) ', t = ' num2str(time)]), 'HorizontalAlignment' ...
            ,'center','VerticalAlignment', 'top');
        
        drawnow;
        
        
end  

computationTime = toc;

