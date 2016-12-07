function L = grassOpt_cg(g,dg,L0,opts)

% grassOpt_cg : solves an optimization problem on the 
% Grassman manifold using a conjugate gradient method.
%
% This function solves the problem:
%
% min g(L)
% s.t. L'L = I (not necessarily square)
%
% where the function g() is invariant with respect to
% an orthonormal transformation of L, i.e. when 
% g(L) = g(LR) for any R'R = I, R square.
%
% Usage : L = grassOpt_cg(g,dg,L0,opts)
%
% Inputs:
%
% g     : function handle for evaluating g(L)
% dg    : function handle for evaluating \nabla g(L)  
% L0    : initial condition for L.  This must be
%         a matrix satisfying the constraint L'L=I.
% opts  : options structure with the following fields:
%         'maxIter', 'relTol' and 'maxStep'.
%
% Outputs:
%
% L     : solution for the above problem
%
% This function is an implementation of the algorithm described
% in the following publications: 
%
% Edelman et. al. 'The geometry of algorithms with orthogonality 
% constraints' SIAM J. Matrix Anal. Appl. 20(2) 303-353).

% Author : P. Goulart - 20 June 2012
%


%intialize convergence test variables
iterCount = 0;
relError  = 1;
valueOld  = g(L0); 

%initialize algorithm parameters
Yk = L0;
dgdY = dg(Yk); 
Gk   = dgdY - Yk*(Yk'*dgdY);
Hk   = -Gk;  

%print the reporting header
printInfo();

%get an initial step-size
tstep = opts.maxStep;

while( abs(relError) > opts.relTol && iterCount < opts.maxIter)
    
    %update counter
    iterCount = iterCount + 1;
            
    %minimize g(L_k(t)) over t
    [U,S,V] = svd(Hk,0);
    
    %do a line search to minimize g(t).  
    k = pi/max(diag(S));  %normalization for geodesic search
    fline = @(t)(g(Yk*V*dcos(S*(k*t))*V' + U*dsin(S*(k*t))*V'));
    [tstep,valueNew] = linesearch(fline,0,tstep,opts.lineSearch,valueOld);
   
    %update L along a geodesic
    Ykp1 = Yk*V*dcos(S.*(k*tstep))*V' + U*dsin(S.*(k*tstep))*V';
        
    %renormalize, just in case...
    Ykp1 = orth(Ykp1);
    
    %find gradient for new L.
    dgdY = dg(Ykp1);
    
    %update G
    Gkp1 = dgdY - Ykp1*(Ykp1'*dgdY);
    
    %parallel transport H and G
    tHk = (-Yk*V*dsin(S.*(tstep*k))+U*dcos(S.*(tstep*k)))*S*V';
    tGk = Gk - (Yk*V*dsin(S.*(tstep*k))+U*(eye(size(S)) - dcos(S.*(tstep*k))))*(U'*Gk);
    
    %new search direction
    normGk2 = norm(Gk,'fro')^2;
    gammak  = sum(sum(Gkp1.*(Gkp1-tGk))) ./ normGk2;
    Hkp1    = -Gkp1 + gammak.*tHk; 
            
    %calculate relative improvement
    relError = (valueOld-valueNew)/abs(valueOld);
    valueOld = valueNew; %for next pass
    
    %print reporting information
    printInfo(iterCount,valueNew,relError);
    
    %If the search direction and gradient are too far
    %apart (arbitarily defined...), then restart
    if(trace(Gkp1'*Hkp1) / trace(Gkp1'*Gkp1) > -0.5)
        Hkp1 = -Gkp1;
        tstep = opts.maxStep;
        fprintf('Reinitializing conjugate gradient method\n');
    end
        
    if(relError <= 0 || tstep == 0)
        %not improving.  Force halt.
        relError = 0;
        L = Yk; %from the previous iteration
    else
        %update L, G, H iterates
        Yk = Ykp1;
        Gk = Gkp1;
        Hk = Hkp1;
        %update the output matrix
        L  = Yk;
    end
    
    %NB: p(n-p) can be really huge, so don't 
    %reset Hk as in Edelman paper.
    
end


if(iterCount == opts.maxIter)
    warning('Didn''t converge')
end



%----------------------------------------------------
%----------------------------------------------------

function printInfo(iterCount,val,relerr)

%print reporting information

if(nargin < 1)
    %print the header
    fprintf(1,'\n\nStarting solver: %s\n\n',mfilename);
    fprintf(1,'iterate  |  objective value  |  relative improvement\n');
    fprintf(1,'----------------------------------------------------\n');
else
    fprintf(1,'%4i     |  %0.6e    |  %0.6e \n',iterCount,val,relerr);
end


   
%----------------------------------------------------
%----------------------------------------------------

function D = dsin(S)

%Computes a diagonal matrix by taking 
%sines of the diagonal elements of the
%input.  This is 'sin' in the sense of 
%Edelman (2.65)

D = diag(sin(diag(S)));





%----------------------------------------------------
%----------------------------------------------------

function D = dcos(S)

%Computes a diagonal matrix by taking 
%sines of the diagonal elements of the
%input.  This is 'cos' in the sense of 
%Edelman (2.65)

D = diag(cos(diag(S)));



       








