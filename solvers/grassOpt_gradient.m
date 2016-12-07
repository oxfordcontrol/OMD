function L = grassOpt_gradient(g,dg,L0,opts)

% grassOpt_gradient : solves an optimization problem on the 
% Grassman manifold using a gradient method.
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
L     = L0;
Llast = L0;

%print the reporting header
printInfo();

%get an initial step-size
tstep = opts.maxStep;

while( abs(relError) > opts.relTol && iterCount < opts.maxIter)
    
    %update counter
    iterCount = iterCount + 1;
            
    %find the component of the gradient tangent to 
    %the Grassman manifold
    %following the method of \S2.5.3 in Edelman
    dgdL   = dg(L);
    G      = dgdL - L*(L'*dgdL);
    
    %search direction is negative gradient
    D = -G;
    
    %compact SVD decomp for D
    [U,S,V] = svd(D,0);
    
    %do a line search to minimize g(t).  
    k = pi/max(diag(S));  %normalization for geodesic search
    fline = @(t)(g(L*V*dcos(S*(k*t))*V' + U*dsin(S*(k*t))*V'));
    [tstep,valueNew] = linesearch(fline,0,tstep,opts.lineSearch,valueOld);
    
    %update L along a geodesic
    L = L*V*dcos(S.*(k*tstep))*V' + U*dsin(S.*(k*tstep))*V';
    
    %renormalize, just in case...
    L = orth(L);
            
    %calculate relative improvement
    relError = (valueOld-valueNew)/abs(valueOld);
    valueOld = valueNew; %for next pass
    
    %print reporting information
    printInfo(iterCount,valueNew,relError);
        
    if(relError <= 0 || tstep == 0)
        %not improving.  Force stop.
        relError = 0;
        L = Llast; %previous iterate
    else
        Llast = L;
    end

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









