function L = lmlOpt_altProject(A,B,L0,opts)

% grassOpt_alternate : finds an approximate solution to the low rank
% matrix approximation problem in the form:
%
% min \|A - XB\|^2
% s.t. X = LML', 
%
% where L is an orthonormal matrix with k columns, and M is 
% a square matrix of size k x k.
%
% Usage : grassOpt_alternate(A,B,L0,opts)
%
% Inputs:
%
% A,B    : data matrices of size p x n
% L0     : initial condition for L.  
% opts   : options structure with the following fields:
%         'maxIter', 'relTol' and 'maxStep'.
%
% Outputs:
%
% L      : optimal solution basis for the problem above,  
%          with L \in R^{p x k} 
%
% The input data matrices (A,B) must both be size
% p x n, and full column rank, with k \le n \le p
%
%  

% Author : P. Goulart - 20 June 2012
%
% This function is an implementation of the alternating 
% projection method described in Algorithm 1 of:
%
% Goulart, Wynn & Pearson, 'Optimal mode decomposition for 
% high-dimensional systems. In 51st IEEE Conference on Decision 
% and Control. Maui, Hawaii, Dec. 2012.  Available at 
% http:\\control.ee.ethz.ch\~goularpa\.


%initialize convergence test variables
iterCount = 0;
relError  = 1;
valueOld  = realmax; 

%initialize algorithm parameters
L = L0;
Llast = L;
k = size(L,2);
%C = orth(A'*L);
C = eye(size(A,2),size(L,2));

%print the reporting header
printInfo();

%norm squared of A.  Calculated so that
%reporting of values will be consistent
%with the overall objective function
nA = trace(A'*A);


while( abs(relError) > opts.relTol && iterCount < opts.maxIter)
    
    %update counter
    iterCount = iterCount + 1;
            
    %Find optimal L' ignoring constraints on C
    L = solve_inner((A*C)',k);
    
    %compute the image of B'*L
    C = orth(B'*L);
    
    %calculate relative improvement
    valueNew = sqrt(nA - getNorm(L,A,C)^2);
    relError = (valueOld - valueNew)/abs(valueOld);
    valueOld = valueNew; %for next pass
    
    %print reporting information
    printInfo(iterCount,valueNew,relError);
    
    if(relError < 0)
        %going backwards.  Force halt.
        relError = 0;
        L = Llast; %from the previous iteration
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

function X = solve_inner(D,k)

% solver for the unconstrained step
%

%Compute a basis for D'
V = orth(D');

%X coincides with the first k right singular vectors
X = V(:,1:k);




%----------------------------------------------------
%----------------------------------------------------

function n = getNorm(L,A,C)

n = norm((L'*A)*C,'fro');







