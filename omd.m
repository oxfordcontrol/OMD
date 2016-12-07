function [L,M] = omd(A,B,k,L0,method,userOpts)

%OMD : solves an optimal low rank approximation problem 
% in the form:
%
% min \|A - XB\|^2
% s.t. X = LML', 
%
% where L is an orthonormal matrix with k columns, and M is 
% a square matrix of size k x k.
%
% Usage : [L,M] = omd(A,B,k,L0,method)
%
% Inputs:
%
% A,B    : data matrices of size p x n
% k      : integer specifying the desired output size.  
% L0     : initial condition for L.  If L0 is not specified,
%          then the initial iterate will be based on the 
%          first k singular vectors of the matrix [A B]
%
% method : Specifies which optimization method to use.
%          It must be one of the following:
%          
%          'alternating' : Uses an approximate method of
%                          alternating directions.  Good for
%                          obtaining a rough solution quickly
%
%          'gradient'    : Uses the gradient descent method
%
%          'conjgrad'    : Uses the conjugate gradient method
%
%          'hybrid'      : Uses the alternating method first
%                          to get a rough solution, followed
%                          by the conjugate gradient method
%                          to achieve local optimality.  This
%                          is usually the best way to get a 
%                          high-precision solution quickly.
%
%          'dmd'         : Returns the dynamic mode decomposition
%
%          Note that if the 'dmd' method is chosen, then no
%          optimization is performed over the input basis L0. 
%          If L0 is not specified, then function returns the 
%          standard dynamic mode decomposition.
%                   
%
% Outputs:
%
% L      : optimal solution basis for the problem above,  
%          with L \in R^{p x k} 
% M      : optimal solution matrix for the problem above,
%          with M \in R^{k x k}
%
% The input data matrices (A,B) must both be size
% p x n, and full column rank, with k \le n \le p
%
%  

% Author : P. Goulart - 20 June 2012
%
% This function is an implementation of the algorithms described
% in the following publications: 
%
% Goulart, Wynn & Pearson, 'Optimal mode decomposition for 
% high-dimensional systems. In 51st IEEE Conference on Decision 
% and Control. Maui, Hawaii, Dec. 2012.  Available at 
% http:\\control.ee.ethz.ch\~goularpa\.
%
% Wynn, Pearson, Ganapathisubramani & Goulart, 'Optimal mode 
% decomposition for unsteady and turbulent flows', June 2012.
% Submitted to Journal of Fluid Mechanics. Available at 
% http:\\control.ee.ethz.ch\~goularpa\  
%
% Edelman et. al. 'The geometry of algorithms with orthogonality 
% constraints' SIAM J. Matrix Anal. Appl. 20(2) 303-353).


%set the default solver options, and merge
%in the user ones
opts.relTol  = 1e-5;
opts.maxIter = 500;
opts.maxStep = 0.1;
opts.lineSearch = [];  %leaves line search defaults to linesearch function
if(nargin == 6)
    opts = setstructfields(opts,userOpts);
end


%Check that a solver was specified
if(nargin < 5 || isempty(method))
    method = 'alternating';
end

%Find an initial L if needed
if(nargin < 4 || isempty(L0)) 
    
    %Use the initial bases of the POD modes
    %Note that in some problems svd([A B],0)
    %might actually be better, but we use
    %B only for consistency with the DMD method
	[L0,~] = svd(B,0);
    L0 = L0(:,1:k);
else
    L0 = orth(L0); %use what was provided
end



switch method
    
    case 'hybrid'
    %----------------- 
        fprintf('\n\nInitiating pre-solve step\n');
        L = lmlOpt_altProject(A,B,L0,opts);
        fprintf('\n\nInitiating refinement step\n');
        [L,M] = omd(A,B,k,L,'conjgrad',opts);
        return;
    
    
    case 'alternating'
    %-----------------    
                
        %call the solver;
        L = lmlOpt_altProject(A,B,L0,opts);
        
    case 'conjgrad'
    %-----------------    
    
        %this method does not work in the special
        %case of 1-D data and k = 1.  In that case,
        %call the gradient method
        if(k == 1 && size(A,1) == 1)
            fprintf('\nWarning: 1-d problem.  Calling gradient method');
            [L,M] = omd(A,B,1,L0,'gradient');
            return;
        end
                
        %objective function and its gradient
        nA = norm(A,'fro')^2;
        g  = @(L)(sqrt(nA-get_g(L,A,B)^2));
        dg = @(L)(-get_dgdL(L,A,B));
        
        %call the solver
        L = grassOpt_cg(g,dg,L0,opts);
               
        
    case 'gradient'
    %----------------- 
                
        %objective function and its gradient
        nA = norm(A,'fro')^2;
        g  = @(L)(sqrt(nA-get_g(L,A,B)^2));
        dg = @(L)(-get_dgdL(L,A,B));
        
        %call the solver
        L = grassOpt_gradient(g,dg,L0,opts);    
        
        
    case 'dmd'
    %-----------------
        
        %calculates the 'dynamic mode decomposition',
        %i.e. fix L to the first k vectors in L0 and
        %then find the corresponding M
        L = L0(:,1:k);
        
    otherwise
    %-----------------
        error('Unknown solver method')
end


%Compute the optimal M and return;
M = getMfromL(A,B,L);




    

%----------------------------------------------------
%----------------------------------------------------

function n = get_g(L,A,B)

%computes the function g(L), where
%g(L) = \|L'*A*C\|^2, where C = orth(B'*L)

C = orth(B'*L);
n = norm((L'*A)*C,'fro');




%----------------------------------------------------
%----------------------------------------------------

function dgdL = get_dgdL(L,A,B)

%computes the gradient of function g(L), where
%g(L) = \|L'*A*C\|^2, where C = orth(B'*L)

BtL = B'*L;
AtL = A'*L;
U = AtL'*BtL;
V = BtL'*BtL;
iVU = V\(U');

%Find the derivative 
dgdL = + 2*(A*(BtL*iVU) + B*(AtL*iVU'))  ...
       - 2*B*(BtL*iVU*iVU');       

%----------------------------------------------------
%----------------------------------------------------

function M = getMfromL(A,B,L)

% getMfromL : Given input data (A,B,L), computes an 
% optimal M as
%
%   M = (L'AB'L)(L'BB'L)^{-1}
%
% Usage: M = getMfromL(A,B,L)


LtB = L'*B;
LtA = L'*A;
M   = (LtA * LtB') / (LtB * LtB');










