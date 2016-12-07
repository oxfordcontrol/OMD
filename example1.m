function [ E_mean G_mean PI_mean ] = example1
% EXAMPLE calculates DMD and OMD eigenvalues, growth rate errors and
% the percentage difference between DMD and OMD when applied to a
% sinusoidal flow. 
%
%   Snapshots are taken from the sinousoidal flow
%
%                 f(x,t) = sin(k*x-w*t)*exp(g*t) + noise
%   
%   with: sampling times (t_i) = [0:dt:dt*Nt]
%         spatial window (x_i) = [0:dx:dx*Nx] 
%
%   Snapshots are arranges into 'before' and 'after' data ensemble matrices
%
%   B = [f(x,t_1) ... f(x,t_Nt)]
%   A = [f(x,t_2) ... f(x,t_Nt+1)]
%
%   DMD and OMD eigenvalues identification techniques are applied to the
%   data ensembles B,A for ranges of temporal frequency 'w' and noise
%   covariance 's'. At each covariance level, 'N' simulations are performed
%   and the average eigenvalues calculated by each method are calculated. 
%
%   Settings:
%
%   Nt = number of snapshots
%   Nx = number of spatial sample points per snapshot
%   dt = time between snapshots
%   dx = spatial sample interval
%
%   w = temporal frequencies
%   s = noise covariances
%   k = spatial wavenumber
%   g = temporal growth rate
%
%   N = number of samples taken at each covariance level
%
%   Outputs:
%
%   E_mean.DMD = average DMD eigenvalues
%   E_mean.OMD = average OMD eigenvalues
%
%   G_mean.DMD = average DMD growth rate errors
%   G_mean.OMD = average OMD growth rate errors
%
%   PI = average percentage difference between OMD and DMD: positive value
%   of PI indicated OMD performs better than DMD. 


%   Author: A. Wynn - 20 June 2012 
%
%   This is an implementation of the setup described in section 4 of
%   the paper:
%
%   Wynn, Pearson, Ganapathisubramani & Goulart, 'Optimal mode 
%   decomposition for unsteady and turbulent flows', June 2012. 
%   Submitted to Journal of Fluid Mechanics. Available at
%   http:\\control.ee.ethz.ch\~goulartpa\
%
%   OMD eigenvalues are calculated using algorithm 1 from the above paper.
%
%   DMD eigenvalues are calculated using the method described in the paper:
%
%   Duke et al., 'An error analysis of the dynamic mode decomposition',
%   Experiments in Fluids 52, 529 (2012).




% settings


% spatial wavenumber and temporal growth rate
k = 1;
g = 1;

% flow sampling settings
Nt = 50;
Nx = 200;
dt = 2*pi/100;
dx = 2*pi/100;

% temporal frequency
w = 2;

% noise covariances
s = 0.1:0.2:1;

% number of samples at each noise level 
N = 250;



% eigenvalue calculations


% holder for averge DMD and OMD eigenvalues
E_mean = struct('DMD',{{}},'OMD',{{}});


% for each w and each s calculate average DMD and OMD eigenvalues for N
% noise samples
for wi=1:length(w)
for si=1:length(s)

E_DMD_tmp = [0;0];
E_OMD_tmp = [0;0];

    % mean of N noise samples at covariance level s
    for n=1:N        
    % get DMD and OMD eigenvalues for each noise sample 
    [E_DMD,E_OMD]  = get_eigs(w(wi),s(si),k,g,Nt,Nx,dt,dx);

    % add to mean
    E_DMD_tmp = E_DMD_tmp + E_DMD/N;
    E_OMD_tmp = E_OMD_tmp + E_OMD/N;
    end

% average DMD and OMD eigenvalues for a given (w,s)
E_mean.DMD{wi,si} = E_DMD_tmp;
E_mean.OMD{wi,si} = E_OMD_tmp;

end
end





% get growth rate errors
G_mean = growth_rate_error(E_mean.DMD,E_mean.OMD,g);


% get percentage difference between OMD over DMD
PI_mean = percent_improve(G_mean.DMD,G_mean.OMD);

fprintf(1,'\n\n');
fprintf(1,'----------------------------------------------------\n');
fprintf(1,'percentage improvement of OMD over DMD:');
fprintf(1,'\n');
PI_mean


% plot results
make_plots(E_mean,G_mean,PI_mean,w,s,g);






% --------------------------------------------
% --------------------------------------------
function [E_DMD,E_OMD] = get_eigs(w,s,k,g,Nt,Nx,dt,dx)
% GET_EIGS Calculates DMD and OMD eigenvalues of snapshot data sampled from
% a siodoidal flow
%   
%   Snapshots are taken from the sinosoidal flow
%
%                 f(x,t) = sin(k*x-w*t)*exp(g*t) + noise
%   
%   with: sampling times (t_i) = [0:dt:dt*Nt]
%         spatial window (x_i) = [0:dx:dx*Nx] 
%
%   Snapshots are arranges into 'before' and 'after' data ensemble matrices
%
%   B = [f(x,t_1) ... f(x,t_Nt)]
%   A = [f(x,t_2) ... f(x,t_Nt+1)] 
%
%   Each snapshot is corrupted with Gaussian noise of covariance s. 
%
%   Inputs:
%
%   Nt = number of snapshots
%   Nx = number of spatial sample points per snapshot
%   dt = time between snapshots
%   dx = spatial sample interval
%
%   w = temporal frequencies
%   s = noise covariances
%   k = spatial wavenumber
%   g = temporal growth rate
%
%   Outputs:
%   
%   E_DMD  = DMD eigenvalues
%   E_OMD  = OMD eigenvalues
%

% generate sineosoidal flow data    
[A,B] = getSine(k,w,g,dt,Nt,dx,Nx,s);
        
% calculate the DMD and OMD eigenvalues
E_DMD = DMD_eigs(A,B,dt,2);
E_OMD = OMD_eigs(A,B,dt,2);






% --------------------------------------------
function [A,B] = getSine(k,w,g,dt,Nt,dx,Nx,s)
% GET_SINE creates before and after data matrices B and A from the
% sinosoidal flow f(x,t)
%
%   snapshot data is sampled from
%
%           f(x,t) := sin(k*x-w*t)*exp(g*t) + noise
%
%   with parameters:
%
%   Nt = number of snapshots
%   Nx = number of spatial sample points per snapshot
%   dt = time between snapshots
%   dx = spatial sample interval
%
%   w = temporal frequency
%   s = noise covariance
%   k = spatial wavenumber
%   g = temporal growth rate
%
%   outputs:
%
%   A = 'after' snapshots
%   B = 'before' snapshots (note A and B overlap)
%
%   columns of A and B are shapshots taken over the temporal window
%   [0:dt:dt*Nt], i.e.
%
%               0 = t_1 < t_2 < ... < t_{Nt+1} = dt*Nt   
%
%   and each snapshot contains velocity information at each point of the
%   spatial window [0:dx:dx*Nx], i.e.
%
%               0 = x_1 < x_2 < ... < x_{Nx+1} = dx*Nx
%
%   That is,
%
%               B = [ u_1 ... u_Nt ],  A = [ u_2 ... u_Nt+1 ] 
%
%   where 
%                        [   f(x_1 ,t_i)   ]
%                u_i =   [        :        ]
%                        [ f(x_{Nx+1},t_i) ]
%


% sample arrays
x = 0:dx:Nx*dx;
t = 0:dt:Nt*dt;

% convert to matrices
[T,X] = meshgrid(t,x);
NU = s.*randn(length(x),length(t));

% sawtooth data set
Y = sin(k.*X - w.*T).*exp(g.*T) + NU;

% data out
B = Y(:,1:end-1); % Before matrix
A = Y(:,2:end);   % After matrix
        


end

% --------------------------------------------
function E_DMD = DMD_eigs(A,B,dt,r)
% DMD_EIGS Calculats r DMD eigenvalues data from data ensembles A 
% and B sampled with timestep dt
%
%   Implements the method described in the paper
%
%   Duke et al., 'An error analysis of the dynamic mode decomposition',
%   Experiments in Fluids 52, 529 (2012).Duke et al. 
%
%   In this implementation of DMD, the singular value matrix is truncated 
%   to contain only values within 10% of the leading singular value. In
%   particular, given snapshot data ensembles, B and A the DMD eigenvalues
%   are the eigenvalues of the matrix 
%
%   Stilde = Ut*A*Vt*Sigt^{-1}
%
%   where B = U*Sig*V' is the compact singular value decomposition of B and
%   Ut, Vt, Sigt are truncated matrices corresponding to the most energetic
%   singular values. Note that this produces the same eigenvalues as the
%   matrix 'S' from Duke et al. equation (2).



[U,Sig,V] = svd(B,0);

% trim to r eigenvectors
Ut   = U(:,1:r);

% calculate dmd eigenvectors 
[L,M] = omd(A,B,r,Ut,'dmd');
eDMD_tmp = eig(M);
eDMD_tmp = sort(eDMD_tmp);

% Calculate DMD eigenvalues
E_DMD = log(eDMD_tmp)./dt;

end

% --------------------------------------------
function  E_OMD  = OMD_eigs(A,B,dt,r)
% OMD_EIGS: calculates r eigenvalues by the OMD method from data ensembles
% A and B sampled with time-step dt
%
%   Calculated r OMD eigenvalues from snapshot data ensembles B and A using
%   algorithm 1 from the paper:
%
%   Wynn, Pearson, Ganapathisubramani & Goulart, 'Optimal mode 
%   decomposition for unsteady and turbulent flows', June 2012. 
%   Submitted to Journal of Fluid Mechanics. Available at
%   http:\\control.ee.ethz.ch\~goulartpa\

% initial guess for L_0
[U,Sig,V] = svd(B,'econ');

% get optimial LML' decomposition using the conjugate-gradient algorithm
[L M]= omd(A,B,r,U(:,1:r),'conjgrad');

% calculate OMD eigenvalues
eOMD_tmp = eig(M);
eOMD_tmp = sort(eOMD_tmp);

% translate eigenvalues to continuous time 
E_OMD = log(eOMD_tmp)./dt;


end


end




% --------------------------------------------
% --------------------------------------------
function G = growth_rate_error( E_DMD , E_OMD , g )
% GROWTH_RATE_ERROR: calculates the growth rate error associated with the
% most unstable eigenvalue calculated by the DMD and OMD methods. 
%
%   If true system has growth rate 'g', the growth rate error associated
%   with an eigenvalue E = a+ib is defined to be 
%
%   growth_rate_error := abs(g-a)/g;
%
%   if error is larger than 'upp' or smaller than 'low' growth rate error
%   is set to zero.
%
%   Inputs:
%
%   E_DMD = DMD eigenvalues
%   E_OMD = OMD eigenvalues
%   g = true growth rate
%
%   Outputs:
%
%   G.DMD = growth rate error of DMD eigenvalues
%   G.OMD = growth rate error of OMD eigenvalues


G = struct('DMD',[],'OMD',[]);

upp = 1;
low = 0.0001;

Nw = size(E_DMD,1);
Ns = size(E_DMD,2);

for i=1:Nw
for j=1:Ns
        
g_DMD = [ real(E_DMD{i,j}(1)) , real(E_DMD{i,j}(2)) ];
g_OMD = [ real(E_OMD{i,j}(1)) , real(E_OMD{i,j}(2)) ];
        
min_err_DMD = min( abs(g - g_DMD(1)), abs(g - g_DMD(2)))/g;
min_err_OMD = min( abs(g - g_OMD(1)), abs(g - g_OMD(2)))/g;


    if ( low <= min_err_DMD ) && ( min_err_DMD <= upp )
    
    G.DMD(i,j) = min_err_DMD;
        
    else
        
    G.DMD(i,j) = 0;
        
    end

    
    if ( low <= min_err_OMD ) && ( min_err_OMD <= upp )
        
    G.OMD(i,j) = min_err_OMD;
            
    else
        
    G.OMD(i,j) = 0;
        
    end
       
        
end
end



end




% --------------------------------------------
% --------------------------------------------
function PI = percent_improve( G_DMD , G_OMD )
% PERCENT_IMPROVE: calculates the percentage difference between growth rate
% errors G_DMD and G_OMD
% 
%   output:
%
%                   PI = ( E_DMD - E_OMD ) / E_DMD
%
%   If both E_DMD and E_OMD are below tolerance, PI := 0.
%   If E_DMD is below tolerance and E_OMD is above tolerance, PI := -Inf.
%

% set tolerance
tol = 0.0001;

% data dimensions
Nw = size(G_DMD,1);
Ns = size(G_DMD,2);

for i=1:Nw
for j=1:Ns
    
    if (G_DMD(i,j) <= tol) && (G_OMD(i,j) <= tol)
    
    PI(i,j) = 0;

    elseif (G_DMD(i,j) <= tol) && (G_OMD(i,j) > tol)
    
    PI(i,j) = -inf;
    
    elseif (G_DMD(i,j) > tol) && (G_OMD(i,j) < tol)
        
    PI(i,j) = 0;
    
    else
      
    PI(i,j) = 100*(G_DMD(i,j) - G_OMD(i,j)) / G_DMD(i,j);
    
    end    

end
end


end







end

