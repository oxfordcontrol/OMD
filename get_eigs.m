function [ E ] = get_eigs( w , s , k , g , Nt , Nx , dt , dx , r , N )

% Calculate eigenvalues corresponding to each pair (w,s) of temporal
% frequency and noise covariance

% lengths of inputs
Nw = length(w);
Ns = length(s);

% cell arrays to hold eigenvalues
E = struct('DMD',{{}},'OMD',{{}});

% calculate DMD and OMD eigenvalues for each pair (w,s) and average over N
% noise samples of each given covariance
for i=1:Nw
    for j=1:Ns
        for l=1:N
        
        % generate data 
        [A,B] = getSine(k,w(i),g,dt,Nt,dx,Nx,s(j));
        
        % DMD eigenvalues
        eDMD(l,:) = DMD_eigs(A,B,dt,r);

        % OMD eigenvalues
        eOMD(l,:) = OMD_eigs(A,B,dt,r);
            
        end
    
    % take average of eigenvalues over N samples
    E.DMD{i,j} = mean(eDMD,1);
    E.OMD{i,j} = mean(eOMD,1);
    
    clear eDMD;
    clear eOMD;
    
    end
 
    
    
end

% --------------------------------------------
function [A,B] = getSine(k,w,g,dt,Nt,dx,Nx,nu)

%{

creates data from the flow

f(x,t) := sin(k*x-w*t)*exp(g*t) + noise

flow settings: 

k  = spatial wavenumber 
w  = temporal frequency
g  = growth rate

sampling settings:

Nt = number of temporal snapshots
Nx = number of spatial gridpoints
dt = temporal stepsize
dx = spatial stepsize

nu = covariance of added Gaussian measurement noise

outputs:

A = 'after' snapshots
B = 'before' snapshots (note A and B overlap)

%}




% sample arrays
x = 0:dx:Nx*dx;
t = 0:dt:Nt*dt;

% convert to matrices
[T,X] = meshgrid(t,x);
NU = nu.*randn(length(x),length(t));

% sawtooth data set
Y = sin(k.*X - w.*T).*exp(g.*T) + NU;

% data out
B = Y(:,1:end-1); % Before matrix
A = Y(:,2:end);   % After matrix
        


end

% --------------------------------------------
function E_DMD = DMD_eigs(A,B,dt,r)
%
% Calculate r DMD eigenvalues data ensembles A 
% and B sampled with timestep dt
%

[U,Sig,V] = svd(B,0);

% trim to r eigenvectors
Ut   = U(:,1:r);
Vt   = V(:,1:r);
Sigt = Sig(1:r,1:r);

% calculate eigenvectors of Stilde
Stilde = Ut'*A*Vt/ Sigt;
eDMD_tmp = eig(Stilde);
eDMD_tmp = sort(eDMD_tmp);

% Calculate DMD eigenvalues
E_DMD = log(eDMD_tmp)./dt;

end

% --------------------------------------------
function  E_OMD  = OMD_eigs(A,B,dt,r)

% OMD_eigs: calculates r eigenvalues by the OMD method from data ensembles
% A and B sampled with time-step dt

% initial guess for L_0
[U,Sig,V] = svd(B,'econ');

% get optimial LML' decomposition using the conjugate-gradient algorithm
[L M C]= OMD_ConjGrad(A,B,r,U(:,1:r));
[L,M]= omd(A,B,r,U(:,1:r),'gradient');

% calculate OMD eigenvalues
eOMD_tmp = eig(M);
eOMD_tmp = sort(eOMD_tmp);

E_OMD = log(eOMD_tmp)./dt;


end

end
