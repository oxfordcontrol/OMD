function Q = orth(A)

%ORTH   Orthogonalization.
%   Q = ORTH(A) is an orthonormal basis for the range of A.
%   That is, Q'*Q = I, the columns of Q span the same space as 
%   the columns of A, and the number of columns of Q is the 
%   rank of A.

% This code is copied from the standard
% matlab function of the same name, but with
% a different tolerance setting.  The matlab
% code does not seem able to identify which 
% columns of the SVD are a basis for the image,
% because there are a large number of singular
% values just over the canned matlab threshold
%
% This seems to be a particular problem when
% finding the image for a projected matrix, i.e.
% something like (I-UU')A, where U is orthonormal

[U,S] = svd(A,0);
[m,n] = size(A);
if m > 1, s = diag(S);
   elseif m == 1, s = S(1);
   else s = 0;
end

%mathworks value
%tol = max(m,n) * max(s) * eps(class(A));

%new value
tol = sqrt(eps);

r = sum(s > tol);
Q = U(:,1:r);


