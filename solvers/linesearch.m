function [tout,fout] = linesearch(f,t1,t2,userOpts,f1)

% LINESEARCH : minimizes a 1-d function over an interval
% using a backstepping procedure
%
% Usage : t = linesearch(f,tmin,tmax)
%
% where f : R -> R is a 1 dimensional function, and 
% t is argmin_t f(t) over the interval.
%

% This function implements an approximate forward/backstep procedure

%set the default search options, and merge
%in the user ones
opts.tol              = 1e-5;
opts.backStepScaling  = 0.5;  %backtracking multiplier
if(nargin == 4)
    opts = setstructfields(opts,userOpts);
end

%small sanity check on options
if(opts.backStepScaling <= 0 || opts.backStepScaling >= 1)
    error('Line search scale factor should be in interval (0,1)');
end

%compute the value at t1 if not provided
if(nargin < 5)
    f1 = f(t1);
end
f2 = f(t2);  %value at nominal step size

%initial search window size
tdiff = t2-t1;

%create a list of points evaluated so far
tList = [t1 t2];
fList = [f1 f2];

%Try a forward scaling search
didForwardSearch = 0;
if(f2 < f1)
    didForwardSearch = true;
    [tF,fF,tFex,fFex] = scalingSearch(t1,t2,f,1./opts.backStepScaling,f2);
    tList = [tList tF tFex];
    fList = [fList fF fFex];
end

%if forward search wasn't used,
%try a backwards scaling search
if(~didForwardSearch)
    while((t2-t1)/tdiff > opts.tol)
        if(f2 < f1)
            %it worked!  Continue for best improvement
            [tB,fB,tBex,fBex] = scalingSearch(t1,t2,f,opts.backStepScaling,f2);
            tList = [tList tB tBex];
            fList = [fList fB fBex];
            break;
        else
            %reduce by factor c and try again
            t2 = t1 + (t2-t1)*opts.backStepScaling;
            f2 = f(t2);
        end
    end
end

%One last check - try to fit a quadratic
%function to the points found so far
tQ = quadraticSearch(tList,fList);
fQ = f(tQ);

%assemble everything and find the approximate minimizer
tList = [tList tQ];
fList = [fList fQ];
    
[fout,idx] = min(fList);
tout = tList(idx);


%----------------------------------------
%----------------------------------------

function [tout,fout,tExtreme,fExtreme] = scalingSearch(t1,t2,f,c,f2Init)

%helper function that keeps scaling until 
%the improvement stops

f2 = f2Init; f2last = realmax;

while(f2 <= f2last)
    f2last = f2; t2last = t2;
    %scale by a factor c and try again
    t2 = t1 + (t2-t1).*c;
    f2 = f(t2);
end

%the approximate optimizer
tout = t2last;
fout = f2last;

%the last point evaluated
tExtreme = t2;
fExtreme = f2;


%----------------------------------------
%----------------------------------------

function tOpt = quadraticSearch(t,f)

%Try to fit a quadratic function and
%find the minimizer
A = [t(:).^2 t(:) ones(size(t(:)))];

if(rank(A)<3)
    tOpt = t(1);
else
    c = A\f(:);
    tOpt = -c(2)/c(1)/2;
end

%further protection against ill conditioning
%and going backwards
if(isinf(tOpt) || isnan(tOpt) || tOpt < min(t))
    tOpt = t(1);
end









