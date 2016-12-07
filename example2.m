function  example2
% EXAMPLE2 plots DMD and OMD growth rate errors and the percentage 
% difference between DMD and OMD using data from Wynn et al. (2012)
%
%   Plots data corresponding to the example of section 4 in the paper:
%
%   Wynn, Pearson, Ganapathisubramani & Goulart, 'Optimal mode 
%   decomposition for unsteady and turbulent flows', June 2012. 
%   Submitted to Journal of Fluid Mechanics. Available at
%   http:\\control.ee.ethz.ch\~goulartpa\
%
%   An equivalent data set can be created by running example1.m with the
%   following settings:
%
%   spatial wavenumber and temporal growth rate
%   k = 1;
%   g = 1;
%
%   flow sampling settings
%   Nt = 50;
%   Nx = 200;
%   dt = 2*pi/100;
%   dx = 2*pi/100;
%
%   temporal frequency 
%   w = 0.6:0.01:1.6;
%
%   noise covariance
%   s = 0.05:0.025:1;
%
%   number of samples at each noise level 
%   N = 1000;


% load data
load jfm_data.mat

% plot results
make_plots(E_mean,G_mean,PI_mean,w,s,g);




end
