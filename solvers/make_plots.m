function make_plots( E , G , PI , w , s , g )
% MAKE_PLOTS plots average DMD and OMD eigenvalues, growth rate errors and
% percentage differences between the two methods
%
%   Inputs:
%
%   E = DMD and OMD eigenvalues
%   G = DMD and OMD growth rate errors
%   PI = percentage difference between DMD and OMD
%
%   Plots:
%
%   i) If PI is smaller than 2x2 then eigenvalues corresponding to
%   frequency w(1) and noise levels s(i) are plotted
%
%   ii) If PI is at least 2x2 contour plots of G.DMD, G.OMD and PI are
%   given
%
%   Author A. Wynn - 20 June 2012 

close all

% PI smaller than 2x2
if (length(s)<2) || (length(w)<2)

figure(1); clf;

hold on

% plot eigenvalues 
for i=1:length(s)
    plot(real(E.DMD{1,i}),imag(E.DMD{1,i}),...
        'o',...
        'MarkerSize',8,...
        'MarkerFaceColor',[1 1 1],...
        'MarkerEdgeColor',[0.1 0.1 0.1]);
    plot(real(E.OMD{1,i}),imag(E.OMD{1,i}),...
        'o',...
        'MarkerSize',8,...
        'MarkerFaceColor',[0.5 0.5 0.5],...
        'MarkerEdgeColor',[0.1 0.1 0.1]);
    
%plot true eigenvalues (on first time to get legend correct)                                
if ( i == 1 )
    plot([g g],[w -w],...
        'o',...
        'MarkerSize',8,...
        'MarkerFaceColor','k',...
        'MarkerEdgeColor','k');
end


end

                                    
% format axes
xlabel( '$\mathrm{Re}(\lambda)$',...
        'Interpreter','latex',...
        'FontSize',16);
    
ylabel('$\mathrm{Im}(\lambda)$',...
        'Interpreter','latex',...
        'FontSize',16);

h = gca;
set(h,...
        'FontSize',16,...
        'Ylim',[-w-0.5,w+0.75],...              
        'Box','Off');

    
% format legend
legend( '$\lambda_i^\mathrm{DMD}$',...
        '$\lambda_i^\mathrm{OMD}$',...
        '$\lambda_\mathrm{true}$');

h = legend;
set(h,...
        'Interpreter','latex',...
        'FontSize',16,...
        'Orientation','Horizontal');

    
    
    

% if PI larger than 2x2 plot contours
else

    
    
%---------------------------------------    
%---------------------------------------
% DMD growth rate errors

% contour levels (log plot)
c = [0.0002 0.0005 0.001 0.002 0.005 0.01];

figure(2); clf;

% plot countor of DMD errors
contourf(s,w,log10(G.DMD))

% colour settings
shading flat
c1=colorbar;
colormap(flipud(gray))
caxis([log10(c(1)-1e-8) log10(c(end)+1e-5)])
set(c1,...
        'YTick',log10(c(1:end)),...
        'YTickLabel',c(1:end))
    
set(gca,...
        'FontSize',16)

% labels
xlabel( 'noise covariance $\sigma$',...
        'Interpreter','Latex',...
        'FontSize',16)
    
ylabel( 'temporal frequency $\omega$',...
        'Interpreter','Latex',...
        'FontSize',16)
    
ylabel(c1,...
        '$\epsilon_{\rm{DMD}}$',...
        'Interpreter','Latex',...
        'FontSize',16);
    
h = gca;
set(h,...
        'Fontsize',16);


    
    
    
%---------------------------------------
%---------------------------------------
% contour plot of OMD growth rate errors

% contour levels (log plot)
c = [0.0002 0.0005 0.001 0.002 0.005 0.01];

figure(3); clf;

% plot countor of DMD errors
contourf(s,w,log10(G.OMD))

% colour settings
shading flat
c1=colorbar;
colormap(flipud(gray))
caxis([log10(c(1)-1e-8) log10(c(end)+1e-5)])
set(c1,...
        'YTick',log10(c(1:end)),...
        'YTickLabel',c(1:end))
    
set(gca,...
        'FontSize',18)

% labels
xlabel( 'noise covariance $\sigma$',...
        'Interpreter','Latex',...
        'FontSize',16)
    
ylabel( 'temporal frequency $\omega$',...
        'Interpreter','Latex',...
        'FontSize',16)
    
ylabel(c1,...
        '$\epsilon_{\rm{OMD}}$',...
        'Interpreter','Latex',...
        'FontSize',16);
    
h = gca;
set(h,...
        'Fontsize',16);



%---------------------------------------
%---------------------------------------
% percentage improvement of OMD over DMD

figure(4); clf;

% contour levels
c = 0:5:35;

% plot contour of percentage improvment
contourf(s,w,PI,c);

% contour settings
shading flat
c1 = colorbar;
colormap(flipud(gray))

% labels
ylabel(c1,...
        '\% improvement $p_\epsilon$',...
        'Interpreter','Latex',...
        'Fontsize',16);
    
xlabel( 'noise covariance $\sigma$',...
        'Interpreter','Latex',...
        'Fontsize',16);
    
ylabel( 'temporal frequency $\omega$',...
        'Interpreter','Latex',...
        'Fontsize',16);

h = gca;
set(h,...
        'Fontsize',16);

end 
    
end