% plot_slices_3D.m
% A MATLAB script used to visualise time slices of the leading spatial
% eigenvector and the maximum of two SEBA vectors produced for the
% switching Double Gyre using the inflated generator. 
% Written by Aleks Badza
% Created: 03/04/2024
% Last Modified: 15/05/2024

set(groot, 'DefaultLineLineWidth', 1, ...
    'DefaultAxesLineWidth', 1, ...
    'DefaultAxesFontSize', 14, ...
    'DefaultTextFontSize', 14, ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultColorbarTickLabelInterpreter', 'latex', ...
    'DefaultAxesTickLabelInterpreter','latex');

filename = "InfGen_Results_SwitchingDoubleGyre.h5";

x_range = h5read(filename,'/x_range');
y_range = h5read(filename,'/y_range');
T_range = h5read(filename,'/T_range');

T_select = T_range(1:2:end);

spacelength = length(x_range)*length(y_range);
N = spacelength*length(T_range);

Eigvecs_all = h5read(filename,'/Eigvecs_Real');
SEBA_all = h5read(filename,'/SEBA');
SEBA_max = max(SEBA_all,[],2);

vecnum = 2;
sebanum = 1;

V_Cutoff = 0.5;
SEBA_Cutoff = 0.35;

figure

ax = nan(1, 2);
tiles = tiledlayout(1, 2);

[T,X,Y] = meshgrid(T_select,x_range,y_range);

ax(1,1) = nexttile;

V = zeros(size(T));

for j = 1:length(T_select)
    
    J = 2*(j-1) + 1;

    ind_first = (J-1)*spacelength+1;
    ind_last = J*spacelength;
    
    v_now = reshape(Eigvecs_all(ind_first:ind_last,vecnum),[length(y_range) length(x_range)])';
    [r,c] = size(v_now);
    
    v_now = v_now(:);
    v_norm = v_now./(max(abs(v_now)));
    inds = find(abs(v_norm)<V_Cutoff);
    v_now(inds) = NaN;
    v_now = reshape(v_now,[r c]);
    
    V(:,j,:) = sqrt(N)*v_now;
    
    v_slices = slice(T,X,Y,V,T_select(j),[],[]);
    v_slices.EdgeColor = 'none';
    hold on
    
end

hold off
view([60 60])
colormap(ax(1,1),bluewhitered)
colorbar
xlabel('$t$')
ylabel('$x$')
zlabel('$y$')

ax(1,2) = nexttile;

S = zeros(size(T));

for j = 1:length(T_select)
    
    J = 2*(j-1) + 1;

    ind_first = (J-1)*spacelength+1;
    ind_last = J*spacelength;
    
    s_now = reshape(SEBA_max(ind_first:ind_last,sebanum),[length(y_range) length(x_range)])';    
    [r,c] = size(s_now);
    
    s_now = s_now(:);
    s_norm = s_now./(max(abs(s_now)));
    inds = find(abs(s_norm)<SEBA_Cutoff);
    s_now(inds) = NaN;
    s_now = reshape(s_now,[r c]);
    
    S(:,j,:) = s_now;
    
    s_slices = slice(T,X,Y,S,T_select(j),[],[]);
    s_slices.EdgeColor = 'none';
    hold on
    
end

hold off
view([60 60])
colormap(ax(1,2),bluewhitered)
colorbar
clim([0 1])
xlabel('$t$')
ylabel('$x$')
zlabel('$y$')

tiles.TileSpacing = 'compact';
 
set(gcf, 'units', 'inches', 'position', [0.6302 2.8385 13.6667 5.8125])
exportgraphics(gcf,'DoubleGyre_3DSlices.png','Resolution',800)