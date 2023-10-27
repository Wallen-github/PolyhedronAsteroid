%--------------------------------------------------------------------------
% Name: plot_pervelNb
%
% Desc: Given sphericity, present a new shape model
%
% Author: Hai-Shuo Wang
% Affiliation: Univercity of Colorado Boulder, CSML
% Time: 09/18/2023
% Version 1.0:
%--------------------------------------------------------------------------
clc;
clear;
close all;
format LONG;
addpath(genpath("./data/"))

info_NB = load('./data/pervel4b_info.txt');
nbR3 = 4;
row_num = 10;
col_num = 10;

for i=1:col_num
    
    per(i) = info_NB(i + (i-1)*row_num,1);
    vel(i) = info_NB(i,2);

    for j=1:row_num

        period_ratio(j,i) = info_NB(j + (i-1)*col_num,3);
        num_loss(j,i) = info_NB(j + (i-1)*col_num,4);
        num_nomove(j,i) = info_NB(j + (i-1)*col_num,5);
        num_shift(j,i) = info_NB(j + (i-1)*col_num,6);
        num_circle(j,i) = info_NB(j + (i-1)*col_num,7);

        info_matrix(j,i) = num_shift(j,i);
        if num_loss(j,i) == nbR3-1
            info_matrix(j,i) = -1;
        end

    end
end

period_ratio = flipud(period_ratio);
num_loss = flipud(num_loss);
num_nomove = flipud(num_nomove);
num_shift = flipud(num_shift);
num_circle = flipud(num_circle);
info_matrix = flipud(info_matrix);
vel = fliplr(vel);

figure
h = heatmap(per,vel,info_matrix,Colormap=pink);
h.Title = 'number of shifting';
h.XLabel = 'perigee (Earth Radius)';
h.YLabel = 'v_\infty (km)';
set(gca,Fontsize=20)

figure
h = heatmap(per,vel,period_ratio,Colormap=pink);
h.Title = 'Period Ratio (P_l/P_0)';
h.XLabel = 'q (Earth Radius)';
h.YLabel = 'v_\infty (km)';
set(gca,Fontsize=20)

% figure
% contourf(peaks)
% colorbar('Ticks',[-5,-2,1,4,7],...
%          'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'})

% V(:,:,1)=period_ratio;
% V(:,:,6)=num_loss;
% V(:,:,11)=num_nomove;
% V(:,:,16)=num_shift;
% xslice = [];   
% yslice = [];
% zslice = [1,6,11,16];
% figure
% slice(V,xslice,yslice,zslice,'nearest');
% hmo.Annotate = true;


% 
% figure
% h = heatmap(per,vel,num_loss);
% h.Title = 'number of loss';
% h.XLabel = 'q (Earth Radius)';
% h.YLabel = 'v_\infty (km)';
% colorbar off
% set(gca,Fontsize=20)
% 
% figure
% h = heatmap(per,vel,num_nomove);
% h.Title = 'number of no move';
% h.XLabel = 'q (Earth Radius)';
% h.YLabel = 'v_\infty (km)';
% colorbar off
% set(gca,Fontsize=20)
% 
% figure
% h = heatmap(per,vel,num_shift);
% h.Title = 'number of shifting';
% h.XLabel = 'q (Earth Radius)';
% h.YLabel = 'v_\infty (km)';
% colorbar off
% set(gca,Fontsize=20)
% 
% figure
% h = heatmap(per,vel,num_circle);
% h.Title = 'number of circling';
% h.XLabel = 'q (Earth Radius)';
% h.YLabel = 'v_\infty (km)';
% colorbar off
% set(gca,Fontsize=20)

