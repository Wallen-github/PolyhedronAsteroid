%--------------------------------------------------------------------------
% Name: plot_pervel2b
%
% Desc: Given sphericity, present a new shape model
%
% Author: Hai-Shuo Wang
% Affiliation: Univercity of Colorado Boulder, CSML
% Time: 09/09/2023
% Version 1.0:
%--------------------------------------------------------------------------
clc;
clear;
close all;
format LONG;
addpath(genpath("./data/"))

info_2B = flipud(load('./data/pervel2b_info_v3.txt')');
peroid = flipud(load('./data/pervel2b_period_v3.txt')');
row_data = load('./data/pervel2b_velinf_v3.txt');
col_data = load('./data/pervel2b_periapsis_v3.txt');

xvalues = col_data(:,1)';
yvalues = flipud(row_data(1,:)');

figure
h = heatmap(xvalues,yvalues,info_2B,Colormap=pink);
h.Title = 'number of shifting';
h.XLabel = 'q (Earth Radius)';
h.YLabel = 'v_\infty (km)';
set(gca,Fontsize=20)

figure
h = heatmap(xvalues,yvalues,peroid,Colormap=pink);
h.Title = 'Period Ratio (P_l/P_0)';
h.XLabel = 'q (Earth Radius)';
h.YLabel = 'v_\infty (km)';
set(gca,Fontsize=20)

% cdata = [45 60 32; 43 54 76; 32 94 68; 23 95 58];
% xvalues = {'Small','Medium','Large'};
% yvalues = {'Green','Red','Blue','Gray'};
% h = heatmap(xvalues,yvalues,cdata);
% 
% h.Title = 'T-Shirt Orders';
% h.XLabel = 'Sizes';
% h.YLabel = 'Colors';