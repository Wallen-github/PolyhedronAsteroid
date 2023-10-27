%--------------------------------------------------------------------------
% Name: plot_inc2b
%
% Desc: Given sphericity, present a new shape model
%
% Author: Hai-Shuo Wang
% Affiliation: Univercity of Colorado Boulder, CSML
% Time: 09/17/2023
% Version 1.0:
%--------------------------------------------------------------------------
clc;
clear;
close all;
format LONG;
addpath(genpath("./data/"))

info = load('./data/inc2b_info_v1.txt');
peroid = load('./data/inc2b_period_v1.txt');
y_data = load('./data/inc2b_inc_v1.txt');
info_2B = [info';peroid'];

xvalues = rad2deg(y_data');
yvalues = {'dr'; 'P_l/P_0'};

figure
h = heatmap(xvalues,yvalues,info_2B);
h.Title = 'C/A results';
h.XLabel = 'inclination (deg)';
% h.YLabel = 'v_\infty (km)';
colorbar off
set(gca,Fontsize=20)