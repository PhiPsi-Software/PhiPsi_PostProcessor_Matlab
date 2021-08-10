function demoDrawLine3d
%DEMODRAWLINE3D Demo of drawLine3d
%
%   Example
%     demoDrawLine3d
%
%   See also
%

% ------
% Author: oqilipo
% Created: 2017-10-13, using R2017b
% Copyright 2017
% clear all; close all; clc; format compact;  format long;

% Add path of source files.
addpath('src_fcw')
addpath('src_geom3d')
addpath('src_meshes3d')
addpath('src_phipsi_post_animate')
addpath('src_phipsi_post_cal')
addpath('src_phipsi_post_main')
addpath('src_phipsi_post_plot')
addpath('src_phipsi_post_read')
addpath('src_phipsi_post_tool')		


		
p0 = [1 2 3];
v1 = [1 0 1];
v2 = [0 -1 1];
line1 = [p0 v1];
lines=[v1 p0;v1 v2];

%% 1
Tools_New_Figure
axis([-10 10 -10 10 -10 10]);
drawLine3d(line1)
drawLine3d(gca,lines)
hold on

%% 2
Tools_New_Figure
axis([-10 10 -10 10 -10 10]);
drawLine3d(line1,'r')
drawLine3d(gca,lines,[1 0 0])
hold on

%% 3
Tools_New_Figure
axis([-10 10 -10 10 -10 10]);
props.Marker='o';
props.MarkerFaceColor='g';
props.MarkerEdgeColor='b';
drawLine3d(line1,props,'Color','m')
drawLine3d(gca,lines,props,'Color','m')
%% 4
Tools_New_Figure
axis([-10 10 -10 10 -10 10]);
drawLine3d(line1,'LineWidth',2,props,'LineStyle','-.')
drawLine3d(gca, lines,'LineWidth',2,props,'LineStyle','-.')
%% 5
Tools_New_Figure
axis([-10 10 -10 10 -10 10]);
drawLine3d(gca,line1,props)
drawLine3d(gca,lines,props)

