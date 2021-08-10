function demoDrawPoint3d
%DEMODRAWPOINT3D Demo of drawPoint3d
%
%   Example
%     demoDrawPoint3d
%
%   See also
%

% ------
% Author: oqilipo
% Created: 2018-01-12, using R2017b
% Copyright 2018

%% Initialisation
ptsx = rand(1,500000);
ptsy = rand(1,500000);
ptsz = rand(1,500000);

pts2 = rand(500000,3);

%% 1
% Tools_New_Figure
% drawPoint3d(ptsx, ptsy, ptsz)
% hold on
% drawPoint3d(gca,pts2)

%% 2
% Tools_New_Figure
% drawPoint3d(pts2,'g.')
% hold on
% drawPoint3d(gca,ptsx, ptsy, ptsz,'r')
%% 3
Tools_New_Figure
props.Marker='o';
props.MarkerEdgeColor='r';
props.MarkerFaceColor='g';
props.LineStyle='none';
drawPoint3d(ptsx, ptsy, ptsz,props,'MarkerSize',9)
hold on
drawPoint3d(gca,pts2,props,'MarkerSize',3)

%% 4
% Tools_New_Figure
% drawPoint3d(ptsx, ptsy, ptsz,'MarkerSize',8,props,'Marker','s')
% hold on
% drawPoint3d(gca, pts2,'MarkerSize',12,props,'Marker','h')
%% 5
% Tools_New_Figure
% drawPoint3d(ptsx, ptsy, ptsz,props)
% hold on
% drawPoint3d(gca,pts2,props)


% Active Figure control widget (2021-08-01)
% Press q to exit.
% Press r (or double-click) to reset to the initial.
fcw(gca); 