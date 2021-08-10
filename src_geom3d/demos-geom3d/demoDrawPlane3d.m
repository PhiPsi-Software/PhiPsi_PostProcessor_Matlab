function demoDrawPlane3d
%DEMODRAWPLANE3D demo for drawPlane3d
%
%   Example
%     demoDrawPlane3d
%
%   See also
%

% ------
% Author: oqilipo
% Created: 2017-10-13, using R2017b
% Copyright 2017

%% Initialisation
p0 = [1 2 3];
v1 = [1 0 1];
v2 = [0 -1 1];
plane1 = [p0 v1 v2];
plane2 = [p0 v1 planeNormal(plane1)];

%% 1
% Tools_New_Figure
% axis([-10 10 -10 10 -10 10]);
% drawPlane3d(plane1)
% drawPlane3d(gca,plane2)

%% 2
% Tools_New_Figure
% axis([-10 10 -10 10 -10 10]);
% drawPlane3d(plane1,'g')
% drawPlane3d(gca,plane2,[0 1 0])

%% 3
Tools_New_Figure
axis([-10 10 -10 10 -10 10]);
props.Marker='o';
props.MarkerEdgeColor='r';
props.MarkerFaceColor='g';
props.FaceAlpha=0.5;
drawPlane3d(plane1,props,'EdgeColor','m')
drawPlane3d(gca,plane2,props,'EdgeColor','m')

%% 4
% Tools_New_Figure
% axis([-10 10 -10 10 -10 10]);
% drawPlane3d(plane1,'FaceColor','y',props,'EdgeColor','b')
% drawPlane3d(gca, plane2,'FaceColor','y',props,'EdgeColor','b')

%% 5
% Tools_New_Figure
% axis([-10 10 -10 10 -10 10]);
% drawPlane3d(plane1,props)
% drawPlane3d(gca,plane2,props)

% Active Figure control widget (2021-08-01)
% Press q to exit.
% Press r (or double-click) to reset to the initial.
fcw(gca); 