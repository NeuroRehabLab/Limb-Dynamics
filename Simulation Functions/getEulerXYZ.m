% function [a] = getEulerZYX(X1,Y1,Z1,X2,Y2,Z2, ... <PROPERTIES>)
% CALCULATE EULER ANGLES BETWEEN ANY TWO 3D BASES
% - define 2 bases consisting of 3D unit vectors 
% - use rotation matrix to find Euler angles
%
% INPUTS ------------------------------------------------------------------
%  X1,Y1,Z1 <numeric [1:3, nSample]> evolution of the X, Y, & Z unit       
%       vectors that form the first basis; 
%  X2,Y2,Z2 <numeric [1:3, nSample]> evolution of the X, Y, & Z unit       
%       vectors that form the second basis; 
%       Euler angles are defined with respect to this basis 
%
% PROPERTIES --------------------------------------------------------------
%  sUnit <string> define output in 'rad' or 'deg' (default: 'rad')
%  bPlot    <boolean> plot or not                 (default: 0)
%  hPlot    <handle> figure handle
%
% OUTPUT ------------------------------------------------------------------
%  a <numeric [nSample, 1:3]> Rotation angles around X, Y, and Z (columns) 
%       per sample (rows); in radians or degrees
%
% NOTE ------------------------------------------------------------------
%   Gimbal Lock problem exists, i.e. getEuler gives a +/- pi results 
%   when abs(a(:,2)) is > 90 deg
%
% EXAMPLE
%   X1 = [1;0;0];
%   Y1 = [0;1;0];
%   Z1 = [0;0;1];
%   X2 = [0;1;0];
%   Y2 = [0;0;1];
%   Z2 = [1;0;0];
%   a = getEulerZYX(X1,Y1,Z1,X2,Y2,Z2,'sUnit','deg')%
% Valeriya Gritsenko ©2011


function [a] = getEulerXYZ(X1,Y1,Z1,X2,Y2,Z2,varargin)
in.bPlot    = false;
in.hPlot    = [];
in.sUnit    = 'rad';
if ~isempty(varargin)
   for i = 1:numel(varargin)/2
      in.(varargin{(i-1)*2+1}) = varargin{i*2};
   end
end

% Euler rotation matrix:  rotation order ZYX, i.e. rotation around X is first,
%   followed by Y, followed by Z
% http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
%   Rx = [1,             0,                   0;...
%        0,             cos(a1), -sin(a1);...
%        0,             sin(a1),  cos(a1)];
%    
%    Ry = [cos(a2), 0,                   sin(a2);...
%        0,             1,                   0;...
%        -sin(a2), 0,                   cos(a2)];
%    
%    Rz = [cos(a3),-sin(a3),       0;...
%        sin(a3), cos(a3),       0;...
%        0,             0,                   1];
%   n2 = Rz*Ry*Rx*n1;

%% Rotate and plot
[nr,nsmpl] = size(X1);

for iSmpl = 1:nsmpl
    %% Solve angles
    %   Solution is based on rotation matrix:
    %   http://www.kwon3d.com/theory/euler/euler_angles.html
    nA = [X1(:,iSmpl), Y1(:,iSmpl), Z1(:,iSmpl)];
    nB = [X2(:,iSmpl), Y2(:,iSmpl), Z2(:,iSmpl)];
    n1 = nA'*nA;
    n2 = nA'*nB;
    R_Sval = n2/n1;
    a2 = asin(R_Sval(3,1));
    if R_Sval(3,3)*cos(a2)>0
        a1 = atan2(-R_Sval(3,2),R_Sval(3,3));
    else
        a1 = atan2(-R_Sval(3,2),R_Sval(3,3))+pi;
    end
    if R_Sval(1,1)*cos(a2)>0
        a3 = atan2(-R_Sval(2,1),R_Sval(1,1));
    else
        a3 = atan2(-R_Sval(2,1),R_Sval(1,1))+pi;
    end  
    if in.bPlot
        % % Check Euler angles
        R_S = [...
            cos(a2)*cos(a3),     cos(a1)*sin(a3)+sin(a1)*sin(a2)*cos(a3),       sin(a1)*sin(a3)-cos(a1)*sin(a2)*cos(a3);...
            -cos(a2)*sin(a3),    cos(a1)*cos(a3)-sin(a1)*sin(a2)*sin(a3),       sin(a1)*cos(a3)+cos(a1)*sin(a2)*sin(a3);...
            sin(a2),             -sin(a1)*cos(a2),                              cos(a1)*cos(a2)];       
        n2c = R_S*n1;
        if isempty(in.hPlot)
            figure;
            in.hPlot = axes;
        else
            axes(in.hPlot)
        end
        plot3([0,n1(1,1)],[0,n1(2,1)],[0,n1(3,1)],'r')
        hold on
        plot3([0,n2(1,1)],[0,n2(2,1)],[0,n2(3,1)],'--r')
        plot3([0,n2c(1,1)],[0,n2c(2,1)],[0,n2c(3,1)],'or')
        plot3([0,n1(1,2)],[0,n1(2,2)],[0,n1(3,2)],'g')
        plot3([0,n2(1,2)],[0,n2(2,2)],[0,n2(3,2)],'--g')
        plot3([0,n2c(1,2)],[0,n2c(2,2)],[0,n2c(3,2)],'og')
        plot3([0,n1(1,3)],[0,n1(2,3)],[0,n1(3,3)],'b')
        plot3([0,n2(1,3)],[0,n2(2,3)],[0,n2(3,3)],'--b')
        plot3([0,n2c(1,3)],[0,n2c(2,3)],[0,n2c(3,3)],'ob')
        grid on
        axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')
        legend('nA','nB','R*nA')
        pause(0.5)
    end
    if strcmp(in.sUnit,'deg')
        a(1:3,iSmpl) = [a1;a2;a3]*180/pi;
    else
        a(1:3,iSmpl) = [a1;a2;a3];
    end
end
