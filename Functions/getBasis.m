% function [vUnitX,vUnitY,vUnitZ] = getBasis(v1,v2, ... <PROPERTIES>)
% CALCULATE EULER ANGLES BETWEEN ANY TWO 3D BASIS
% - define 2 bases consisting of 3D unit vectors 
% - use quaternion matrix to find Euler angles
%
% INPUTS ------------------------------------------------------------------
%  v1 <numeric [1:3, nSample]> temporal evolution of the 1st column vector        
%       1st basis axis will be colinear with this vector 
%  v2 <numeric [1:3, nSample]> temporal evolution of the 2nd column vector 
%       2nd basis axis will lie in the v1 v2 plane
%       3rd basis axis will be perpendicular to the v1 v2 plane 
%
% PROPERTIES --------------------------------------------------------------
%  sAxis <string> define which basis axis to use for v1 and v2 directions
%                (default: v1 determines X, v2 detemines Y)
%  bPlot    <boolean> plot or not            (default: 0)
%  hPlot    <handle>  figure handle
%  iSmpl    <numeric> sample to plot         (default: 1)
%
% OUTPUT ------------------------------------------------------------------
%  vUnitX,vUnitY,vUnitZ <numeric [1:3, nSample]> Basis of unit vectors 
%       X, Y, and Z (rows) per sample (columns)
%
% NOTE ------------------------------------------------------------------
%   Gimbal Lock problem exists, i.e. getEuler gives a +/- pi results 
%   when abs(a(:,2)) is > 90 deg
%
% EXAMPLE
%   v1 = [3;0;0];
%   v2 = [0;0;3];
%   [vUnitX,vUnitY,vUnitZ] = getBasis(v1,v2,'sAxis','xz');
%
% Valeriya Gritsenko ©2011


function [vUnitX,vUnitY,vUnitZ] = getBasis(v1,v2,varargin)

in.bPlot    = false;
in.hPlot    = [];
in.sAxis    = 'xy';
in.iSmpl    = 1;
if ~isempty(varargin)
   for i = 1:numel(varargin)/2
      in.(varargin{(i-1)*2+1}) = varargin{i*2};
   end
end
nAbs1 = sqrt(v1(1,:).^2+v1(2,:).^2+v1(3,:).^2);
nAbs2 = sqrt(v2(1,:).^2+v2(2,:).^2+v2(3,:).^2);
if strcmp(in.sAxis(1),'x') || strcmp(in.sAxis(1),'X')
    vUnitX(1,:) = v1(1,:)./nAbs1; %X coordinate/vector length
    vUnitX(2,:) = v1(2,:)./nAbs1; %Y coordinate/vector length
    vUnitX(3,:) = v1(3,:)./nAbs1; %Z coordinate/vector length
    if strcmp(in.sAxis(2),'y') || strcmp(in.sAxis(2),'Y')
        %create unit vectors coliear with Y
        vUnitY(1,:) = v2(1,:)./nAbs2;
        vUnitY(2,:) = v2(2,:)./nAbs2;
        vUnitY(3,:) = v2(3,:)./nAbs2;       
        %create unit vector coliear with Z
        vUnitZ = cross(vUnitX,vUnitY); %cross product between Z and X vectors, creates perpendicular vector
        vUnitZ(1,:) = vUnitZ(1,:)./sqrt(vUnitZ(1,:).^2+vUnitZ(2,:).^2+vUnitZ(3,:).^2);
        vUnitZ(2,:) = vUnitZ(2,:)./sqrt(vUnitZ(1,:).^2+vUnitZ(2,:).^2+vUnitZ(3,:).^2);
        vUnitZ(3,:) = vUnitZ(3,:)./sqrt(vUnitZ(1,:).^2+vUnitZ(2,:).^2+vUnitZ(3,:).^2);
        %ensure perpendicularity (orthogonality) between all unit vectors in the basis
        vUnitY = cross(vUnitZ,vUnitX); 
    elseif strcmp(in.sAxis(2),'z') || strcmp(in.sAxis(2),'Z')
        %create unit vectors coliear with Z
        vUnitZ(1,:) = v2(1,:)./nAbs2;
        vUnitZ(2,:) = v2(2,:)./nAbs2;
        vUnitZ(3,:) = v2(3,:)./nAbs2;       
        %create unit vector coliear with Z
        vUnitY = cross(vUnitZ,vUnitX); %cross product between Z and X vectors, creates perpendicular vector
        vUnitY(1,:) = vUnitY(1,:)./sqrt(vUnitY(1,:).^2+vUnitY(2,:).^2+vUnitY(3,:).^2);
        vUnitY(2,:) = vUnitY(2,:)./sqrt(vUnitY(1,:).^2+vUnitY(2,:).^2+vUnitY(3,:).^2);
        vUnitY(3,:) = vUnitY(3,:)./sqrt(vUnitY(1,:).^2+vUnitY(2,:).^2+vUnitY(3,:).^2);
        %ensure perpendicularity (orthogonality) between all unit vectors in the basis
        vUnitZ = cross(vUnitX,vUnitY); 
    end
elseif strcmp(in.sAxis(1),'y') || strcmp(in.sAxis(1),'Y')
    vUnitY(1,:) = v1(1,:)./nAbs1; %X coordinate/vector length
    vUnitY(2,:) = v1(2,:)./nAbs1; %Y coordinate/vector length
    vUnitY(3,:) = v1(3,:)./nAbs1; %Z coordinate/vector length
    if strcmp(in.sAxis(2),'x') || strcmp(in.sAxis(2),'X')
        %create unit vectors coliear with Y
        vUnitX(1,:) = v2(1,:)./nAbs2;
        vUnitX(2,:) = v2(2,:)./nAbs2;
        vUnitX(3,:) = v2(3,:)./nAbs2;       
        %create unit vector coliear with Z
        vUnitZ = cross(vUnitX,vUnitY); %cross product between X and Y vectors, creates perpendicular vector
        vUnitZ(1,:) = vUnitZ(1,:)./sqrt(vUnitZ(1,:).^2+vUnitZ(2,:).^2+vUnitZ(3,:).^2);
        vUnitZ(2,:) = vUnitZ(2,:)./sqrt(vUnitZ(1,:).^2+vUnitZ(2,:).^2+vUnitZ(3,:).^2);
        vUnitZ(3,:) = vUnitZ(3,:)./sqrt(vUnitZ(1,:).^2+vUnitZ(2,:).^2+vUnitZ(3,:).^2);
        %ensure perpendicularity (orthogonality) between all unit vectors in the basis
        vUnitX = cross(vUnitY,vUnitZ); 
    elseif strcmp(in.sAxis(2),'z') || strcmp(in.sAxis(2),'Z')
        %create unit vectors coliear with Z
        vUnitZ(1,:) = v2(1,:)./nAbs2;
        vUnitZ(2,:) = v2(2,:)./nAbs2;
        vUnitZ(3,:) = v2(3,:)./nAbs2;       
        %create unit vector coliear with Z
        vUnitX = cross(vUnitY,vUnitZ); %cross product between Z and X vectors, creates perpendicular vector
        vUnitX(1,:) = vUnitX(1,:)./sqrt(vUnitX(1,:).^2+vUnitX(2,:).^2+vUnitX(3,:).^2);
        vUnitX(2,:) = vUnitX(2,:)./sqrt(vUnitX(1,:).^2+vUnitX(2,:).^2+vUnitX(3,:).^2);
        vUnitX(3,:) = vUnitX(3,:)./sqrt(vUnitX(1,:).^2+vUnitX(2,:).^2+vUnitX(3,:).^2);
        %ensure perpendicularity (orthogonality) between all unit vectors in the basis
        vUnitZ = cross(vUnitX,vUnitY); 
    end
elseif strcmp(in.sAxis(1),'z') || strcmp(in.sAxis(1),'Z')
    vUnitZ(1,:) = v1(1,:)./nAbs1; %X coordinate/vector length
    vUnitZ(2,:) = v1(2,:)./nAbs1; %Y coordinate/vector length
    vUnitZ(3,:) = v1(3,:)./nAbs1; %Z coordinate/vector length
    if strcmp(in.sAxis(2),'x') || strcmp(in.sAxis(2),'X')
        %create unit vector coliear with X
        vUnitX(1,:) = v2(1,:)./nAbs2;
        vUnitX(2,:) = v2(2,:)./nAbs2;
        vUnitX(3,:) = v2(3,:)./nAbs2;       
        %create unit vector coliear with Y
        vUnitY = cross(vUnitZ,vUnitX); %cross product between X and Y vectors, creates perpendicular vector
        vUnitY(1,:) = vUnitY(1,:)./sqrt(vUnitY(1,:).^2+vUnitY(2,:).^2+vUnitY(3,:).^2);
        vUnitY(2,:) = vUnitY(2,:)./sqrt(vUnitY(1,:).^2+vUnitY(2,:).^2+vUnitY(3,:).^2);
        vUnitY(3,:) = vUnitY(3,:)./sqrt(vUnitY(1,:).^2+vUnitY(2,:).^2+vUnitY(3,:).^2);
        %ensure perpendicularity (orthogonality) between all unit vectors in the basis
        vUnitX = cross(vUnitY,vUnitZ); 
    elseif strcmp(in.sAxis(2),'y') || strcmp(in.sAxis(2),'Y')
        %create unit vector coliear with Y
        vUnitY(1,:) = v2(1,:)./nAbs2;
        vUnitY(2,:) = v2(2,:)./nAbs2;
        vUnitY(3,:) = v2(3,:)./nAbs2;       
        %create unit vector coliear with X
        vUnitX = cross(vUnitY,vUnitZ); %cross product between Y and Z vectors, creates perpendicular vector
        vUnitX(1,:) = vUnitX(1,:)./sqrt(vUnitX(1,:).^2+vUnitX(2,:).^2+vUnitX(3,:).^2);
        vUnitX(2,:) = vUnitX(2,:)./sqrt(vUnitX(1,:).^2+vUnitX(2,:).^2+vUnitX(3,:).^2);
        vUnitX(3,:) = vUnitX(3,:)./sqrt(vUnitX(1,:).^2+vUnitX(2,:).^2+vUnitX(3,:).^2);
        %ensure perpendicularity (orthogonality) between all unit vectors in the basis
        vUnitY = cross(vUnitZ,vUnitX);
    end
end
if in.bPlot
    if isempty(in.hPlot)
        figure;
        in.hPlot = axes;
    else
        axes(in.hPlot)
    end
    plot3([0,v1(1,in.iSmpl)],[0,v1(2,in.iSmpl)],[0,v1(3,in.iSmpl)],'k')
    hold on
    plot3([0,v2(1,in.iSmpl)],[0,v2(2,in.iSmpl)],[0,v2(3,in.iSmpl)],'--k')
    plot3([0,vUnitX(1,in.iSmpl)],[0,vUnitX(2,in.iSmpl)],[0,vUnitX(3,in.iSmpl)],'r')
    plot3([0,vUnitY(1,in.iSmpl)],[0,vUnitY(2,in.iSmpl)],[0,vUnitY(3,in.iSmpl)],'g')
    plot3([0,vUnitZ(1,in.iSmpl)],[0,vUnitZ(2,in.iSmpl)],[0,vUnitZ(3,in.iSmpl)],'b')
    grid on
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend('v1','v2',in.sAxis(1),in.sAxis(2))
end
