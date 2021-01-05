% function getEuler_ARM(nData,nRate,nLength,<PROPERTIES>)
% CALCULATE EULER ANGLES OF SHOULDER, ELBOW & WRIST OF ARM
% - uses simMechanics model getEuler_3j.mdl to derive XYZ euler angles
%   from marker positions
%
% INPUTS ------------------------------------------------------------------
%  nData   <numeric> marker data [nSmpl,nMarker*3], where
%           nSmpl are samples, nMarker are X Y Z coordinates per 
%           Shoulder, Elbow, Wrist, and Hand marker
%           NOTE that first 3 marker columns represent offset and are
%           subtracted from the rest of the coordinates
%  nRate   <numeric> signal sampling frequency
%  nDataW  <numeric> marker data [nSmpl,nMarker*3], where
%           nSmpl are samples, nMarker are X Y Z coordinates per 
%           two wrist markers
% OUTPUTS -----------------------------------------------------------------
%
%  nSEul <numeric> shoulder X Y Z angles in rad
%  nEEul <numeric> elbow Y angle in rad
%  nWEul <numeric> wrist Y Z angles in rad
%
% PROPERTIES --------------------------------------------------------------
%  bPlot   <logical> plot stick figure with basis, default is false
%
% EXAMPLE
%   [nSEul,nEEul,nWEul] = getEuler_ARM(nData,480,nDataW,'bPlot',1)
% Valeriya Gritsenko ©2012

function [nEulerS,nEulerE,nEulerW] = getEuler_ARM_140205(nPosRaw,nRate,nDataW,varargin)
in.bPlot    = false;
if ~isempty(varargin)
    for i = 1:numel(varargin)/2
        in.(varargin{(i-1)*2+1}) = varargin{i*2};
    end
end
%% INITIALIZE Arm Model Parameters
[nSmpl,nMarker] = size(nPosRaw);
[nSmplW,nMarkerW] = size(nDataW);
if nSmpl~=nSmplW
    nAdd = nSmpl-nSmplW;
    nDataW(end:end+nAdd,:) = NaN;
end

tTime       = [1:nSmpl]/nRate;
nDOF = nMarker-3;

nPos(1:length(nPosRaw),1:nDOF) = NaN;
nPos(:,[1,4,7])   = nPosRaw(:,[4,7,10])-[nPosRaw(:,1),nPosRaw(:,1),nPosRaw(:,1)]; % Subtract shoulder motion
nPos(:,[2,5,8])   = nPosRaw(:,[5,8,11])-[nPosRaw(:,2),nPosRaw(:,2),nPosRaw(:,2)];
nPos(:,[3,6,9])   = nPosRaw(:,[6,9,12])-[nPosRaw(:,3),nPosRaw(:,3),nPosRaw(:,3)];
% bNaN1    = isnan(nPos(:,1));
% bNaN2    = isnan(nPos(:,4));
% bNaN3    = isnan(nPos(:,7));
% bNaN     = bNaN1 | bNaN2 | bNaN3;
% nKin(bNaN,:,:) = [];
% tTime(bNaN) = [];

if length(tTime)>3
    %%
    if in.bPlot
        figure
        subplot(3,1,1)
        plot(tTime,nPos(:,1),'r',tTime,nPos(:,2),'g',tTime,nPos(:,3),'b')
        title('Pos')
        ylabel('Elbow')
        subplot(3,1,2)
        plot(tTime,nPos(:,4),'r',tTime,nPos(:,5),'g',tTime,nPos(:,6),'b')
        ylabel('Wrist')
        subplot(3,1,3)
        plot(tTime,nPos(:,7),'r',tTime,nPos(:,8),'g',tTime,nPos(:,9),'b')
        ylabel('Hand')
        legend('X','Y','Z')
    end
    %% ORIGIN BASIS
    a3 = 0*pi/180;
    Rz = [cos(a3),  -sin(a3),       0;...
          sin(a3),   cos(a3),       0;...
          0,           0,           1];

    bOrigin1(1:3,1:nSmpl,1:3) = 0;
    bOrigin1(1,:,1) = 1;
    bOrigin1(2,:,2) = 1;
    bOrigin1(3,:,3) = 1;
    bOrigin(:,:,1) = Rz*bOrigin1(:,:,1);
    bOrigin(:,:,2) = Rz*bOrigin1(:,:,2);
    bOrigin(:,:,3) = Rz*bOrigin1(:,:,3);
    %% WRIST & ELBOW joint angles

    vZ = nPos(:,1:3)'-nPos(:,4:6)'; %Radius-Ulna Z
    vTemp1 = nPos(:,1:3)'; %-Z axis vector: from shoulder to elbow
    vTemp2 = nPos(:,4:6)'; %vector from shoulder to wrist joint
    vX = cross(vTemp1,vTemp2); %vector orthogonal to the shoulder-elbow-wrist plane
    [vUnitX,vUnitY,vUnitZ] = getBasis(vZ,vX,'sAxis','zx','bPlot',0);
    vBasis_RU(:,:,1) = vUnitX;                                                                                                 
    vBasis_RU(:,:,2) = vUnitY;
    vBasis_RU(:,:,3) = vUnitZ;

    vZ = -nPos(:,1:3)'; %Humerus Z
    [vUnitX,vUnitY,vUnitZ] = getBasis(vZ,vX,'sAxis','zx','bPlot',0);
    vBasis_H(:,:,1) = vUnitX;                                                                                                 
    vBasis_H(:,:,2) = vUnitY;
    vBasis_H(:,:,3) = vUnitZ;
   
    
    nEulerE = getEulerXYZ(vBasis_RU(:,:,1),vBasis_RU(:,:,2),vBasis_RU(:,:,3),...
        vBasis_H(:,:,1),vBasis_H(:,:,2),vBasis_H(:,:,3),'sUnit','rad');
    
    iGymbal = find(nEulerE(1,:)>pi);
    if ~isempty(iGymbal) %correction for gymbal lock: rotation of basis
        nRotAng = -pi/2;
         nR = [1 0 0 ;...
             0 cos(nRotAng) -sin(nRotAng);...
             0 sin(nRotAng) cos(nRotAng)];
        vBasis_RU2 = vBasis_RU;
        for iSmpl = iGymbal(1):iGymbal(end)
            tmp = nR*reshape(vBasis_RU(:,iSmpl,:),3,3)';
            vBasis_RU2(:,iSmpl,1) = tmp(1,:)';
            vBasis_RU2(:,iSmpl,2) = tmp(2,:)';
            vBasis_RU2(:,iSmpl,3) = tmp(3,:)';
        end
        nEulerE1 = getEulerXYZ(vBasis_RU2(:,:,1),vBasis_RU2(:,:,2),vBasis_RU2(:,:,3),...
            vBasis_H(:,:,1),vBasis_H(:,:,2),vBasis_H(:,:,3),'sUnit','rad');
        nEulerE(1,iGymbal) = nEulerE1(1,iGymbal)-nRotAng;
        nEulerE(2,iGymbal) = nEulerE1(2,iGymbal);
        nEulerE(3,iGymbal) = nEulerE1(3,iGymbal);
    end
    %% SHOULDER ANGLE
    nEulerS = getEulerXYZ(vBasis_H(:,:,1),vBasis_H(:,:,2),vBasis_H(:,:,3),...
        bOrigin(:,:,1),bOrigin(:,:,2),bOrigin(:,:,3),'sUnit','rad');
    iGymbal = find(nEulerS(3,:)>pi);
    if ~isempty(iGymbal) %correction for gymbal lock: rotation of basis
        nRotAng = pi/2;
         nR = [cos(nRotAng) -sin(nRotAng) 0;...
               sin(nRotAng) cos(nRotAng)  0;...
               0 0 1];
        vBasis_H2 = vBasis_H;
        for iSmpl = iGymbal(1):iGymbal(end)
            tmp = nR*reshape(vBasis_H(:,iSmpl,:),3,3)';
            vBasis_H2(:,iSmpl,1) = tmp(1,:)';
            vBasis_H2(:,iSmpl,2) = tmp(2,:)';
            vBasis_H2(:,iSmpl,3) = tmp(3,:)';
        end
        nEulerS1 = getEulerXYZ(vBasis_H2(:,:,1),vBasis_H2(:,:,2),vBasis_H2(:,:,3),...
            bOrigin(:,:,1),bOrigin(:,:,2),bOrigin(:,:,3),'sUnit','rad');
        nEulerS(1,iGymbal) = nEulerS1(1,iGymbal);
        nEulerS(2,iGymbal) = nEulerS1(2,iGymbal);
        nEulerS(3,iGymbal) = nEulerS1(3,iGymbal)-nRotAng;
    end
    %% WRIST joint angles
    
    vZ = nPos(:,4:6)'-nPos(:,7:9)';
    vX = nDataW(:,4:6)'-nDataW(:,1:3)';
    [vUnitX,vUnitY,vUnitZ] = getBasis(vZ,vX,'sAxis','zx','bPlot',0);
    vBasis_W(:,:,1) = vUnitX;                                                                                                 
    vBasis_W(:,:,2) = vUnitY;
    vBasis_W(:,:,3) = vUnitZ;
    
    nEulerW = getEulerXYZ(vBasis_W(:,:,1),vBasis_W(:,:,2),vBasis_W(:,:,3),...
        vBasis_RU(:,:,1),vBasis_RU(:,:,2),vBasis_RU(:,:,3),'sUnit','rad');
    
    %% PLOT kinematics results
    if in.bPlot
        tTime       = [1:nSmpl]/nRate;
        [hFig, hPlot] = setPlot('nRow',3,'nCol',1,'sAnnotation',...
            'Euler angles');
        plot(hPlot(1),tTime,nEulerS(1,:)*180/pi,'r',...
            tTime,nEulerS(2,:)*180/pi,'g',...
            tTime,nEulerS(3,:)*180/pi,'b')
        legend(hPlot(1),'Eul_X','Eul_Y','Eul_Z')
        ylabel(hPlot(1),'Shoulder, deg')
        plot(hPlot(2),tTime,nEulerE(1,:)*180/pi,'r',...
            tTime,nEulerE(2,:)*180/pi,'g',...
            tTime,nEulerE(3,:)*180/pi,'b')
        ylabel(hPlot(2),'Elbow, deg')
        plot(hPlot(3),tTime,nEulerW(1,:)*180/pi,'r',...
            tTime,nEulerW(2,:)*180/pi,'g',...
            tTime,nEulerW(3,:)*180/pi,'b')
        xlabel(hPlot(3),'Time, s')
        ylabel(hPlot(3),'Wrist, deg')
    end
else
    nEulerS(1:3,1:length(tTime)) = NaN;
    nEulerE(1:3,1:length(tTime)) = NaN;
    nEulerW(1:3,1:length(tTime)) = NaN;
end
