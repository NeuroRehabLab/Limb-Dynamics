% function MODEL_getTorques(nEuler,nPosRaw,nRate,nLength,nWeight,<PROPERTIES>)                                                                               
% CALCULATE EULER ANGLES OF SHOULDER, ELBOW & WRIST OF ARM                                                              
% - uses simMechanics model getEuler_3j.mdl to derive XYZ euler angles                                                                                       
%   from marker positions      
%                                                                                                                            
% INPUTS ------------------------------------------------------------------                                                  
%  nEuler   <numeric> Euler angles [nSmpl,nAngle], where nSmpl are samples,
%           nAngle are joint angles: 3 shoulder, 1 elbow, & 3 wrist
% nPosRaw    <numeric> marker data [nSmpl,9], where 
%           nSmpl are samples, nMarker are X Y Z coordinates of 
%           elbow, wrist, & hand with shoulder at [0,0,0]
%           NOTE shoulder needs to be subtracted from hand coordinates
%  nIniPos <numeric> initial arm posture [1,9]: marker coordinates of
%           elbow XYZ, wrist XYZ, and hand XYZ with shoulder subtracted
%  nRate   <numeric> signal sampling frequency
%  nLength <numeric> lengths of humerus, radius/ulna, and hand in meters
%  nWeight <numeric> subject weight in kg
%                                                                                                                            
% OUTPUTS -----------------------------------------------------------------
%  nOut  <cell> consists of the following variables:
%           nSTor <numeric> [nSmpl,3] shoulder X Y Z muscle torques in Nm/rad
%           nETor <numeric> [nSmpl,1] elbow Y muscle torque in Nm/rad
%           nWTor <numeric> [nSmpl,2] wrist Y Z muscle torques in Nm/rad
%           nSTorG <numeric> [nSmpl,3] shoulder X Y Z gravity torques in Nm/rad
%           nETorG <numeric> [nSmpl,1] elbow Y gravity torque in Nm/rad
%           nWTorG <numeric> [nSmpl,2] wrist Y Z gravity torques in Nm/rad
%           nSTorT <numeric> [nSmpl,3] shoulder X Y Z total 1j torques in Nm/rad
%           nETorT <numeric> [nSmpl,1] elbow Y total 1j torque in Nm/rad
%           nWTorT <numeric> [nSmpl,2] wrist Y Z total 1j torques in Nm/rad
%           nSTorIT <numeric> [nSmpl,3] shoulder X Y Z interaction torques in Nm/rad
%           nETorIT <numeric> [nSmpl,1] elbow Y interaction torque in Nm/rad
%           nWTorIT <numeric> [nSmpl,2] wrist Y Z interaction torques in Nm/rad
%           nSTorT0 <numeric> [nSmpl,3] shoulder X Y Z interaction+1j torques in Nm/rad
%           nETorT0 <numeric> [nSmpl,1] elbow Y interaction+1j torque in Nm/rad
%           nWTorT0 <numeric> [nSmpl,2] wrist Y Z interaction+1j torques in Nm/rad
%
% PROPERTIES --------------------------------------------------------------                                                  
%  bPlot   <logical> plot stick figure with basis, default is false 
%                                                                                                                            
% EXAMPLE                                                                                                                    
%   [nOut] = MODEL_getTorques(nEuler,nPosRaw,480,[0.3,0.3,0.2],70)                                                                         
% Valeriya Gritsenko ©2012                                                                                                  
                                                                                                                             
function [nOut] = MODEL_getTorques(nEuler,nPosRaw,nRate,nLength,nSubjWeight,varargin)                                                                                     
in.bPlot    = false;     
if ~isempty(varargin)                                                                                                        
    for i = 1:numel(varargin)/2                                                                                              
        in.(varargin{(i-1)*2+1}) = varargin{i*2};                                                                            
    end                                                                                                                      
end    
%%
nPos(1:length(nPosRaw),1:9) = NaN;
nPos(:,[1,4,7])   = nPosRaw(:,[4,7,10])-[nPosRaw(:,1),nPosRaw(:,1),nPosRaw(:,1)]; % Subtract shoulder motion
nPos(:,[2,5,8])   = nPosRaw(:,[5,8,11])-[nPosRaw(:,2),nPosRaw(:,2),nPosRaw(:,2)];
nPos(:,[3,6,9])   = nPosRaw(:,[6,9,12])-[nPosRaw(:,3),nPosRaw(:,3),nPosRaw(:,3)];
Ekin = nPos(:,1:3);
Wkin = nPos(:,4:6);
Hkin = nPos(:,7:9);
%% INITIALIZE Arm Model Parameters                                                                             
nRateSim    = nRate;                                                                                             
[nSmpl,nDOF] = size(nEuler);
g           = 9.81; %Gravity, m/s^2                                                                            
tTime       = [1:nSmpl]/nRateSim;  
bNaN    = isnan(nEuler(:,1));
nEuler(bNaN,:) = [];
tTime(bNaN) = []; 
nKin(1:length(nEuler),1:3,1:nDOF) = NaN; 
for iSignal = 1:nDOF
    nTmp    = splineKin(butterfilt(nEuler(:,iSignal),nRate,2,'nOrder',2),nRateSim);
%     nTmp    = splineKin(nEuler(:,iSignal),nRateSim);
    bOut    = nTmp(:,2)>2 | nTmp(:,2)<-2;
    nTmp(bOut,2) = NaN;
    nTmp(:,2)    = interpolate(nTmp(:,2));
    if sum(isnan(nTmp(:,2)))~=0
       nTmp(isnan(nTmp(:,2)),2) = nanmean(nTmp(:,2));
    end
    iOut = find(bOut);
    if ~isempty(iOut) && iOut(end) == length(bOut)
        iOut(end) = [];
    end
    nTmp(iOut+1,1) = NaN;
    nTmp(:,1)    = interpolate(nTmp(:,1));
    if sum(isnan(nTmp(:,1)))~=0
       nTmp(isnan(nTmp(:,1)),1) = nanmean(nTmp(:,1));
    end
    bOut    = nTmp(:,3)>10 | nTmp(:,3)<-10;
    nTmp(bOut,3) = NaN;
    nTmp(:,3)    = interpolate(nTmp(:,3));
    if sum(isnan(nTmp(:,3)))~=0
       nTmp(isnan(nTmp(:,3)),3) = nanmean(nTmp(:,3));
    end
    nKin(:,:,iSignal) = nTmp;
%     nKin(:,2,iSignal) = butterfilt(nKin(:,2,iSignal),nRate,6,'nOrder',2);
%     nKin(:,3,iSignal) = butterfilt(nKin(:,3,iSignal),nRate,6,'nOrder',2);
end                                                                                                                               
nKin = real(nKin);
InKinS1.time                 = tTime;                                                                           
InKinS1.signals.values       = nKin(:,1,1);%-nKin(1,1,1);                                                                           
InKinS1.signals.values(:,2)  = nKin(:,2,1);                                                                      
InKinS1.signals.values(:,3)  = nKin(:,3,1);                                                                      
InKinS1.signals.dimensions   = 3;                                                                               
InKinS2.time                 = tTime;                                                                           
InKinS2.signals.values       = nKin(:,1,2);%-nKin(1,1,2);                                                                           
InKinS2.signals.values(:,2)  = nKin(:,2,2);                                                                      
InKinS2.signals.values(:,3)  = nKin(:,3,2);                                                                      
InKinS2.signals.dimensions   = 3; 
InKinS3.time                 = tTime;                                                                           
InKinS3.signals.values       = nKin(:,1,3);%-nKin(1,1,3);                                                                           
InKinS3.signals.values(:,2)  = nKin(:,2,3);                                                                      
InKinS3.signals.values(:,3)  = nKin(:,3,3);                                                                      
InKinS3.signals.dimensions   = 3; 
InKinE2.time                 = tTime;                                                                          
InKinE2.signals.values       = nKin(:,1,4);%-nKin(1,1,4);                                                                         
InKinE2.signals.values(:,2)  = nKin(:,2,4);                                                                    
InKinE2.signals.values(:,3)  = nKin(:,3,4);                                                                                                                                         
InKinE2.signals.dimensions   = 3; 
InKinW1.time                 = tTime;                                                                           
InKinW1.signals.values       = nKin(:,1,5);%-nKin(1,1,5);                                                                           
InKinW1.signals.values(:,2)  = nKin(:,2,5);                                                                      
InKinW1.signals.values(:,3)  = nKin(:,3,5);                                                                      
InKinW1.signals.dimensions   = 3;    
InKinW2.time                 = tTime;                                                                           
InKinW2.signals.values       = nKin(:,1,6);%-nKin(1,1,6);                                                                           
InKinW2.signals.values(:,2)  = nKin(:,2,6);                                                                      
InKinW2.signals.values(:,3)  = nKin(:,3,6);                                                                      
InKinW2.signals.dimensions   = 3;     
InKinW3.time                 = tTime;                                                                           
InKinW3.signals.values       = nKin(:,1,7);%-nKin(1,1,7);                                                                           
InKinW3.signals.values(:,2)  = nKin(:,2,7);                                                                      
InKinW3.signals.values(:,3)  = nKin(:,3,7);                                                                      
InKinW3.signals.dimensions   = 3;     
%%
if in.bPlot
    figure
    subplot(3,3,1)
    plot(tTime,nKin(:,1,1)*180/pi,'r',tTime,nKin(:,1,2)*180/pi,'g',tTime,nKin(:,1,3)*180/pi,'b')
    title('Ang')
    ylabel('Shoulder')
    subplot(3,3,2)
    plot(tTime,nKin(:,2,1)*180/pi,'r',tTime,nKin(:,2,2)*180/pi,'g',tTime,nKin(:,2,3)*180/pi,'b')
    title('Ang Vel')
    subplot(3,3,3)
    plot(tTime,nKin(:,3,1)*180/pi,'r',tTime,nKin(:,3,2)*180/pi,'g',tTime,nKin(:,3,3)*180/pi,'b')
    title('Ang Acc')
    subplot(3,3,4)
    plot(tTime,nKin(:,1,4)*180/pi,'g')
    ylabel('Elbow')
    subplot(3,3,5)
    plot(tTime,nKin(:,2,4)*180/pi,'g')
    subplot(3,3,6)
    plot(tTime,nKin(:,3,4)*180/pi,'g')
    subplot(3,3,7)
    plot(tTime,nKin(:,1,5)*180/pi,'r',tTime,nKin(:,1,6)*180/pi,'g',tTime,nKin(:,1,7)*180/pi,'b')
    ylabel('Wrist')
    subplot(3,3,8)
    plot(tTime,nKin(:,2,5)*180/pi,'r',tTime,nKin(:,2,6)*180/pi,'g',tTime,nKin(:,2,7)*180/pi,'b')
    subplot(3,3,9)
    plot(tTime,nKin(:,3,5)*180/pi,'r',tTime,nKin(:,3,6)*180/pi,'g',tTime,nKin(:,3,7)*180/pi,'b')
    legend('X','Y','Z')
end                                                                                                              
%% MORPHOMETRY (from Winters)                                                                                  
Ls      = nLength(1);                                  
Le      = nLength(2);                                  
Lw      = nLength(3);                                     
viscS         = 0;% 0.35;                                                                                     
viscE         = 0; %0.1;                                                                                      
%elbow: hand                                                                                                   
wrw = 0.006;     %segment weight/total body weight                                                             
cmprw = 0.506; %center fo mass (CM)/segment length @ Proximal                                                  
cmdrw = 0.494; %CM/segment length @ Distal                                                                     
rcmw = 0.297;   %Radius of Gyration/segment length @ Center fo Gravity                                         
rpw = 0.587;      %Radius of Gyration/segment length @ Proximal                                                
rdw = 0.577;      %Radius of Gyration/segment length @ Distal                                                  
rhow = 1.16;      %Density                                                                                     
%elbow: forearm                                                                                                
wre = 0.016;     %segment weight/total body weight                                                             
cmpre = 0.430; %center fo mass (CM)/segment length @ Proximal                                                  
cmdre = 0.570; %CM/segment length @ Distal                                                                     
rcme = 0.303;   %Radius of Gyration/segment length @ Center fo Gravity                                         
rpe = 0.526;      %Radius of Gyration/segment length @ Proximal                                                
rde = 0.647;      %Radius of Gyration/segment length @ Distal                                                  
rhoe = 1.13;      %Density                                                                                     
%shoulder: upper arm                                                                                           
wrs = 0.028;     %segment weight/total body weight                                                             
cmprs = 0.436; %center fo mass (CM)/segment length @ Proximal                                                  
cmdrs = 0.564; %CM/segment length @ Distal                                                                     
rcms = 0.322;   %Radius of Gyration/segment length @ Center fo Gravity                                         
rps = 0.542;      %Radius of Gyration/segment length @ Proximal                                                
rds = 0.645;      %Radius of Gyration/segment length @ Distal                                                  
rhos = 1.07;      %Density                                                                                     
%% SEGMENT INERTIA                                                                                             
ms = nSubjWeight*wrs;   %upper arm mass                                                                
rs = Ls*cmprs;          %distance to CM                                                              
me = nSubjWeight*wre;   %forearm mass                                                                          
re = Le*cmpre;          %distance to CM                                                                         
mw = nSubjWeight*wrw;   %hand mass                                                                          
rw = Lw*cmpre;          %distance to CM                                                                         
                                                                                                               
%inertia approximated by cyllinder with constant density rotated around CM                                     
% xRod = @(mass,height,radius)[...                                                                               
%     mass*height^2/2 ...                                                                       
%     mass*height^2/12+mass*radius^2/4 ...                                                                       
%     mass*radius^2/12++mass*radius^2/4 ...                                                                                        
%     ]'*[1 1 1].*eye(3);                                                                                        
zRod = @(mass,height,radius)[...                                                                             
    mass*height^2/12+mass*radius^2/4 ...                                                                     
    mass*height^2/12+mass*radius^2/4 ...                                                                     
    mass*radius^2/2 ...                                                                                      
    ]'*[1 1 1].*eye(3);                                                                                      
% yRod = @(mass,height,radius)[...                                                                             
%     mass*height^2/12+mass*radius^2/4 ...                                                                     
%     mass*radius^2/2 ...                                                                                      
%     mass*height^2/12+mass*radius^2/4 ...                                                                     
%     ]'*[1 1 1].*eye(3);                                                                                      
% zEllipsoid = @(mass,radiusX,radiusY,radiusZ)[...                                                                             
%     mass*(radiusY^2+radiusZ^2)/5 ...                                                                     
%     mass*(radiusX^2+radiusZ^2)/5 ...                                                                     
%     mass*(radiusX^2+radiusY^2)/5 ...                                                                     
%     ]'*[1 1 1].*eye(3);                                                                                      
%                                                                                                              
                                                                                                               
Is   = zRod(ms,Ls,.03); %(short radius=0.05m)
Ie   = zRod(me,Le,.03); %(short radius=0.03m)
Iw   = zRod(mw,Lw,.03); %(short radius=0.03m)
%% SEGMENT POSITIONS & ICs                                                                                     
Spos    = [0,0,0];                                                                                         
% Epos    = nIniPos(1,1:3);  
% EposCM  = Epos*rs/sqrt(Epos(1)^2+Epos(2)^2+Epos(3)^2);
% Wpos    = nIniPos(1,4:6);     %wrist CC        
% WposL   = Wpos - Epos;
% WposCM  = WposL*re/sqrt(WposL(1)^2+WposL(2)^2+WposL(3)^2);
% Hpos    = nIniPos(1,7:9);     %hand CC        
% HposL   = Hpos - Wpos;
% HposCM  = HposL*rw/sqrt(HposL(1)^2+HposL(2)^2+HposL(3)^2);
Epos    = [0,0,-Ls];
EposCM  = [0,0,-rs];
Wpos    = [0,0,-(Le+Ls)];     %wrist IC
WposL   = Wpos - Epos;
WposCM  = [0,0,-re];
Hpos    = [0,0,-(Lw+Le+Ls)];    %hand IC
HposL   = Hpos - Wpos;
HposCM  = [0,0,-rw];
nIC    = nKin(1,1:2,1);
nIC(2,:)    = nKin(1,1:2,2);
nIC(3,:)    = nKin(1,1:2,3);
nIC(4,:)    = nKin(1,1:2,4);
nIC(5,:)    = nKin(1,1:2,5);
nIC(6,:)    = nKin(1,1:2,6);
nIC(7,:)    = nKin(1,1:2,7);
%% calculate inertia tensor for initial arm position
% bOrigin = [1,0,0;0,1,0;0,0,1];
% vZ = Epos';
% vY = bOrigin(:,2);
% [vUnitX,vUnitY,vUnitZ] = getBasis(vZ,vY,'sAxis','zy','bPlot',0);  
% bH = [vUnitX,vUnitY,vUnitZ];
% nRotMtrx = bH*bOrigin;
% IsR = nRotMtrx*Is*nRotMtrx'; %rotated Is
% vZ = WposL';
% vY = bOrigin(:,2);
% [vUnitX,vUnitY,vUnitZ] = getBasis(vZ,vY,'sAxis','zy','bPlot',0);  
% bRU = [vUnitX,vUnitY,vUnitZ];
% nRotMtrx = bRU*bOrigin;
% IeR = nRotMtrx*Ie*nRotMtrx'; %rotated Ie
% vZ = HposL';
% vY = bOrigin(:,2);
% [vUnitX,vUnitY,vUnitZ] = getBasis(vZ,vY,'sAxis','zy','bPlot',0);  
% bW = [vUnitX,vUnitY,vUnitZ];
% nRotMtrx = bW*bOrigin;
% IwR = nRotMtrx*Iw*nRotMtrx'; %rotated Ie
%% RUN SIMULATION
sModel = 'KINARM_Dyn_Inv';                                                                                                                                                                                             
try
    simMode = get_param(sModel, 'SimulationMode');
catch
    load_system(sModel)
    simMode = get_param(sModel, 'SimulationMode');
end
assignin('base','g',g)
assignin('base','iStep',1/nRateSim)
assignin('base','tSimEnd',tTime(end))
assignin('base','tSimStart',tTime(1))
assignin('base','InKinS1',InKinS1)
assignin('base','InKinS2',InKinS2)
assignin('base','InKinS3',InKinS3)
assignin('base','InKinE2',InKinE2)
assignin('base','InKinW1',InKinW1)
assignin('base','InKinW2',InKinW2)
assignin('base','InKinW3',InKinW3)
assignin('base','nIC',nIC)
assignin('base','Hpos',Hpos)
assignin('base','Epos',Epos)
assignin('base','Spos',Spos)
assignin('base','Wpos',Wpos)
assignin('base','HposCM',HposCM)
assignin('base','EposCM',EposCM)
assignin('base','WposCM',WposCM)
assignin('base','mw',mw)
assignin('base','Iw',Iw)
assignin('base','me',me)
assignin('base','Ie',Ie)
assignin('base','ms',ms)
assignin('base','Is',Is)
% assignin('base','IsR',IsR)
% assignin('base','IeR',IeR)
% assignin('base','IwR',IwR)
assignin('base','IsR',Is)
assignin('base','IeR',Ie)
assignin('base','IwR',Iw)
disp('simulaitons initialized ...')                                                                            

sim(sModel)    
%% Muscle torques in presence of gravity and interaction torques                                                                                                             
% set_param(sModel, 'StopTime', '10')                                                                                                                                                                                         
% simConfig = getActiveConfigSet(sModel);                                                                      
% simConfig.Components(1)                                                                                      
% get_param(sModel, 'FixedStep')                                                                               
% get_param(sModel, 'StopTime')                                                                                
%                                                                                                              
% set_param(sModel, 'FixedStep', 1/nRateSim)                                                                   
% set_param(sModel, 'Solver','ode4');                                                                          
%   
tTimeSim    = STor1Sim.time;
nSTor       = STor1Sim.signals.values;                                                               
nSTor(:,2)  = STor2Sim.signals.values;                                                               
nSTor(:,3)  = STor3Sim.signals.values;                                                               
nETor       = ETor2Sim.signals.values;                                                               
nWTor       = WTor1Sim.signals.values;                                                               
nWTor(:,2)  = WTor2Sim.signals.values;                                                               
nWTor(:,3)  = WTor3Sim.signals.values;                                                               
tTime       = [1:nSmpl]/nRate;  
nSTor = interp1(tTimeSim,nSTor,tTime,'nearest');
nETor = interp1(tTimeSim,nETor,tTime,'nearest');
nWTor = interp1(tTimeSim,nWTor,tTime,'nearest');
nSTor = butterfilt(nSTor,nRate,6,'nOrder',2);
nETor = butterfilt(nETor,nRate,6,'nOrder',2);
nWTor = butterfilt(nWTor,nRate,6,'nOrder',2);
nSTor(bNaN,:) = NaN;
nETor(bNaN) = NaN;
nWTor(bNaN,:) = NaN;
%% simulate without gravity
assignin('base','g',0)
sim(sModel)    
tTimeSim    = STor1Sim.time;
nSTor0       = STor1Sim.signals.values;                                                               
nSTor0(:,2)  = STor2Sim.signals.values;                                                               
nSTor0(:,3)  = STor3Sim.signals.values;                                                               
nETor0       = ETor2Sim.signals.values;                                                               
nWTor0       = WTor1Sim.signals.values;                                                               
nWTor0(:,2)  = WTor2Sim.signals.values;                                                               
nWTor0(:,3)  = WTor3Sim.signals.values;                                                               
nSTor0 = interp1(tTimeSim,nSTor0,tTime,'nearest');
nETor0 = interp1(tTimeSim,nETor0,tTime,'nearest');
nWTor0 = interp1(tTimeSim,nWTor0,tTime,'nearest');
nSTor0 = butterfilt(nSTor0,nRate,6,'nOrder',2);
nETor0 = butterfilt(nETor0,nRate,6,'nOrder',2);
nWTor0 = butterfilt(nWTor0,nRate,6,'nOrder',2);
nSTor0(bNaN,:) = NaN;
nETor0(bNaN) = NaN;
nWTor0(bNaN,:) = NaN;
nSTorG = nSTor-nSTor0; %Passive gravitational torques
nETorG = nETor-nETor0;
nWTorG = nWTor-nWTor0;
HkinSim = interp1(tTimeSim,HkinSim,tTime,'nearest');
HkinSim(bNaN,:) = NaN;
WkinSim = interp1(tTimeSim,WkinSim,tTime,'nearest');
WkinSim(bNaN,:) = NaN;
EkinSim = interp1(tTimeSim,EkinSim,tTime,'nearest');
EkinSim(bNaN,:) = NaN;
%% simulate without ITs & without G
Ls = Ls+Le+Lw;
Le = Le+Lw;
rs = (ms*rs+me*re+mw*rw)/(ms+me+mw);
re = (me*re+mw*rw)/(me+mw);
ms = ms+me+mw;
me = me+mw;
%inertia approximated by cyllinder with constant density rotated around CM                                                                                                                                                                                                                                       
Is   = zRod(ms,Ls,.03); %(short radius=0.03m)                                                                  
Ie   = zRod(me,Le,.03); %(short radius=0.03m)     
% % EposCM  = Epos*rs/sqrt(Epos(1)^2+Epos(2)^2+Epos(3)^2);
% % WposCM  = WposL*re/sqrt(WposL(1)^2+WposL(2)^2+WposL(3)^2);
% % HposCM  = HposL*rw/sqrt(HposL(1)^2+HposL(2)^2+HposL(3)^2);
% EposCM  = [0,0,-rs];
% WposCM  = [0,0,-re];
% HposCM  = [0,0,-rw];
% 
% assignin('base','me',me)
% assignin('base','Ie',Ie)
% assignin('base','ms',ms)
% assignin('base','Is',Is)
% assignin('base','HposCM',HposCM)
% assignin('base','EposCM',EposCM)
% assignin('base','WposCM',WposCM)
% 
% sModel = 'NoITs_Dyn_Inv';                                                                                                                                                                                             
% sim(sModel)    
% tTimeSim    = STor1Sim.time;
% nSTorT       = STor1Sim.signals.values;                                                             
% nSTorT(:,2)  = STor2Sim.signals.values;                                                               
% nSTorT(:,3)  = STor3Sim.signals.values;                                                               
% nETorT       = ETor2Sim.signals.values;                                                               
% nWTorT       = WTor1Sim.signals.values;                                                               
% nWTorT(:,2)  = WTor2Sim.signals.values;                                                               
% nWTorT(:,3)  = WTor3Sim.signals.values;                                                               
% nSTorT = interp1(tTimeSim,nSTorT,tTime,'nearest');
% nETorT = interp1(tTimeSim,nETorT,tTime,'nearest');
% nWTorT = interp1(tTimeSim,nWTorT,tTime,'nearest');
% nSTorT = butterfilt(nSTorT,nRate,6,'nOrder',2);
% nETorT = butterfilt(nETorT,nRate,6,'nOrder',2);
% nWTorT = butterfilt(nWTorT,nRate,6,'nOrder',2);
% nSTorT(bNaN,:) = NaN;
% nETorT(bNaN) = NaN;
% nWTorT(bNaN,:) = NaN;
% nSTorIT = nSTor0-nSTorT;
% nETorIT = nETor0-nETorT;
% nWTorIT = nWTor0-nWTorT;

% EkinCM  = [Ekin(:,1)*rs./sqrt(Ekin(:,1).^2+Ekin(:,2).^2+Ekin(:,3).^2),...
%     Ekin(:,2)*rs./sqrt(Ekin(:,1).^2+Ekin(:,2).^2+Ekin(:,3).^2),...
%     Ekin(:,3)*rs./sqrt(Ekin(:,1).^2+Ekin(:,2).^2+Ekin(:,3).^2)];
% WkinL   = Wkin - Ekin;
% WkinCM  = [WkinL(:,1)*re./sqrt(WkinL(:,1).^2+WkinL(:,2).^2+WkinL(:,3).^2),...
%     WkinL(:,2)*re./sqrt(WkinL(:,1).^2+WkinL(:,2).^2+WkinL(:,3).^2),...
%     WkinL(:,3)*re./sqrt(WkinL(:,1).^2+WkinL(:,2).^2+WkinL(:,3).^2)];
% HkinL   = Hkin - Wkin;
% HkinCM  = [HkinL(:,1)*rw./sqrt(HkinL(:,1).^2+HkinL(:,2).^2+HkinL(:,3).^2),...
%     HkinL(:,2)*rw./sqrt(HkinL(:,1).^2+HkinL(:,2).^2+HkinL(:,3).^2),...
%     HkinL(:,3)*rw./sqrt(HkinL(:,1).^2+HkinL(:,2).^2+HkinL(:,3).^2)];
% SAcc    = reshape(nKin(:,3,1:3),nSmpl,3);
% EAcc    = nKin(:,3,4);
% WAcc    = reshape(nKin(:,3,5:7),nSmpl,3);
% EradCM = sqrt(EkinCM(:,2).^2+EkinCM(:,3).^2);
% EradCM(:,2) = sqrt(EkinCM(:,1).^2+EkinCM(:,3).^2);
% EradCM(:,3) = sqrt(EkinCM(:,1).^2+EkinCM(:,2).^2);
% 
% Ts = Is.*SAcc + cross(EkinCM,ms*(EradCM.*SAcc));
% Te = Ie.*EAcc + cross(WkinCM,me*WkinCM*EAcc);
% Tw = Iw.*WAcc + cross(HkinCM,mw*HkinCM*WAcc);

%% PLOT kinematics results  
%get goodness of reconstruction
tTimePos = [1:length(Hkin)]/nRate;
nDel = (tTime(end)-tTimePos(end))/2;
tTimePos = tTimePos + nDel;
nRMS1 = sqrt(nanmean((HkinSim(:,1)-Hkin(:,1)).^2));
nRMS2 = sqrt(nanmean((HkinSim(:,2)-Hkin(:,2)).^2));
nRMS3 = sqrt(nanmean((HkinSim(:,3)-Hkin(:,3)).^2));
if in.bPlot
    [hFig, hPlot] = setPlot('nRow',3,'nCol',1,'sAnnotation',...
        'joint torques');
    plot(hPlot(1),tTime,nSTor(:,1),'*r',...
        tTime,nSTor(:,2),'*g',...
        tTime,nSTor(:,3),'*b',...
        tTime,nSTorG(:,1),'--r',...
        tTime,nSTorG(:,2),'--g',...
        tTime,nSTorG(:,3),'--b',...
        tTime,nSTor0(:,1),'r',...
        tTime,nSTor0(:,2),'g',...
        tTime,nSTor0(:,3),'b')
    legend(hPlot(1),'Mtor_X','Mtor_Y','Mtor_Z','Gtor_X','Gtor_Y','Gtor_Z','0tor_X','0tor_Y','0tor_Z')
    ylabel(hPlot(1),'Shoulder, Nm/rad')
    plot(hPlot(2),tTime,nETor,'*r',...
        tTime,nETorG,'--r',...
        tTime,nETor0,'r')
    ylabel(hPlot(2),'Elbow, Nm/rad')
    plot(hPlot(3),tTime,nWTor(:,1),'*r',...
        tTime,nWTor(:,2),'*g',...
        tTime,nWTor(:,3),'*b',...
        tTime,nWTorG(:,1),'--r',...
        tTime,nWTorG(:,2),'--g',...
        tTime,nWTorG(:,3),'--b',...
        tTime,nWTor0(:,1),'r',...
        tTime,nWTor0(:,2),'g',...
        tTime,nWTor0(:,3),'b')
    xlabel(hPlot(3),'Time, s')
    ylabel(hPlot(3),'Wrist, Nm/rad')

    [hFig, hPlot] = setPlot('nRow',3,'nCol',1,'sAnnotation',...                                                    
        'Verification of position simulation');
    plot(hPlot(1,1),tTime,HkinSim(:,1),'--r',...
        tTime,HkinSim(:,2),'--b',...
        tTime,HkinSim(:,3),'--g',...
        tTimePos,Hkin(:,1),'r',...
        tTimePos,Hkin(:,2),'b',...
        tTimePos,Hkin(:,3),'g')
    legend('Sim_X','Sim_Y','Sim_Z','Exp_X','Exp_Y','Exp_Z')
    xlabel(hPlot(1,1),'Time, s')
    ylabel(hPlot(1,1),'Hand position, m')
    title(hPlot(1,1),['rms X ',num2str(nRMS1),', ','rms Y ',num2str(nRMS2),', ','rms Z ',num2str(nRMS3),', '])

    plot(hPlot(2,1),tTime,WkinSim(:,1),'--r',...
        tTime,WkinSim(:,2),'--b',...
        tTime,WkinSim(:,3),'--g',...
        tTimePos,Wkin(:,1),'r',...
        tTimePos,Wkin(:,2),'b',...
        tTimePos,Wkin(:,3),'g')
    xlabel(hPlot(2,1),'Time, s')
    ylabel(hPlot(2,1),'Wrist position, m')

    plot(hPlot(3,1),tTime,EkinSim(:,1),'--r',...
        tTime,EkinSim(:,2),'--b',...
        tTime,EkinSim(:,3),'--g',...
        tTimePos,Ekin(:,1),'r',...
        tTimePos,Ekin(:,2),'b',...
        tTimePos,Ekin(:,3),'g')
    xlabel(hPlot(3,1),'Time, s')
    ylabel(hPlot(3,1),'Elbow position, m')
end  
%% outputs
nOut{1} = nSTor(:,1);
nOut{2} = nSTor(:,2);
nOut{3} = nSTor(:,3);
nOut{4} = nETor';
nOut{5} = nWTor(:,1);
nOut{6} = nWTor(:,2);
nOut{7} = nWTor(:,3);
nOut{8} = nSTorG(:,1);
nOut{9} = nSTorG(:,2);
nOut{10} = nSTorG(:,3);
nOut{11} = nETorG';
nOut{12} = nWTorG(:,1);
nOut{13} = nWTorG(:,2);
nOut{14} = nWTorG(:,3);
% nOut{15} = nSTorT(:,1);
% nOut{16} = nSTorT(:,2);
% nOut{17} = nSTorT(:,3);
% nOut{18} = nETorT';
% nOut{19} = nWTorT(:,1);
% nOut{20} = nWTorT(:,2);
% nOut{21} = nWTorT(:,3);
% nOut{22} = nSTorIT(:,1);
% nOut{23} = nSTorIT(:,2);
% nOut{24} = nSTorIT(:,3);
% nOut{25} = nETorIT';
% nOut{26} = nWTorIT(:,1);
% nOut{27} = nWTorIT(:,2);
% nOut{28} = nWTorIT(:,3);
nOut{29} = nSTor0(:,1);
nOut{30} = nSTor0(:,2);
nOut{31} = nSTor0(:,3);
nOut{32} = nETor0';
nOut{33} = nWTor0(:,1);
nOut{34} = nWTor0(:,2);
nOut{35} = nWTor0(:,3);
nOut{36} = (nRMS1+nRMS2+nRMS3)/3;