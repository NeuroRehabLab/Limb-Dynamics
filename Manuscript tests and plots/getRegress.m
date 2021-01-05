%__________________________________________________________________________
% function [B,R,P,hTitle] = getRegress(x,y, ...<PROPERTIES>)
% STANDARD REGRESSION WITH STREAMLINED DB INTERFACE
%
% INPUT -------------------------------------------------------------------
%  [x,y]    <numeric> vectors
%
% PROPERTIES ---------------------- e.g. function(...,'Property',value,...)
%  bPlot       <boolean> plot or not            (default: 0)
%  hPlot       <handle> axes handle             (default; setPlot)
%  idTrialList <numeric> optional matching idTrialList for troubleshooting data
%  cData       <numeric> RGB color of regression lines [nTrial x 3] 
%                                               (default: [0 0 0] = black)  
%  nMarkerSize <numeric> size of plotted markers (default: 6)
%  nMarkerEdgeColor <numeric or string>         (default: 'none')
%  nMarkerFaceColor <numeric or string>         (default: 'none')
%  nLineWidth  <numeric> line width             (default: 1)
%  nAxisLim    <numeric> axis limits            (default: empty, eg [0 1 0 1])
%
% OUTPUT ------------------------------------------------------------------
%  B     <numeric> [nOffset, nSlope]
%  R     <numeric> R-value, Pearson's coef.
%  P     <numeric> P-value of fit (test P<ALPHA)
%  hTitle <handle> handle to stats title
%
% GUI - GRAPHICAL USER INTERFACE ------------------------------------------
%  1) RIGHT-CLICK on AXES to EXPORT DATA TO WORKSPACE
%  2) RIGHT-CLICK on DATA POINTS to JUMP TO monitorTrial 
% (only if idTrialList is defined & monitorTrial is open)
%
% EXAMPLE -----------------------------------------------------------------
%  ROT = @(x) [cos(x) sin(x);-sin(x) cos(x)];
%  XY = [2*randn(1,100);1*randn(1,100)];
%  XY = 2*ROT(-pi/3)*XY+2; 
%  getRegress(XY(1,:),XY(2,:),'bPlot',1,'cData',[1 0 0]) % plot X vs Y in red
%
%   [B,R,P,hTitle] = getRegress(X,Y,...
%       'bPlot' ,1,...
%       'hPlot' ,hPlot(1),...
%       'nMarkerEdgeColor' ,[1 1 1]);
%    axis square
%    axis([0 .1 0 .5])
%    ylabel('Relative Contribution (nu)')
%    xlabel('Target dexterity (au)')
%    title_addstring(hTitle,'Pathway U1') % add additional line to title
% 
% See also: getOutliers, setPlot, setLegend, regress, axisequalsquare, title_addstring

% Sergiy Yakovenko © 2008-2013

% FUTURE FEATURES

function [B,R,P,hTitle] = getRegress(x,y,varargin)
if nargin<2, B=[];R=[];P=[]; return; end
if isempty(x) || isempty(y), B=[];R=[];P=[]; return; end
x = double(x); 
y = double(y);

% Assign defaults and inputs
in.bPlot             = 0;
in.hPlot             = [];
in.cData             = [0 0 0];
in.nMarkerFaceColor   = in.cData;
in.idTrialList       = [];
in.nMarkerSize       = 6;
in.nMarkerEdgeColor  = 'none';
in.nLineWidth        = 1;
in.nAxisLim          = [];
in.bPlotDensity      = 0;
if ~isempty(varargin)
   for i = 1:numel(varargin)/2
      in.(varargin{(i-1)*2+1}) = varargin{i*2};
   end
end
if isempty(in.hPlot) && in.bPlot==1
   [~, in.hPlot] = setPlot; 
end
if numel(x)>size(in.cData,1),
   in.cData = ones(numel(x),1)*in.cData;
end
hTitle = [];

% Discard NaNs & Check numel
bNAN = isnan(x) | isnan(y);
x = x(~bNAN);
y = y(~bNAN);

if ~isempty(in.idTrialList), 
   in.idTrialList = in.idTrialList(~bNAN);
   in.cData = in.cData(~bNAN,:);
end
if numel(x)<3 || numel(y)<3
   disp('WARNING:getRegress: Less than 3 events found')
   B=[];R=[];P=[];
   return
end

if size(y,2)>size(y,1), y=y';x=x';end
[B,bint,r,rint,stats] = regress(y,[ones(length(x),1),x]); % test model y=B1+B2*x+err
% B     - coefs of regression
% BINT  - 95% confidence intervals for B
% R     - R-square statistic
% RINT  - a matrix of intervals that can be used to diagnose outliers.
% If RINT(i,:) does not contain zero,
% then the i-th residual is larger than would be expected, at the 5%
% significance level.  This is evidence that the I-th observation is an
% outlier.
% STATS - a vector containing, in the following order, the R-square statistic,
% the F statistic and p value for the full model, and an estimate of
% the error variance.

R2 = stats(1); % R*R
P = round(10000*stats(3))/10000; % default precision is .01%
R = sign(B(2))*sqrt(R2); % sign of slope * sqrt of Rsq

if in.bPlot
   % SET AXES SETTINGS & GUIs   
   axes(in.hPlot)
   cmenuAxes = uicontextmenu;
   uimenu(cmenuAxes, 'Label','SAVE DATA to WORKSPACE', 'Callback', @SaveDataWork);
   uimenu(cmenuAxes, 'Label','SAVE DATA to WORKSPACE AS ...', 'Callback', @SaveAsDataWork);
   set(in.hPlot,'uicontextmenu',cmenuAxes)
   hold on
   dx          = [min(x),max(x)];
   UserData.B  = B; 
   UserData.R  = R; 
   UserData.R2 = R2; 
   UserData.P  = P;
   UserData.X  = x;
   UserData.Y  = y;
   
   if B(2)>=0
      sTitle = {['Y=',num2str(sprintf('%0.3g',UserData.B(1))),'+',...
         num2str(sprintf('%0.3g',UserData.B(2))),'*X',...
         ' P=',num2str(sprintf('%0.3g',UserData.P))],...
         [' R=',num2str(sprintf('%0.3g',UserData.R)),...
         ' R2=',num2str(sprintf('%0.3g',UserData.R2))]...
         };
   else
      sTitle = {['Y=',num2str(sprintf('%0.3g',UserData.B(1))),...
         num2str(sprintf('%0.3g',UserData.B(2))),'*X',...
         ' P=',num2str(sprintf('%0.3g',UserData.P))],...
         [' R=',num2str(sprintf('%0.3g',UserData.R)),...
         ' R2=',num2str(sprintf('%0.3g',UserData.R2))]...
         };
   end
   plot(dx,B(1)+B(2)*dx,'-',...
      'UserData'     , UserData,...
      'ButtonDown'   , @onLineClick,...
      'LineWidth'    , in.nLineWidth,...
      'Color'        , in.cData(1,:),...
      'Tag'          , 'data')
   if in.bPlotDensity
      if isempty(in.nAxisLim)
         plotDensity(x,y,'bPlot',1,'nBin',100,'hPlot',in.hPlot);
      else
         plotDensity(x,y,'bPlot',1,'nBin',100,'hPlot',in.hPlot,'nAxisLim',in.nAxisLim);
      end
   end
   if ~isempty(in.idTrialList) % plot interactive elements
      nTrial = numel(in.idTrialList);
      cmenu = zeros(nTrial,1);
      for i = 1:nTrial
         sScript  = ['monitorTrial(''set_idTrial'',',num2str(in.idTrialList(i)),');']; % jump to trial
         cmenu(i) = uicontextmenu;
         uimenu(cmenu(i), 'Label',num2str(in.idTrialList(i)), 'Callback', sScript);
         plot(x(i),y(i)       ,'o',...
            'UIContextMenu'   , cmenu(i),...
            'MarkerEdgeColor' ,in.nMarkerEdgeColor,...
            'MarkerFaceColor' ,in.nMarkerFaceColor(i,:),...
            'MarkerSize'      ,in.nMarkerSize)
      end
   else
      plot(x,y,'o',...
         'MarkerEdgeColor' ,in.nMarkerEdgeColor,...
         'MarkerFaceColor' ,in.nMarkerFaceColor(1,:),...
         'MarkerSize'      ,in.nMarkerSize)
   end
   if ~isempty(in.nAxisLim)
      axis(in.nAxisLim);
   end
   set(in.hPlot,'FontSize',9)
   hTitle = title(sTitle,'FontSize',9);
end

% AXIS LEFT-CLICK --------------------------------------------------------
function onLineClick(hObject, ~, ~, varargin)
hLine       = hObject;
sLineStyle  = get(hLine,'LineStyle');
if strcmp(sLineStyle,'--')
   set(hLine,'LineStyle','-');
else
   hLineList   = findobj('parent',gca,'tag','data');
   UserData    = get(hLine,'UserData');
   if UserData.B(2)>0
      sTitle = {['Y=',num2str(UserData.B(1)),'+',num2str(UserData.B(2)),'*X',...
         ' P=',num2str(sprintf('%0.3g',UserData.P))],...
         [' R=',num2str(sprintf('%0.3g',UserData.R)),...
         ' R2=',num2str(sprintf('%0.3g',UserData.R2))]...
         };
   else
      sTitle = {['Y=',num2str(UserData.B(1)),num2str(UserData.B(2)),'*X',...
         ' P=',num2str(sprintf('%0.3g',UserData.P))],...
         [' R=',num2str(sprintf('%0.3g',UserData.R)),...
         ' R2=',num2str(sprintf('%0.3g',UserData.R2))]...
         };
   end
   title(sTitle)
   if ~isempty(hLineList),
      set(hLineList,'LineStyle','-');
      set(hLine,'LineStyle','--');
   end
end


% AXIS RIGHT-CLICK --------------------------------------------------------
function SaveDataWork(~, ~, ~, varargin)
% HARD HATS. WORK IN PROGRESS 
hLineList = findobj('parent',gca,'tag','data');
if ~isempty(hLineList)
   for iLine = 1:numel(hLineList)
      UserData(iLine) = get(hLineList(iLine),'userdata');
   end
else
   UserData = [];
end
assignin('base','UserData',UserData)
disp('UserData = ')
disp(UserData)

function SaveAsDataWork(~, ~, ~, varargin)
% HARD HATS. WORK IN PROGRESS
% ASK FOR VARNAME
sName    = 'SAVE DATA AS ...';
sPrompt  = {'ENTER VARIABLE NAME TO BE ASSIGNED TO DATA IN WORKSPACE:'};
nLine    = 1;
sDefault = {'UserData'};
sVar     = inputdlg(sPrompt,sName,nLine,sDefault);
if isempty(sVar)
   sVar = 'UserData';
else
   sVar = sVar{1};
end
   
hLineList = findobj('parent',gca,'tag','data');
if ~isempty(hLineList)
   for iLine = 1:numel(hLineList)
      UserData(iLine) = get(hLineList(iLine),'userdata');
   end
else
   UserData = [];
end
assignin('base',sVar,UserData)
disp([sVar,' = '])
disp(UserData)


% HISTORY
% 09-Sep-2013 added hTitle output


