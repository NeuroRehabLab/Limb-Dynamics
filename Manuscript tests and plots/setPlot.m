%__________________________________________________________________________
% function [hFig, hPlot, hText] = setPlot(<PROPERTIES>)
% SET A STANDARD REPORT FIGURE
%  
% PROPERTIES ---------------------- e.g. function(...,'Property',value,...)
%  nRow        <num> number of axis columns              (default: 1)
%  nCol        <num> number of axis rows                 (default: 1)
%  sAnnotation <str> figure annotation 
%  cAnnotation <num> RGB color of annotation             (default: [1 0 0])
%  sTickDir    <str> tick direction 'in' / 'out'         (default: 'out')
%  sInterpeter <str> 'none' / 'tex'                      (default: 'none')
%  cColor      <str/3 x num> axis bacground color        (default: 'none')
%  sBox        <str> set axes box 'on' / 'off'           (default: 'off')
% 
% OUTPUTS -----------------------------------------------------------------
%  hFig        <handle> figure handle
%  hPlot       <handle> axis handles, [1:nRow, 1:nCol] 
%  hText       <handle> annotation handle
%
% EXAMPLE ----------------------------------------------------------------- 
% [hFig, hPlot]  = setPlot('nRow',2,'nCol',3,'sAnnotation','REPORT');
% hPlot          = reshape(hPlot,1,6); 
% 
% NOTE --------------------------------------------------------------------
% 
% See also: 
%  sci_plot_xy - plot any X,Y values with UI connection to monitorTrial
%  plotSignalXY- plot any two signals with stats
%  plotBar     - bar plots with streamlined interface.
%  plotRaster  - plot a raster report with colour-coded auxiliary events
%  plotTrial   - plot signal with fixed/editable events
%  setLegend   - set legend description anywhere in plot 

% Sergiy Yakovenko © 2012-2013

function [hFig, hPlot, hText] = setPlot(varargin)
% INPUT -------------------------------------------------------------------
% Assign default properties
in.nRow        = 1;
in.nCol        = 1;
in.sAnnotation = '';
in.cAnnotation = [1 0 0];
in.nPosition   = [0 .9 1 .1];
in.cEdgeColor  = 'none';
in.sInterpeter = 'none';   % or 'tex'
in.sTickDir    = 'in';     % 'in' / 'out'
in.cColor      = 'none';   % [1 1 1]
in.sBox        = 'off';    % 'on'/'off'
% Assign inputs 
if ~isempty(varargin)
   for i = 1:numel(varargin)/2
      in.(varargin{(i-1)*2+1}) = varargin{i*2};
   end
end

% JAVA MEMORY HANDLING -- PREVENT OUT OF MEMORY JAVA ERROR ----------------
java_memory_check(0,0);

% PLOT --------------------------------------------------------------------
hFig     = figure;
hText    = annotation('textbox',in.nPosition,'EdgeColor',in.cEdgeColor);
set(hText,'String',in.sAnnotation,'Color',in.cAnnotation,'Tag','hGUI','Interpreter',in.sInterpeter)

% SET FIG GUIs
hMenuFig = uicontextmenu;
set(hMenuFig, 'tag','uimenu_reportfig')
uimenu(hMenuFig, 'Label','PRINT PDF landscape', 'Callback', @PrintPDFlandscape);
uimenu(hMenuFig, 'Label','PRINT PDF portrait',  'Callback', @PrintPDFportrait);
uimenu(hMenuFig, 'Label','PRINT EPS landscape', 'Callback', @PrintEPSlandscape);
uimenu(hMenuFig, 'Label','PRINT EPS portrait',  'Callback', @PrintEPSportrait);
uimenu(hMenuFig, 'Label','EXPORT to WORKSPACE', 'Callback', @SaveDataWork);
set(hText,'uicontextmenu',hMenuFig)

hPlot = zeros(in.nRow, in.nCol);
nPlot = 0;
for iCol = 1:in.nCol
   for iRow = 1:in.nRow
      
      nPlot    = nPlot+1;
      nCol     = nPlot - in.nCol*floor((nPlot-1)/in.nCol);
      nRow     = floor((nPlot-1)/in.nCol)+1;
      hPlot(nRow,nCol) = subplot(in.nRow, in.nCol, nPlot);
      set(hPlot(nRow,nCol),...
         'TickDir',in.sTickDir,...
         'Color'  ,in.cColor,...
         'Box'    ,in.sBox)
      
      % HARD HATS
%       nPos     = get(hPlot(nRow,nCol),'Position');
%       nPos     = [nPos(1:2)+nPos(3:4) .02 .01];
%       hButton  = uicontrol(...
%          'Style'     , 'radiobutton',...
%          'Units'     , 'normalized',...
%          'Position'  , nPos,...
%          'String'    , '',...
%          'UserData'  , hPlot(nRow,nCol), ...
%          'Callback'  , @shrink_expand, ...
%          'Value'     , true);
      
%       title([num2str(nRow),':',num2str(nCol)])

   end
end

% ScalePlot(hPlot)
% shrink_expand(hButton) % run to initialize plot

% SUBFUNCTIONS ------------------------------------------------------------
function PrintPDFlandscape(hObject,~,~, varargin)
hFig        = get(get(hObject,'parent'),'parent');
handles     = guihandles(hFig);
sAnnotation = get(handles.hGUI,'string');
if iscell(sAnnotation), sAnnotation = sAnnotation{1}; end
printpdf(hFig,['FIG ',sAnnotation,'.pdf'],'landscape')


function PrintPDFportrait(hObject,~,~, varargin)
hFig        = get(get(hObject,'parent'),'parent');
handles     = guihandles(hFig);
sAnnotation = get(handles.hGUI,'string');
if iscell(sAnnotation), sAnnotation = sAnnotation{1}; end
printpdf(hFig,['FIG ',sAnnotation,'.pdf'],'portrait')


function PrintEPSlandscape(hObject,~,~, varargin)
hFig        = get(get(hObject,'parent'),'parent');
handles     = guihandles(hFig);
sAnnotation = get(handles.hGUI,'string');
if iscell(sAnnotation), sAnnotation = sAnnotation{1}; end
printeps(hFig,['FIG ',sAnnotation,'.eps'],'landscape')

function PrintEPSportrait(hObject,~,~, varargin)
hFig        = get(get(hObject,'parent'),'parent');
handles     = guihandles(hFig);
sAnnotation = get(handles.hGUI,'string');
if iscell(sAnnotation), sAnnotation = sAnnotation{1}; end
printeps(hFig,['FIG ',sAnnotation,'.eps'],'portrait')

% NEW FEATURE. HARD HATS.
function shrink_expand(hObject,varargin)
hFig        = get(hObject,'Parent');
hAxes       = findobj('Parent',hFig,'Type','axes');
nAxes       = numel(hAxes);
hScale      = findobj('Parent',hFig,'Type','uicontrol','Style','radiobutton');

% hAxes & hScale are in the reverse order of creation and require flipping
hAxes       = hAxes(end:-1:1);
hScale      = hScale(end:-1:1);
x0          = .1; % control width for header & footer

% FIND MAXIMIZED & MINIMIZED PANELS
bMax        = zeros(1,nAxes);
nPos        = zeros(4,nAxes);
for i = 1:nAxes
   bMax(i)     = get(hScale(i),'Value');
   nPos(1:4,i) = get(hAxes(i),'Position');
end
nMin        = .05;
nMax        = (1 - x0 - nMin*sum(bMax==0) )/sum(bMax==1);
if nMax==Inf, nMax=0; end

% SET POSITION OF PANELS
for i = 1:nAxes
   
   if ~bMax(i)
      h = nMin - x0/10;
      set(hAxes(i),'visible','off'); % hide plot
   else
      h = nMax - x0/10;
      set(hAxes(i),'visible','on'); % show plot
   end
   x = nPos(1,i);
   y = x0/2 + nMin*sum( bMax((1:nAxes)>i)==0 ) + nMax*sum( bMax((1:nAxes)>i)==1 );
   w = nPos(3,i);
   set(hAxes(i),'Position',[x y w h]); % reposition axes
   set(hScale(i),'Position',[x+w y+h-.01 .02 .01]); % reposition radiobutton

end

function SaveDataWork(hObject, ~, ~, varargin)
% Find all data lines
hLineList   = findobj('tag','data','type','line');

% Find all that come from the current figure
nLine       = numel(hLineList);
bLine       = zeros(nLine,1);
for iLine = 1:nLine
    hAxes = hLineList(iLine).Parent;
    bLine(iLine) = hAxes.Parent==gcf;
end
hLineList   = hLineList(logical(bLine));

if ~isempty(hLineList)
   for iLine = 1:numel(hLineList)
      UserData(iLine) = hLineList(iLine).UserData;
   end
else
   UserData=[];
end
assignin('base','UserData',UserData)
disp('UserData = ')
disp(UserData)

function ScalePlot(hPlot)
npos = get(hPlot,'Position');
for iPlot = 1:numel(hPlot)
    set(hPlot(iPlot),'Position',[npos(iPlot,1:3), mean(npos(:,4))])
end

















