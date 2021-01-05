%%_________________________________________________________________________
% function [idList,bList,sFieldOutput] = getMeta(sMeta,Qry,sFieldOutput)
% GET INFORMATION FROM META TABLES
%
% INPUT -------------------------------------------------------------------
%  sMeta    <char> name of meta table e.g. 'metaTrial'
%              or
%           <struc> structure table, e.g. sMeta.sField = 1:3;
%
%  Qry      <struc> contains WHERE query values
%        The fields values can be single or multiple, so that
%        multiple values of a field are combined with OR, 
%        and values across fields are combined with AND.
%     EXAMPLE:
%        getMeta('metaTrial',qry('bTrial',1,'idTrial',1:10)) 
%           is equivalent to a statement:
%           bTrial==1 AND (idTrial==1 OR idTrial==2 OR ... idTrial==10)
%           which returns id of trials from 1 to 10 with bTrial==1
%     LOGICAL OPERATORS ARE SUPPORTED FOR NUMERICAL FIELDS:
%        use 'x' to denote a field variable e.g. 
%        ...'nObstacleX','x~=0 & x<10',...
%        (note: NaN is not supported in the logical expressions)
%     MASK *string* IS SUPPORTED FOR CELL-STRING FIELDS:
%        ...'sName','*Dr.Who*',...
%     NOTE if isempty(Qry) then return all requested values in specified 
%        fields (default: id field), e.g. idTrial = getMeta('metaTrial','')
%
%  sFieldOutput <char> or <cell char> output fields (default: id<meta>)
%     [] --- select all fields
%     When an array of fields is given the output is a structure array with
%     the specified fields.
%
% OUTPUT ------------------------------------------------------------------
%  idList   <numeric> IDs of meta table or <numeric or cell> 
%        if sFieldOutput is specified 
%  bList    <boolean> selected values from meta table
%
% EXAMPLES ----------------------------------------------------------------
%
%  Table.name     = {'Peter','Jane','Spock','Kate'};
%  Table.id       = 1:4;
%  Table.bLogical = [false false true false];
%  getMeta(Table,qry('bLogical',true),'name')
%  >> ans = 'Spock'
%
%  getMeta('metaTrial',qry('idTrial',5)) 
%  >> ans = 5
%
% NOTE --------------------------------------------------------------------
%  Examples of custom anonymous functions using getMeta() & qry()
%     SESSION  = @(sQRY,xField) getMeta('metaSession',qry(sQRY),xField); 
%     TRIAL    = @(sQRY,xField) getMeta('metaTrial',qry(sQRY),xField);
%     idTrial  = TRIAL({'bTrial',1,'idSession',1},'idTrial')
%
% WARNING -----------------------------------------------------------------
%  Order of elements returned by a query follows ID field order. When order
%  matters, specify sFieldOutput as a cell array of required fields
%
% See also: 
% getMeta()         % get information from META-tables
% getEvent()        % get events 
% getSignal()       % get information from SIGNAL-tables
% setMeta()         % set information in META-tables 
% setEvent()        % add events
% setSignal()       % add signal data
% setSession()      % helps to set sessions / subjects / dates

% Sergiy Yakovenko © 2008-2012

% FUTURE FEATURES
% - Multiple queries for sMeta,Qry input arrays

function [idList,bList,sFieldOutput] = getMeta(sMeta,Qry,sFieldOutput)
global dbData 
if nargin==0 
   idList = dbData;
   return; 
end

if ischar(sMeta) % DATABASE DATA TABLE
   
   if nargin==1
      % CHECK THAT TABLE IS IN dbData
      if isfield(dbData,sMeta) % FIELD EXISTS
         idList = dbData.(sMeta);
      else
         idList = [];
      end
      return; 
   end
   if ~isfield(dbData,sMeta); idList =[]; return; end
   if nargin<3
      
      [idList,bList,sFieldOutput] = qryTable(dbData.(sMeta),Qry);
      
   else
      
      if isempty(sFieldOutput) % SELECT ALL AVAILABLE FIELDS
         sFieldOutput = fieldnames(dbData.(sMeta));
      end
      [idList,bList,sFieldOutput] = qryTable(dbData.(sMeta),Qry,sFieldOutput);
   end
   
else % GENERIC DATA TABLE
   
   if nargin<3
      [idList,bList,sFieldOutput] = qryTable(sMeta,Qry);
   else
      
      if isempty(sFieldOutput) % SELECT ALL AVAILABLE FIELDS
         sFieldOutput = fieldnames(sMeta);
      end
      [idList,bList,sFieldOutput] = qryTable(sMeta,Qry,sFieldOutput);
   end   
   
end





















