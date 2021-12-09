clear all;

%% Definition of available stimuli
% labels of stimulation parameters
labelsLight = {'noLight','light_','light__','light___'};
labelsPulses = {'x100ms','x500ms','x1s','x3s','x10s'};

% stimulation protocol structure
p = struct; 
p = setProtocolStructByNameWidthPeriodPulses(p,'x100ms',0.1,3,50); % POZOR JE TO 100ms
p = setProtocolStructByNameWidthPeriodPulses(p,'x500ms',0.5,5,100); %500ms 5s period 100pulses
p = setProtocolStructByNameWidthPeriodPulses(p,'x1s',1,5,100);
p = setProtocolStructByNameWidthPeriodPulses(p,'x3s',1,5,5);
p = setProtocolStructByNameWidthPeriodPulses(p,'x10s',10,17.8,1);



%% get Erikas .xlsx with information about each block of stimulation 
fp_Tstim='Z:\Petranova\SST, CHR2 H134R\table_slicestim_SST_CHR2_H134R.xlsx';

%  ImportOptions_ = detectImportOptions(fp_Tstim, 'NumHeaderLines', 1);
%  TblockXls = readtable(fp_Tstim,ImportOptions_);
%  save('TblockXls','TblockXls'); % save for testing purpose

load('TblockXls','TblockXls');

% %%  get xls table statistics
% Tstats = table('Size',[5 4],'VariableTypes',repmat({'double'},1,4),'VariableNames',labelsLight,'RowNames',labelsPulses);
% 
% for ip=1:numel(labelsPulses)
%     pulseType = labelsPulses{ip};
%     for il=1:numel(labelsLight)
%         lightType = labelsLight{il};
%         contains = strcmpi(TblockXls.(pulseType),'x') & strcmpi(TblockXls.(lightType),'x') & strcmpi(TblockXls.notes,''); % only "x" and only with empty "note" column
%         Tstats{pulseType,lightType}  = numel(find(contains));
%     end
%     
% end
% 



%% get a list with all .smrx files and path to them
p_EEG = 'Z:\Petranova\SST, CHR2 H134R\recordings'; %D:\tempErikaSST
p_EEG = 'D:\tempErikaSST\recordings'; %
flist = dir([p_EEG '\**\*.smrx']);  % file list

p_block_signal = 'D:\tempErikaSST\block_signal__by_IDblock';



Tblock = getTblockTemplate(1000);
Tied = getTiedTemplate(50000);


%% Support functions
function Tblock = getTblockTemplate(Nrows)
% this can be standalone function
% definition of a table of blocks
% this table is the main data type in the analysis
% a row in the table represents a "block" - a part of EEG with parameters:
varTypes = {'double','categorical','categorical','categorical','cell',  'double' ,'logical','categorical', 'categorical',    'double',  'double' ,'logical' };
varNames = {'IDblock', 'FileName',     'date',       'slice'    , 'signal','fs',     'valid',       'light',       'pulse',      'startSec','endSec'  ,'done' };
Tblock = table('Size',[Nrows,numel(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
end

function Tied = getTiedTemplate(Nrows)
% this can be standalone function
% definition of a table of blocks
% this table is the main data type in the analysis
% a row in the table represents a "block" - a part of EEG with parameters:
varTypes = {'double', 'double',          'double',      'double',      'double',           'double',  'double', 'double'  , 'double'};
varNames = {'IDied',        'IDblock',    'startInBlockInd', 'endInBlockInd', 'lastStimBefore',   'HFOwidth_ms','IEDamp',  'HFOfreq','HFOpwr'};
Tied = table('Size',[Nrows,numel(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
end


function p_structure = setProtocolStructByNameWidthPeriodPulses(p_structure, name,widthSec,periodSec,pulses)
        % define stimulation sequence
        p_structure.(name).pulses = pulses;
        p_structure.(name).widthSec = widthSec;
        p_structure.(name).periodSec = periodSec;
        p_structure.(name).totalDurationSec = pulses*periodSec; % total length of stimualtion 
        p_structure.(name).percentHi=(pulses*widthSec)/(pulses*periodSec); % percent of time being hi state
         
end