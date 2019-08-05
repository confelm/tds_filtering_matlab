% setparameters
% script to set default values for all parameters used in analysis

%% define directory parameters
[parameters.dir.code,pathtype] = getmfilepath; % get the live path project is evaluating
options4box.Resize = 'on'; 
options4box.WindowStyle='normal';
dir2use = inputdlg({'Please provide "UserDirectory" [see step 1, page 3 of our Online Appendix] : '},'User Input',[1 100],{''},options4box);
if isempty(dir2use) || ~ischar(dir2use{1,1})
    error('User Input is not Correct')
else
    parameters.dir.project = dir2use{1,1};
end
if ~isequal(parameters.dir.project(end),pathtype)
    parameters.dir.project =  [parameters.dir.project,pathtype]; 
end
clear dir2use options4box
parameters.dir.processed.exchangeRates = [parameters.dir.project,'data',pathtype,'processed',pathtype,'exchange_rates',pathtype]; 
parameters.dir.processed.equities      = [parameters.dir.project,'data',pathtype,'processed',pathtype,'equities',pathtype];
parameters.dir.tmp.equities            = [parameters.dir.project,'data',pathtype,'temp',pathtype,'equities',pathtype];
if ~exist(parameters.dir.processed.equities,'dir')
    mkdir(parameters.dir.processed.equities)
end
if ~exist(parameters.dir.processed.exchangeRates,'dir')
    mkdir(parameters.dir.processed.exchangeRates)
end
if ~exist(parameters.dir.tmp.equities,'dir')
    mkdir(parameters.dir.tmp.equities)
end
% define raw data directories to load data from 
parameters.dir.raw.exchangeRates = [parameters.dir.project,'data',pathtype,'raw',pathtype,'exchange_rates',pathtype]; 
parameters.dir.raw.equities = [parameters.dir.project,'data',pathtype,'raw',pathtype,'equities',pathtype]; 
%% define filtering specific parameters
parameters.filters.adjustmentCutoff.price = 0.05; % Cutoff to price adjusted data differentials (tds vs self adjusted). Default is 0.05
parameters.filters.adjustmentCutoff.volume = 0.1; % Cutoff for Volume adjusted Data Differentials. Default is 0.1
parameters.filters.filteringBase = 'price'; % 'price' or 'ret' for tds.  Default is 'price'. 
parameters.filters.nDaydead = 10; % dead Price Cuttoff (days)- Consider dead if more than x days constant price at the end of sample. Default : 10
parameters.filters.nDayStale = 30; % stale-illiquid price cutoff - Consider missing if price is stale for more than 30 days in a row. Default : 30
parameters.filters.retcutoff = 1; % extreme daily return reversal Cuttoff (leg 1). Default : 1
parameters.filters.retLeadcutoff = -0.5; % extreme daily return reversal Cuttoff (leg 2). Default is -0.5
parameters.filters.lowerVolcutoff = 1*10^-6; % lower volatility cutoff level - This is in sample and tries to catch outliers stocks with problematic returns -. Default is 1*10^-6
parameters.filters.upperVolcutoff = 0.4; % upper volatility Cuttoff level - This is in sample and tries to catch outliers stocks with problematic returns -. Default is 0.4
parameters.filters.priceLevelcutoff = 20; % Percentage of Price Cutoff. Default : 20 / Instead of a pennyThreashold in dollars
parameters.filters.validStockDay = 0.005; % or 0.01. Percentage of stocks valid per day. Below day is classified as non actually trading one. Default = 0.005
parameters.filters.pennyThreashold = 0; % can be set to 1$ or 5$ or 0 to surpass applying any usual Price Level Filter. Default is 0
parameters.filters.longTermIlliquidCutoff = 0.1; % long Term illiquidity cuttoff. Default : 0.1
parameters.filters.shortTermIlliquidCutoff = 0.95; % short Term illiquidity cuttoff. Default : 0.9 
parameters.filters.extremePerformanceCutoff = 0.98; % extreme in Sample Performance. Default : 0.98
parameters.filters.shortPeriodCutoff = 120; % Short Period Validity Threashold. Default : 120
parameters.filters.extraShortCutoff = 20; % Extremely short Period Validity. Default : 20 
parameters.filters.minStocksForDay = 0; % Minimum allowed stocks per day. Default :0
parameters.filters.minDaysForStock = 0; % Minimum allowed daily observations per stock. Default 0 
parameters.filters.commonTracSet = {'ORD';'FULLPAID';'ORDSUBR'}; % Set of TDS common type identifiers (TRAD,TRAC) This should be updated with the results of the filtering03_getAllStockTypes.m. 
parameters.filters.xetraThreashold = 0.2; % Percentage price differentials between Germany and XETRA duplicated (used in the specific handling of Germany data). 
parameters.filters.currencyExceptions = {'zimbabwe';'ecuador';'venezuela';'arab emirates';'bahrain';'palestine';'russian federation'}; % Countries with an exception in their home currencies 
parameters.filters.emuCurrencySymbols = {'AS';'BF';'M';'DR';'FF';'L';'LF';'FL';'PE';'DM';'EP';...
                    [char(163),'E'];'CY';'KK'; ['M',char(163)];'TO';'EK';'LV';'LT'}; % Define Symbols for historical EMU Countries currencies
%% define variable specific parameters 
parameters.variables.timeseries = {'af';'ax';'mv';'nosh';'p';'ri';'up';'WC03501'}; % variables used in analysis
parameters.variables.static = {'DSCD';'BDATE';'ENAME';'EXMNEM';'GEOGN';'ISIN';'ISINID';'LOC';'PCUR';'TRAC';'TRAD';'TYPE'}; 
parameters.variables.currencyNeutral = {'af';'ax';'nosh'}; % variables set in parameters.variables.timeseries, which are not affected by currency 
parameters.variables.pricing.adjusted = {'p'}; % Defines the pricing variables. For these variables adjusted data will be re-created
parameters.variables.pricing.unadjsuted = {'up'}; % Define the relevant unadjusted oricing series datatypes which correspond to the parameters.variables.pricing.adjusted variables
parameters.variables.variables2filter = {'ret';'retx';'p';'po';'ri'}; % this set of variables will be filtered
parameters.variables.stockIdentifier = 'DSCD'; % this will define the stocksUse variable
parameters.variables.marketIdentifier = 'GEOGN'; % this defines the underlying market 
parameters.variables.exchangeIdentifier = 'EXMNEM'; % this identifies the underlying exchange
parameters.variables.string4missing = '#NA'; % define the string used for missing data in raw data of variables
%
%% Embeeded Functions
function [pathmfile,pathtype] = getmfilepath()
% This function gets the path of the executing mfile
filenamePath = mfilename('fullpath'); 
if isunix
    pathtype = '/'; 
else
    pathtype = '\'; 
end
pathmfile = filenamePath(1:find(ismember(filenamePath,pathtype),1,'last'));
end

%% 2019