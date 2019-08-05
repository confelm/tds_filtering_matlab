% process_raw_tds 
% code processes the raw tds data as downloaded based on internet appendix of Landis and Skouras
% (2019). Code is automatically called by the main batfile handling process
clear;close all;clc; 
%% define parameters 
setparameters
%% process tds exchange rates data
% read excel files for the exchange rates data - Sheet 1 sterling rates
% and Sheet2 euro rates. 
% check existency of raw data - user may have not organized data properly - 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Processing raw tds data...')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(datestr(now))
disp(' ')
if ~exist(parameters.dir.raw.exchangeRates,'dir')
    warning('Cannot find data for exchange rates......check readme.m file')
    return
end
[ret,stocksUse,tUse,currencyUse] = readxlsxfiles_tds([parameters.dir.raw.exchangeRates,'exchange_rates.xlsx'],1,'timeseries',parameters.variables.string4missing);  %#ok<ASGLU>
save([parameters.dir.processed.exchangeRates,'sterling_rates.mat'],'ret','stocksUse','tUse','currencyUse','-v7.3')
clear ret stocksUse tUse currencyUse
% read euro related rates from Sheet2
[ret,stocksUse,tUse,currencyUse] =  readxlsxfiles_tds([parameters.dir.raw.exchangeRates,'exchange_rates.xlsx'],2,'timeseries',parameters.variables.string4missing);  %#ok<ASGLU>
% save exchange rates data to directory 
save([parameters.dir.processed.exchangeRates,'euro_rates.mat'],'ret','stocksUse','tUse','currencyUse','-v7.3')
clear ret stocksUse tUse currencyUse
% download the uk2emu rates data. Add into the online appendix if needed
%% process raw tds equities data 
% get all raw files for equities
% All saved in the temp data directory to be removed later. 
if ~exist(parameters.dir.raw.equities,'dir')
    warning('Cannot find raw xlsx data for equities...check readme.m file')
    return
end
% read all files
[filesList] = getfilenames(parameters.dir.raw.equities,'xlsx'); 
if isempty(filesList)
    warning('Cannot find raw xlsx data for equities...check readme.m file')
    return
end
% get the index of all folders for markets 
[tdsListIndex,tdsListName] = cellfun(@(x) getfolderindex(x,pathtype),filesList,'UniformOutput',false); 
tdsListIndex = cell2mat(tdsListIndex);
staticIndex = cell2mat(cellfun(@(x) contains(x,'static.xlsx'),filesList,'UniformOutput',false));
nList = max(tdsListIndex); % get total number of markets to be processed. 
tdsListNameUnique = tdsListName(staticIndex); % get only once for each block
% loop on all lists
tdsMarketList = cell(nList,1); % preallocate variable with market names 
tdsMarketListName = cell(length(tdsListName),1); 
errorIndex = false(length(filesList),1); % initiate an error handler for the process 
for iList =1:nList
  disp('---------------------')
  disp('')
  % define savedirectory for market 
  savedirectory4list = [parameters.dir.tmp.equities,tdsListNameUnique{iList},pathtype]; 
  if ~exist(savedirectory4list,'dir');mkdir(savedirectory4list);end % create dir in the first iteration 
  % loop on all datatypes
  % process static data first 
  files4list = filesList(tdsListIndex==iList); 
  fileName4static  = filesList{tdsListIndex==iList&staticIndex}; 
  files2process4list = setxor(files4list,fileName4static); 
  % process static data at the start of process. 
  if isempty(fileName4static)
      disp(['static data not available for list no : ',num2str(iList),'...'])
      disp('abort processing list....')
      errorIndex(tdsListIndex==iList) = true; % add files to errorneous processed list so user to know at the end
      continue 
  else
      disp(['Processing File : ',fileName4static])
      try
          [static,headers] = readxlsxfiles_tds(fileName4static,1,'static',parameters.variables.string4missing);
      catch errorlog
          disp('Error reading static.xlsx file')
          disp(errorlog.message)
          continue
      end
      ixStock = find(ismember(headers,parameters.variables.stockIdentifier));
      ixMarket = find(ismember(headers,parameters.variables.marketIdentifier)); 
      ixExchange = find(ismember(headers,parameters.variables.exchangeIdentifier)); 
      if isempty(ixStock)|| isempty(ixMarket)
          disp('Invalid parameters for Stock or Market Identifiers')
          disp('abort processing list....')
          errorIndex(tdsListIndex==iList) = true; 
          continue
      else
          stocksUse = static(:,ixStock);
          % Convert to cellstring 
          %nanIndexStock = false(length(stocksUse),1); 
          for istock=1:length(stocksUse)
              % Replace to cellstring variables : 
             if ~ischar(stocksUse{istock})
                 stocksUse{istock} = num2str(stocksUse{istock}); 
             end
          end
          clear istock
          % analysis for exchange
          exchangeUse = static(:,ixExchange); 
          exchangeTabulate = tabulate(exchangeUse); 
          exchangeName = exchangeTabulate{find(cell2mat(exchangeTabulate(:,2))==max(cell2mat(exchangeTabulate(:,2)))),1};  %#ok<FNDSB>
          isxetra = isequal(exchangeName,'XETRA'); 
          % analysis for market 
          marketUse = static(:,ixMarket);  
          marketTabulate = tabulate(marketUse); 
          marketName = marketTabulate{find(cell2mat(marketTabulate(:,2))==max(cell2mat(marketTabulate(:,2)))),1};  %#ok<FNDSB>
          isgermany = isequal(marketUse,'germany'); 
          if isxetra && ~isgermany;warning('Error with XETRA exchange');continue;end
          % apply different names for germany 
          if isxetra && isgermany
              marketName = 'xetra'; 
          elseif ~isxetra && isgermany
              marketName = 'germanyexXetra'; 
          end
          tdsMarketList{iList} = marketName; % This will later give all sublists for a specific market
          tdsMarketListName(tdsListIndex==iList,1) = {marketName};
          disp('')
          disp(['Processig Market/Exchange : ',marketName,'/',exchangeName])
          disp('')
          save([savedirectory4list,'static.mat'],'static','headers','stocksUse','marketName','isxetra','-v7.3')
          clear *Tabulate ixStock ixMarket ixExchange static header stocksUse exchangeUse exchangeName
      end
  end
  % process all other datatypes for the market
  h = waitbar(0,['processing raw data for ',marketName,'... list # : ',tdsListNameUnique{iList}]);
  variablesProcessed4list = {'static'}; 
  for idatatype = 1:length(files2process4list)
      variableIndex = files2process4list{idatatype,1}(find(ismember(files2process4list{idatatype},pathtype),1,'last')+1:end-5);
      disp(['Processing file : ',files2process4list{idatatype}])
      try
          [variable,stocksUse,tUse,currencyUse] = readxlsxfiles_tds(files2process4list{idatatype},1,'timeseries',parameters.variables.string4missing);
      catch errorlog
          disp('ERROR reafing .xlsx file')
          disp(errorlog.message)
          continue
      end
      disp(['variable : ',variableIndex,' lines : ',num2str(size(variable,1)),' columns : ',num2str(size(variable,2))])
      disp(['stocksUse : ',variableIndex,' lines : ',num2str(length(stocksUse))])
      disp(['tUse : ',variableIndex,' lines : ',num2str(length(tUse))])
      eval([variableIndex,' = variable;'])
      save([savedirectory4list,variableIndex,'.mat'],variableIndex,'stocksUse','tUse','currencyUse','-v7.3')
      variablesProcessed4list  = vertcat(variablesProcessed4list,variableIndex);  %#ok<AGROW>
      eval(['clear ',variableIndex,' stocksUse tUse currencyUse variableIndex'])
      waitbar(idatatype/length(files2process4list))
  end
  clear idatatype
  close(h)
  clear h
  % Update Indexes for market by List 
  % check if all variables are processed 
  disp('---------------------')
  if length(files2process4list)<length(parameters.variables.timeseries)-1
      disp(['Missing Variables : ',parameters.variables.timeseries(~ismember(parameters.variables.timeseries,variablesProcessed4list))])
  else
      disp(['no missing variables for ',marketName])
  end
  disp('---------------------')
  disp(['process all markets progress : ',num2str(floor((iList/nList)*100)),'%'])
  disp('---------------------')
  clear variablesProcessed4list files4list fileName4static files2process4list marketName savedirectory4list
end
clear iList
% report problematic processed files
if sum(errorIndex)>0
    disp('-----------------------------------------------')
    disp(['# download errors : ',num2str(sum(errorIndex))])
    disp(['problematic lists : ',filesList(errorIndex)])
    disp('-----------------------------------------------')
else
    disp('-----------------------------------------------')
    disp('all files properly processed...')
    disp('-----------------------------------------------')
end
clear errorIndex nList tdsListNameUnique
%% produce data for each market / concatenate all lists / produce variables 
%  get a unique number of markets
parameters.tdsMarketList = unique(tdsMarketList); % keep the final list of markets to be processed by tds
% loop on tds markets 
for imarket = 1:length(parameters.tdsMarketList)
    %% prepare static data for market
    marketUse = parameters.tdsMarketList{imarket}; 
    list4market = tdsListIndex(ismember(tdsMarketListName,marketUse)); 
    listname4market = tdsListName(ismember(tdsMarketListName,marketUse)); 
    files4market = filesList(ismember(tdsMarketListName,marketUse));
    staticfiles4market = filesList(staticIndex&ismember(tdsMarketListName,marketUse)); 
    nfile2concatenate = length(unique(list4market)); 
    % loop on submarkets to concatenate
    for isubmarket =1:nfile2concatenate
        staticData = matfile(strrep([staticfiles4market{isubmarket}(1:end-5),'.mat'],'raw','temp')); % since the matfiles on the temporary data directory. 
        if isubmarket==1
              % identify static data directory 
              static = staticData.static; 
              headers = staticData.headers; 
              stocksUse_initial = staticData.stocksUse; 
              marketName = staticData.marketName; 
              % create directory for market
              marketDirectory = [parameters.dir.tmp.equities,lower(marketName),pathtype]; 
              if ~exist(marketDirectory,'dir'); mkdir(marketDirectory);end
        else
              static = vertcat(static,staticData.static);  %#ok<AGROW>
              if ~isequaln(headers,staticData.headers); error('Static Fields error');end
              stocksUse_initial = vertcat(stocksUse_initial,staticData.stocksUse);  %#ok<AGROW>
        end
    end
    % replace stocks with an "#NA" code in final variables. 
    isProblematicStock = ismember(stocksUse_initial,parameters.variables.string4missing); 
    stocksUse = stocksUse_initial(~isProblematicStock,1); 
    static = static(~isProblematicStock,:); 
    % save static data for market in the above directory with name as in country 
    save([marketDirectory,'static.mat'],'static','stocksUse','headers','marketName','stocksUse_initial','isProblematicStock','-v7.3'); 
    clear isubmarket static
    %% prepare time series data for all submarkets
    for idatatype = 1:length(parameters.variables.timeseries)
        variableIndex = parameters.variables.timeseries{idatatype};
        parameters.variables.currencyConversionIndex(idatatype,1) = ~ismember(variableIndex,parameters.variables.currencyNeutral); 
        datatypeIndex = cell2mat(cellfun(@(x) contains(x,[pathtype,variableIndex,'.xlsx']),files4market,'UniformOutput',false));
        files4datatype = files4market(datatypeIndex);    
        % loop on sublists (i.e. market having more than 1 list)
        for isubmarket =1:length(files4datatype) % This gives the number of submarkets:
            % need to define this in earlier parts of the code - at the exact same way.
            variableData = matfile(strrep([files4datatype{isubmarket}(1:end-5),'.mat'],'raw','temp')); % since the matfiles on the temporary data directory.
            if isubmarket==1
                % identify static data directory
                eval(['variable = variableData.',variableIndex,';']);
                stocksUse_variable = variableData.stocksUse;
                currencyUse_variable = variableData.currencyUse; 
            else
                eval(['variable = horzcat(variable,variableData.',variableIndex,');'])
                stocksUse_variable = vertcat(stocksUse_variable,variableData.stocksUse); %#ok<AGROW>
                currencyUse_variable = vertcat(currencyUse_variable,variableData.currencyUse); %#ok<AGROW>
            end
            tUse_variable = variableData.tUse; % Then Check if same to all variables. 
            clear variableData 
            % preallocate variable to check 
            if isubmarket==1
                tUseAll_submarkets= nan(length(tUse_variable),length(files4datatype)); 
            end
            tUseAll_submarkets(:,isubmarket) = tUse_variable; 
        end
        clear isubmarket
        % check lists have the same tUse
        if size(unique(tUseAll_submarkets','rows'),1)==1
            disp('tUse conformed for all lists in specific datatype and market')
            tUse = tUse_variable; 
            clear tUse_variable tUseAll_submarkets 
        else
            disp('tUse error among lists for market...exiting')
            continue 
        end
        % preallocate second step of checks 
        if idatatype==1
            tUse_byDatatype = nan(length(tUse),length(parameters.variables.timeseries)); 
        end
        tUse_byDatatype(:,idatatype) = tUse; 
        % relace stocks with an "#NA" Code in final variables
        isProblematicStock_variable = ismember(stocksUse_variable,parameters.variables.string4missing);
        if ~isequal(length(stocksUse_initial),length(stocksUse_variable))
            disp(['problematic data for market : ',marketUse,' for datatype : ',variableIndex,' skipping...'])
            continue
        else
            % Check that not nan stocks are the same 
            stocks2check_variable = stocksUse_variable(~isProblematicStock); 
            stocks2remove4check = isProblematicStock_variable(~isProblematicStock); 
            disp('surpassing extra checks for stocks')
            %{
            if ~isequal(stocksUse(~stocks2remove4check),stocks2check_variable(~stocks2remove4check))
                disp(['problematic stock conforming from lists for datatype and market : ',marketUse,' for datatype : ',variableIndex,' skipping...'])
                continue
            end
            %}
            clear stocks2check_variable stocks2remove4check
        end
        % conform data to static data
        eval([variableIndex,' = variable(:,~isProblematicStock);'])
        clear variable
        currencyUse = currencyUse_variable(~isProblematicStock,1); 
        currencyUse = cellfun(@(x) convert2cellstr(x,parameters.variables.string4missing),currencyUse,'UniformOutput',false);% Make Currency Use a cesstr variable
        stocksUse = stocksUse_variable(~isProblematicStock); 
        % save static data for market in the above directory with name as in country
        if parameters.variables.currencyConversionIndex(idatatype,1)
            save([marketDirectory,variableIndex,'_rawProcessed.mat'],variableIndex,'stocksUse','tUse','currencyUse');
        else
            % save raw variable with its final name and without currencyUse
            % variable. 
            save([marketDirectory,variableIndex,'.mat'],variableIndex,'stocksUse','tUse');
        end        
        eval(['clear ',variableIndex])   
        clear variableIndex dirs4datatype files4datatype datatypeIndex tUse ...
            stocksUse currencyUse currencyUse_variable...
            stocksUse_variable isProblematicStock_variable currencyConversionIndex
    end
    % check stocksUse and tUse among datatypes for the same market 
    if size(unique(tUse_byDatatype','rows'),1)==1
        disp('tUse conformed for all datatypes and lists in market')
        tUse = tUse_byDatatype(:,1); 
        clear tUse_byDatatype
    else
        disp('tUse error conforming among datatypes or lists for market...exiting')
        continue
    end
    %% translate all variables to home currency - use currency use for this - 
    % all are translated either they are filtered or not 
    % Make an analysis based on currencyUse and apply restrictions of
    % Landis and Skouras (2019). 
    disp('-----------------------')
    disp('Currency Conversions...')
    disp('-----------------------')
    exchangeRateFile = matfile([parameters.dir.processed.exchangeRates,'sterling_rates.mat']);
    euroRateFile = matfile([parameters.dir.processed.exchangeRates,'euro_rates.mat']);
    counter = 0 ; 
    for idatatype = 1:length(parameters.variables.timeseries)
        variableIndex = parameters.variables.timeseries{idatatype}; 
        % Correct for currencies or do nothing with the datatype in cases
        % no corrections are applied. 
        if parameters.variables.currencyConversionIndex(idatatype,1)
            % For the first datatype analyze currencies and exceptions
            % after accounting for missing data
            matfile2process = matfile([marketDirectory,variableIndex,'_rawProcessed.mat']); % get Currencies from unadjusted prices
            if counter==0
                counter = 1; 
                %read matfile with datatypes
                matfile2process4currency = matfile([marketDirectory,'up_rawProcessed.mat']); % get Currencies from unadjusted prices
                % Compute some stats on Currency Use
                currency2use = matfile2process4currency.currencyUse;
                % Get currencies
                currencyAnalysis = tabulate(currency2use);
                currencyAnalysis = currencyAnalysis(~ismember(currencyAnalysis(:,1),parameters.variables.string4missing),:);
                parameters.market.name{imarket,1} = lower(marketUse);
                homeCurrency = currencyAnalysis(cell2mat(currencyAnalysis(:,2))==max(cell2mat(currencyAnalysis(:,2))),1);
                homeCurrency = homeCurrency{1,1};
                parameters.market.homeCurrency{imarket,1} = homeCurrency;
                disp(['Home Currency : ',parameters.market.homeCurrency{imarket,1}])
                clear homeCurrency
                % define historical currency - for EMU in countries 
                if isequal(parameters.market.homeCurrency,'E')
                    currencyAnalysisNew = currencyAnalysis(~ismember(currencyAnalysis(:,1),'E'),:);
                    historicalCurrency = currencyAnalysisNew(cell2mat(currencyAnalysisNew(:,2))==max(cell2mat(currencyAnalysisNew(:,2))),1);
                    historicalCurrency = historicalCurrency{1,1};
                    parameters.market.historicalCurrency{imarket,1} = historicalCurrency;
                    parameters.market.isEmu(imarket,1) = true; 
                else
                    parameters.market.historicalCurrency{imarket,1} = ''; % empty for non emu in countries
                    parameters.market.isEmu(imarket,1) = false; 
                end
                disp(['Historical Currency (only for EMU) : ',parameters.market.historicalCurrency{imarket,1}])
                clear currencyAnalysisNew historicalCurrency 
                % Add secondary currency for some exceptions in markets
                currencyAnalysisFinal = currencyAnalysis(~ismember(currencyAnalysis(:,1),parameters.market.homeCurrency{imarket,1}),:);
                if ~isempty(currencyAnalysisFinal)
                    alternativeCurrency = currencyAnalysisFinal(cell2mat(currencyAnalysisFinal(:,2))==max(cell2mat(currencyAnalysisFinal(:,2))),1);
                    alternativeCurrency = alternativeCurrency{1,1};
                    parameters.market.alternativeCurrency{imarket,1} = alternativeCurrency;
                else
                    parameters.market.alternativeCurrency{imarket,1} = '';
                end
                disp(['Alternative Currency : ',parameters.market.alternativeCurrency{imarket,1}])
                clear alternativeCurrency currencyAnalysis 
            end
            % Then correct currencies based on the home currency
            eval(['variable=matfile2process.',variableIndex,';'])
            currencyVariable = matfile2process.currencyUse; 
            tUse = matfile2process.tUse; 
            stocksUse =matfile2process.stocksUse; 
            parameters2use.dir_exchangeRates = parameters.dir.processed.exchangeRates; 
            parameters2use.currencyConversionIndex = parameters.variables.currencyConversionIndex(idatatype,1);
            parameters2use.homeCurrency = parameters.market.homeCurrency; 
            parameters2use.historicalCurrency = parameters.market.historicalCurrency; 
            parameters2use.emuCurrencySet = parameters.filters.emuCurrencySymbols ; 
            parameters2use.isEmu = parameters.market.isEmu; 
            parameters2use.variableIndex = variableIndex; 
            parameters2use.string4missing = parameters.variables.string4missing; 
            parameters2use.priceVariables = parameters.variables.pricing; 
            [variableNew,currencyNew]=currencyCorrection(variable,currencyVariable,tUse,parameters2use);   
            % Check for Currency Correrctions - Need to adjust this code
            % also to new setup
            if isequal(setdiff(unique(currencyNew),{parameters2use.string4missing}),parameters.market.homeCurrency) || isequal(unique(currencyNew),{parameters2use.string4missing})
                disp('All Currencies Succesfully Corrected')
            else
                error('All Currencies are not properly corrected')
            end
            % dynamically rename variable 
            eval([variableIndex,'=variableNew;'])
            currencyUse=currencyNew;
            disp(['Saving Raw Conformed Variable : ',variableIndex])
            save([marketDirectory,variableIndex,'.mat'],...
                variableIndex,'currencyUse','stocksUse','tUse','-v7.3')
            % Compute be variable
            if isequal(variableIndex,'WC03501')
                disp('')
                disp('Produce Book value variable...')
                be = WC03501; 
                save([marketDirectory,'be.mat'],...
                     'be','currencyUse','stocksUse','tUse','-v7.3') 
                clear be WC03501
            end
            % Compute return variable 
            if isequal(variableIndex,'ri')
                disp('')
                disp('Produce returns variable...')
                ret = [nan(1,length(stocksUse));-1+(ri(2:end,:)./ri(1:end-1,:))]; 
                save([marketDirectory,'ret.mat'],...
                     'ret','currencyUse','stocksUse','tUse','-v7.3') 
                clear ret
            end
            clear variable variableNew currencyVariable stocksUse tUse parameters2use
        else
            % do nothing
        end
        % Note - We do not remove any stocks at this point. All stocks will
        % be removed after the application of all filters in stocks. 
        clear variableIndex matfile2process 
    end
    clear counter
    %% create new datatypes/variables which will be used in analysis
    % load static data
    disp('Create datatypes from static data information....')
    static2use = load([marketDirectory,'static.mat']); 
    tUse = load([marketDirectory,'ri.mat'],'tUse'); tUse=tUse.tUse; 
    % get indexes of fields
    [~,iyTRAC] =ismember('TRAC',static2use.headers);
    [~,iyTRAD] =ismember('TRAD',static2use.headers);
    [~,iyENAME] =ismember('ENAME',static2use.headers);
    [~,iyPCUR] =ismember('PCUR',static2use.headers);
    [~,iyGEOGN] =ismember('GEOGN',static2use.headers);
    [~,iyEXMNEM] =ismember('EXMNEM',static2use.headers);
    [~,iyLOC] =ismember('LOC',static2use.headers);
    [~,iyISIN] =ismember('ISIN',static2use.headers);
    [~,iyISINID] =ismember('ISINID',static2use.headers);
    % define variables by converting into cellstrings
    trac = cellfun(@(x) convert2cellstr(x,parameters.variables.string4missing),static2use.static(:,iyTRAC),'UniformOutput',false); clear iyTRAC
    trad = cellfun(@(x) convert2cellstr(x,parameters.variables.string4missing),static2use.static(:,iyTRAD),'UniformOutput',false); clear iyTRAD
    ename = cellfun(@(x) convert2cellstr(x,parameters.variables.string4missing),static2use.static(:,iyENAME),'UniformOutput',false); clear iyENAME
    pcur = cellfun(@(x) convert2cellstr(x,parameters.variables.string4missing),static2use.static(:,iyPCUR),'UniformOutput',false); clear iyPCUR
    geogn = cellfun(@(x) convert2cellstr(x,parameters.variables.string4missing),static2use.static(:,iyGEOGN),'UniformOutput',false); clear iyGEOGN
    exmnem = cellfun(@(x) convert2cellstr(x,parameters.variables.string4missing),static2use.static(:,iyEXMNEM),'UniformOutput',false); clear iyEXMNEM
    loc = cellfun(@(x) convert2cellstr(x,parameters.variables.string4missing),static2use.static(:,iyLOC),'UniformOutput',false); clear iyLOC
    isin = cellfun(@(x) convert2cellstr(x,parameters.variables.string4missing),static2use.static(:,iyISIN),'UniformOutput',false); clear iyISIN
    isinid = cellfun(@(x) convert2cellstr(x,parameters.variables.string4missing),static2use.static(:,iyISINID),'UniformOutput',false); clear iyISINID
    stocksUse = static2use.stocksUse;
    % create exchangeCode variable
    allVenuesInfo = tabulate(exmnem); 
    allVenuesInfo = allVenuesInfo(~ismember(allVenuesInfo(:,1),parameters.variables.string4missing),:);
    % isXetra = (isequal(marketUse,'Germany') && isequal(allVenuesInfo{1,1},'XETRA')); 
    [~,locInInfo] = ismember(exmnem,allVenuesInfo(:,1)); 
    exchangeCode = repmat(locInInfo',length(tUse),1);
    parameters.markets.exchanges{imarket,1} = allVenuesInfo(:,1); % Keep the exchange classifications here 
    isMainVenue = exchangeCode==1; 
    clear allVenuesInfo locInInfo exmnem
    % create home currency variable - asset trades in home currency
    allCurrencyInfo = tabulate(pcur); 
    allCurrencyInfo = allCurrencyInfo(~ismember(allCurrencyInfo(:,1),parameters.variables.string4missing),:);
    [~,locInInfo] = ismember(pcur,allCurrencyInfo(:,1)); 
    currencyCode = repmat(locInInfo',length(tUse),1);
    parameters.markets.currencies{imarket,1} = allCurrencyInfo(:,1); % Keep the exchange classifications here
    % define currencies to use 
    if parameters.market.isEmu || ismember(marketUse,parameters.filters.currencyExceptions)
        isHomeCurrency = currencyCode==1||currencyCode==2;
    else
        isHomeCurrency = currencyCode==1;
    end
    clear allCurrencyInfo locInInfo pcur
    % create share Type variable - Actually the opposite we used to do in
    % previous version
    shareType = ismember(trac,vertcat(parameters.filters.commonTracSet,parameters.variables.string4missing));    
    shareType = repmat(shareType',length(tUse),1); 
    isCommon_trac = shareType == 1; 
    clear trac trad 
    % create GEOGN variable
    isDomestic = ismember(geogn,marketUse);
    isDomestic = repmat(isDomestic',length(tUse),1); 
    clear geogn
    % create isin variable - Keep it as it is 
    % isin = repmat(isin',length(tUse),1); 
    % create major variable just in case we need this
    isMajor = ismember(isinid,'P'); 
    isMajor = repmat(isMajor',length(tUse),1); 
    % process name variable - define non common based on name and
    % crossListings - Here load this as Table and not csv -  
    abbList = table2array(readtable('abbreviations_list.csv')); 
    abbs.markets = abbList(:,1); 
    abbs.abbs = abbList(:,2); 
    abbs.exceptions = abbList(:,3); 
    abbs.crossList = abbList(:,4); 
    % get the data for the specific market
    [ix,iy] = ismember(marketUse,upper(abbs.markets));
    if ix>0
        exceptions = abbs.exceptions(iy);
        crossListings = abbs.crossList(iy);
        abbs = abbs.abbs(iy);
    else
        disp(['Market : ',market2check,' Not in Abbreaviations List... skipping this market'])
        continue
    end
    clear ix iy
    % define variable based on its ename abbreaviations 
    isNonCommon_abbs = false(length(stocksUse),1);
    if ~isempty(abbs{1,1})
        blacklists = textscan(abbs{1,1},'%q',200,'delimiter',',','whitespace',''); % This will Keep exact whitespaces
        blacklists = blacklists{:};
        if ~isempty(exceptions{:})
             exceptions = textscan(exceptions{:},'%q',20,'delimiter',','); % Leave it simple for the moment 
        end
        exceptions = exceptions{:};
        % Remove Stocks Containing any of the Abbreviations
        abbreviationUsed = cell(length(stocksUse),1);
        abbEffectCounter = zeros(length(blacklists),1);
        for iStock = 1:length(stocksUse) 
            for iAbb = 1:length(blacklists)
                isInEname = strfind(ename{iStock},blacklists{iAbb});
                if ~isempty(isInEname)
                    % Check if an Exception
                    if ~isempty(exceptions)
                        exceptionCounter =false(length(exceptions),1);
                        for iException=1:length(exceptions)
                            isException = strfind(ename{iStock},exceptions{iException});
                            if ~isempty(isException)
                                exceptionCounter(iException)=true;
                            end
                        end
                        % No Exceptions were found
                        if sum(exceptionCounter)==0 
                            abbEffectCounter(iAbb) = abbEffectCounter(iAbb)+1;
                            isNonCommon_abbs(iStock,1) = true;
                            abbreviationUsed{iStock}=blacklists{iAbb};
                        end       
                    else    
                        abbEffectCounter(iAbb) = abbEffectCounter(iAbb)+1;
                        isNonCommon_abbs(iStock,1) = true;
                        abbreviationUsed{iStock}=blacklists{iAbb};
                    end  
                end
                clear isInEname
            end
        end
    else % For the cases no Abbreviations are available
         blacklists =cell(1,1);
         abbEffectCounter = [];
    end
    clear iAbb istock abbEffectCounter
    isNonCommon_abbs = repmat(isNonCommon_abbs',length(tUse),1); 
    % Get the cross listings 
    isRemoved_crossList = false(length(stocksUse),1);
    if ~isnan(crossListings{:})
        crossListings = textscan(crossListings{:},'%q',200,'delimiter',',');
        crossListings = crossListings{:};
        % Remove Stocks Containing any of the Abbreviations
        for iStock = 1:length(stocksUse)
            for iAbb = 1:length(crossListings)
                isCrossListing = strfind(ename{iStock},crossListings{iAbb});
                if ~isempty(isCrossListing)   
                    if ~isequal(marketUse,'Germany') || isxetra
                        isRemoved_crossList(iStock) = true;
                    end
                end
                clear isCrossListing
            end
        end
        clear iAbb iStock
    end
    clear crossListings
    % save all variables 
    disp('saving static data based variables...')
    save([marketDirectory,'exchangeCode.mat'],...
        'exchangeCode','isxetra','isRemoved_crossList','stocksUse','tUse','isMainVenue','-v7.3')
    clear exchangeCode isxetra isRemoved_crossList isMainVenue
    save([marketDirectory,'shareType.mat'],...
        'shareType','isCommon_trac','isNonCommon_abbs','stocksUse','tUse','-v7.3')   % all common stock based variables
    clear shareType isCommon_trac isNonCommon_abbs
    save([marketDirectory,'currencyCode.mat'],...
        'currencyCode','isHomeCurrency','stocksUse','tUse','-v7.3')
    clear currencyCode isHomeCurrency 
    save([marketDirectory,'isin.mat'],...
        'isin','stocksUse','tUse','-v7.3')
    clear isin 
    save([marketDirectory,'major.mat'],...
        'isMajor','stocksUse','tUse','-v7.3')
    clear isMajor
    save([marketDirectory,'geogn.mat'],...
        'isDomestic','stocksUse','tUse','-v7.3')
    clear isDomestic isin* loc
    %% Filter Data for each market 
    disp('-----------------------')
    disp('...filtering data......')
    disp('-----------------------')
    % define filtered data directory 
    saveDirectory = [parameters.dir.processed.equities,lower(marketUse),pathtype];
    if ~exist(saveDirectory,'dir')
        mkdir(saveDirectory)
    end
    % produce filtered returns 
    load([marketDirectory,'ri.mat'],'ri','stocksUse','tUse')
    if isequal(parameters.filters.filteringBase,'ret')
        [ret,mv,tUse,stocksUse]=filters_tds(parameters,marketDirectory,'ret');
        [ri,~,~,~]=filters_tds(parameters,marketDirectory,'ri');
    else
        [ri,mv,tUse,stocksUse]=filters_tds(parameters,marketDirectory,'ri');
        ret = [nan(1,length(stocksUse));-1+(ri(2:end,:)./ri(1:end-1,:))];
    end
    % discard stocks after filtering 
    isAllNanStocks=mean(isnan(ret))==1;
    isAllNanDays=mean(isnan(ret),2)==1;   
    % Remove Completely Empty Entries :
    ret=ret(~isAllNanDays,~isAllNanStocks);
    ri=ri(~isAllNanDays,~isAllNanStocks);
    stocksUse=stocksUse(~isAllNanStocks);
    tUse=tUse(~isAllNanDays);
    save([saveDirectory,'ret.mat'],'ret','stocksUse','tUse','-v7.3')
    save([saveDirectory,'ri.mat'],'ri','stocksUse','tUse','isAllNanStocks','isAllNanDays','-v7.3')
    clear ret ri
    % apply correction to marketcap also 
    mv=mv(~isAllNanDays,~isAllNanStocks);
    save([saveDirectory,'mv.mat'],'mv','stocksUse','tUse','-v7.3')
    clear mv
    % filter and align data for all other datatypes except ri and mv 
    datatypes2filter = setxor(vertcat(parameters.variables.timeseries,{'static';'be'}),{'ri','mv','WC03501'});
    for idatatype = 1:length(datatypes2filter)
        variableIndex = datatypes2filter{idatatype};
        if ~isequal(variableIndex,'static')
            [isPricing,ixLocPricing] = ismember(variableIndex,parameters.variables.pricing.adjusted);
            % if variable is pricing recreate its adjusted data and apply filters
            if isPricing
                af = matfile([marketDirectory,'af.mat']);
                unadjustedVariable = matfile([marketDirectory,parameters.variables.pricing.adjusted{ixLocPricing},'.mat']);
                eval(['uvariable = unadjustedVariable.',parameters.variables.pricing.adjusted{ixLocPricing},';']);
                adjustedVariable = af.af.*uvariable;
                eval([variableIndex,'= adjustedVariable;']);
                clear uvariable unadjustedVariable adjustedVariable
                save([saveDirectory,variableIndex,'.mat'],variableIndex,'stocksUse','tUse','-v7.3')
                clear af
            end
            clear isPricing ixLocPricing
            % filter some variables of interest or simply align data
            if isequal(variableIndex,'p')
                [p,~,~,~]=filters_tds(parameters,marketDirectory,'p');
                % create other variables
                retx=[nan(1,size(p,2));-1+p(2:end,:)./p(1:end-1,:)];
                retx=retx(~isAllNanDays,~isAllNanStocks);
                save([saveDirectory,'retx.mat'],'retx','stocksUse','tUse','-v7.3')
                p=p(~isAllNanDays,~isAllNanStocks);
                save([saveDirectory,'p.mat'],'p','stocksUse','tUse','-v7.3')
                clear retx p
            elseif isequal(variableIndex,'po')
                [po,~,~,~]=filters_tds(parameters,marketDirectory,'po');
                po=po(~isAllNanDays,~isAllNanStocks);
                save([saveDirectory,'po.mat'],'p','stocksUse','tUse','-v7.3')
                clear po
            else
                % load matfile
                variable = load([marketDirectory,variableIndex,'.mat']);
                % align to returns variable
                eval([variableIndex,'=variable.',variableIndex,'(~isAllNanDays,~isAllNanStocks);'])
                % save variable
                save([saveDirectory,variableIndex,'.mat'],variableIndex,'stocksUse','tUse','-v7.3')
                eval(['clear ',variableIndex])
            end
        else
            variable = load([marketDirectory,variableIndex,'.mat']);
            static = variable.static(~isAllNanStocks,:);
            headers = variable.headers;
            save([saveDirectory,variableIndex,'.mat'],variableIndex,'stocksUse','headers','-v7.3')
            clear static headers
        end
        clear variable
    end
    % clear all leftover variables
    clear idatatype isAllDaysNan isAllStockNan files4market filesList headers...
         staticfiles4market marketUse marketName...
         stocksUse_initial list4market listname4market marketDirectory...
         isProblematicStock variableIndex         
end
clear imarket nfile2concatenate tdsMarketList 
%% apply other actions 
% handle filtered data for german market 

%% remove the temporary data directories and delete temp data
rmdir(parameters.dir.tmp.equities,'s')

%%%%%%%%%%%%%%%%%%%%%%%%
%% Embeeded Functions %%
%%%%%%%%%%%%%%%%%%%%%%%%
%% read excel files from tds 
% works if user has followed the steps provided in the online appendix 
function [data,nargout1,nargout2,nargout3] = readxlsxfiles_tds(fileName,sheetName,requestType,string4missing)
if nargout<4
    string4missing = '#NA'; 
end
[~,~,raw]=xlsread(fileName,sheetName);
if isequal(requestType,'timeseries')
    % faster to process by stock
    data = cell2mat(cellfun(@(x) convert_data(x),raw(4:end,2:end),'UniformOutput',false)); % Change this from raw data by transforming "#NA" to nan. - this was wrongly applied here -
    clear istock raw2process
    nargout1 = cellfun(@(x) convert_stockName(x,string4missing),raw(2,2:end),'UniformOutput',false);  % Stocks as Cell
    nargout1 = nargout1'; % get stocks by removing the formula 
    nargout2 = datenum(raw(4:end,1)); % Convert dates to datenumbers
    nargout3 = raw(3,2:end)'; % Currencies as Cell
elseif isequal(requestType,'static')
    if nargout>2
        warning('Too many outputs requested for static files')
        nargout2 = []; 
        nargout3 = []; 
    end
    % Check for empty fields
    ixEndFields = find(cell2mat(cellfun(@(x) ~ischar(x),raw(1,:),'UniformOutput',false))~=1,1,'last'); 
    data = raw(2:end,1:ixEndFields);
    nargout1 = raw(1,1:ixEndFields); 
end
end
%% convert stock names
function [stockName] = convert_stockName(stockNameInFunction,string4missing)
ixEndFunction = max(strfind(stockNameInFunction,'(')); 
if ~isempty(ixEndFunction)
    stockName = stockNameInFunction(6:ixEndFunction-1);
else
    stockName = string4missing; 
end
end
%% convert data in cells to doubles 
function [dataAsNumber] = convert_data(dataAsCellentry)
if ~isnumeric(dataAsCellentry) || isempty(dataAsCellentry)
    dataAsNumber = nan; 
elseif isnumeric(dataAsCellentry) && ~isempty(dataAsCellentry)
    dataAsNumber = dataAsCellentry; 
end
end
%% convert nan to strings 
function [cellAscellstr] = convert2cellstr(cellVariable,missingString)
     if isnan(cellVariable)
         cellAscellstr = missingString; 
     else
         cellAscellstr=cellVariable;
     end
end
%% get filenames of files to specific dir
function [fileList] = getfilenames(directory,filetype)
% [fileList] = getfilenames(directory,filetype);
% readFileNames function returns all directories for the files of specific
% filetype in a directory and all of its subdirectories. 
%
% example : 
% 
% [matfileNames] = readFileNames(directory,'.mat');
DirContents=dir(directory);
fileList=[];
extList=filetype;
if ispc % or (strcmpi(computer,'PCWIN') || strcmpi(computer,'PCWIN64'))
    NameSeperator='\';
elseif isunix % or (strcmpi(computer,'GLNX86') || strcmpi(computer,'GLNXA64'))
    NameSeperator='/';
end
% dir contents of folders 
for i=1:numel(DirContents)
    if(~(strcmpi(DirContents(i).name,'.') || strcmpi(DirContents(i).name,'..')))
        if(~DirContents(i).isdir)
            extension=DirContents(i).name(end-length(filetype)+1:end);
            if(numel(find(strcmpi(extension,extList)))~=0)
                fileList=cat(1,fileList,{[directory,NameSeperator,DirContents(i).name]});
            end
        else
            getlist=getfilenames([directory,NameSeperator,DirContents(i).name],filetype);
            fileList=cat(1,fileList,getlist);
        end
    end
end
end
%% get index of files in the directory based on market numbering
function [folderIndex,folderName] = getfolderindex(filename,pathtype)
% get the index from folder directory
    ixpathtype = strfind(filename,pathtype);
    folderName = filename(ixpathtype(end-1)+1:ixpathtype(end)-1);
    folderIndex = str2double(folderName);
end
%% Currency Correction 
function [dataInHomeCurrency,currencyVector]= currencyCorrection(dataRaw,currencyVectorRaw,dateVector,parameters2use)
% Apllies to data in market expressed in a currency other than the home
% currency 
% preallocate output variables 
dataInHomeCurrency = dataRaw;
currencyVector = currencyVectorRaw;
variableIndex = parameters2use.variableIndex; 
disp('---------------------------')
disp('Currency Conversion Report')
disp('---------------------------')
tabulate(currencyVectorRaw)
disp('---------------------------')
% Do something only if variable needs conversion by default 
if parameters2use.currencyConversionIndex
    % load EMU in rates
    if parameters2use.isEmu
        eRate = load([parameters2use.dir_exchangeRates,'euro_rates.mat']);
        emuConversionRatesNames = eRate.stocksUse;
        emuConversionRates = eRate.ret; 
        emuConversionRatesCurrency = eRate.currencyUse; 
        emuConversionRatesDates = eRate.tUse; 
        emuConversionRatesStocks = eRate.stocksUse; clear eRate %#ok<NASGU>
        [isHomeCountry,targetHomeColn] = ismember(parameters2use.historicalCurrency,emuConversionRatesCurrency);
        disp(['Converting From Currency : ',emuConversionRatesNames{1,targetHomeColn}])
        if ~isHomeCountry
            disp('Home Currency Conversion Rates could not be found')
        end
        homeRate = emuConversionRates(:,targetHomeColn);
        [homeConversionRate,homeDatesConformed] = singledateCut(homeRate,emuConversionRatesDates,dateVector);
        if isequal(homeDatesConformed,dateVector)&& isequal(size(homeConversionRate,1),size(dataRaw,1))
            disp('Conformation Succesfull')
        else
            disp('Error in Rates Dates Conformation')
        end
    end
    % Procced to Currency Correction
   if parameters2use.isEmu
        disp(' ')
        disp([ '# No Home Currency Stocks: ', num2str(size(currencyVectorRaw,1)- size(strmatch(parameters2use.homeCurrency,currencyVectorRaw,'exact'),1)- size(strmatch(parameters2use.historicalCurrency,currencyVectorRaw,'exact'),1)-size(strmatch({''},currencyVectorRaw,'exact'),1)) ])
        for i=1:size(currencyVectorRaw,1)
            if~isequal(currencyVectorRaw(i,1),{parameters2use.string4missing})
                if isequal(currencyVectorRaw(i,1),parameters2use.homeCurrency)
                    dataInHomeCurrency(:,i)= dataRaw(:,i);% or element by element multiplication with ones
                elseif isequal(currencyVectorRaw(i,1),parameters2use.historicalCurrency)
                    dataInHomeCurrency(:,i)= dataRaw(:,i)./ homeConversionRate; %Here put the fixed ireversible rate to euro.
                    currencyVector(i,1) = parameters2use.homeCurrency;
                elseif isequal(isequal(currencyVectorRaw(i,1),parameters2use.homeCurrency),0) && isequal(isequal(currencyVectorRaw(i,1),parameters2use.historicalCurrency),0)
                    fromCurrencyIndex=currencyVectorRaw(i,1);
                    toCurrencyIndex = parameters2use.historicalCurrency; % For Possible options type help currencyConverter
                    dataInHomeCurrency(:,i)= currencyConverter(dataRaw(:,i),dateVector',fromCurrencyIndex,toCurrencyIndex,variableIndex,parameters2use);
                    currencyVector(i,1) = parameters2use.homeCurrency;
                end
            end
        end
        disp([ '# of stocks with All days NaN:', num2str(sum(mean(isnan(dataInHomeCurrency))==1))  ])
        disp([ 'No Home Currency Indicators:', num2str(size(currencyVectorRaw,1)- size(strmatch(parameters2use.homeCurrency,currencyVectorRaw,'exact'),1)- size(strmatch(parameters2use.historicalCurrency,currencyVectorRaw,'exact'),1)-size(strmatch({''},currencyVectorRaw,'exact'),1))])
        disp([ 'No Home Currency Corrected Indicators:', num2str(size(currencyVector,1)- size(strmatch(parameters2use.homeCurrency,currencyVector,'exact'),1)- size(strmatch(parameters2use.historicalCurrency,currencyVector,'exact'),1)-size(strmatch({''},currencyVector,'exact'),1))])
        disp(' ')      
    elseif ~parameters2use.isEmu
        disp(' ')
        disp([ '# No Home Currency Stocks: ', num2str(size(currencyVectorRaw,1)- size(strmatch(parameters2use.homeCurrency,currencyVectorRaw,'exact'),1)-size(strmatch({''},currencyVectorRaw,'exact'),1))])
        index =zeros(size(currencyVectorRaw,1),1);
        for i=1:size(currencyVectorRaw,1)
            if~isequal(currencyVectorRaw(i,1),{parameters2use.string4missing})
                if isequal(currencyVectorRaw(i,1),parameters2use.homeCurrency)
                    dataInHomeCurrency(:,i)= dataRaw(:,i);% or element by element multiplication with ones
                elseif ~isequal(currencyVectorRaw(i,1),parameters2use.homeCurrency)
                    index(i)=1; 
                    fromCurrencyIndex=currencyVectorRaw(i,1);
                    toCurrencyIndex = parameters2use.homeCurrency; % For Possible options type help currencyConverter
                    dataInHomeCurrency(:,i)= currencyConverter(dataRaw(:,i),dateVector',fromCurrencyIndex,toCurrencyIndex,variableIndex,parameters2use);
                    currencyVector(i,1) = parameters2use.homeCurrency;
                    disp(['Correction non home : ', currencyVectorRaw{i,1},' vs : ',parameters2use.homeCurrency{1,1}])
                end
            end
        end
        disp([ '# of stocks with All days NaN:', num2str(sum(mean(isnan(dataInHomeCurrency))==1))  ])
        disp([ 'No Home Currency Indicators:', num2str(size(currencyVectorRaw,1)- size(strmatch(parameters2use.homeCurrency,currencyVectorRaw,'exact'),1)-size(strmatch({''},currencyVectorRaw,'exact'),1))])
        disp([ 'No Home Currency Corrected Indicators:', num2str(size(currencyVector,1)- size(strmatch(parameters2use.homeCurrency,currencyVector,'exact'),1)-size(strmatch({''},currencyVector,'exact'),1))])
        disp(' ')   
   end
else
    % For Ratios Only transform all Currencies to Identifier.
    for i=1:size(currencyVectorRaw,1)
        if~isequal(currencyVectorRaw(i,1),{parameters2use.string4missing})
            if ~isequal(currencyVectorRaw(i,1),parameters2use.homeCurrency)
                currencyVector(i,1) = parameters2use.homeCurrency;
            end
        end
    end
end
disp('')
end
%% Currency Converter
function [xInToCurrency,conversionRate]=currencyConverter(xInFromCurrency,t,fromCurrencyIndex,toCurrencyIndex,nameOfX,parameters2use)
% Fuction to convert Series to all possible DTS currencies.
% 
% xInFromCurrency: data to be multiplied by appropriate exchange rate
% nameOfX: name of xInFromCurrency, i.e. 'p' if p, which is used to handle UK
% case which is different due to penny denomination
% {from,to}CurrencyIndex: use codes from list that follows
%
% isEmu: true/false scalar stating whether the fromCurrencyIndex refers to a country that has joined EMU 
% In such cases, the fromCurrencyIndex should be equal to MARKETPARAMETERS.precurIdentifier;
% otherwise it should be fromCurrencyIndex=MARKETPARAMETERS.curIdentifier;
%
% For all cases the Thompson Reuters daily Closing Spot rates on Sterling
% are used. So all cross rates are computed via UK Pound rate. 
%
% Example:
%
% toCurrencyIndex = {'U$'};
%
% [upDollar] = currencyConverter(up,t,fromCurrencyIndex,toCurrencyIndex,isEmu,'up');
% 
%
% Important Notes : 
%
% Requires PARAMETERS files for each market in workspace (MARKETPARAMETERS
% is used as global. Also user should follow our TDS directions -in docs- for encoding
% as some currency symbols to be handled (this will be removed in a future release.)
%
% For the case of UK price TDS data, as they are expressed in UK
% pennies we divide series with 100 before converting.
% Datatypes that need this conversion are defined is 
% VARIABLES.priceVariables variable. 
%
% If Currency is not identifiable (currency abbreaviation not matching with
% list of currencies symbols in currency files) code sets unkown currency stock's data to
% nan. So check if Datastream re-codes currencies abbreaviations. 
%
% The Cayman Islands rate is missing from Reuters UK Pound series, so we
% use instead the US to CD and US to UK rate to compute it. 
% 3 similar cases : BELIZE, BRUNEI and CFA FRANC BENIN are ignored. 
%
% Check VB currency if is an abbreviation change.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Landis Conrad, 2016. 
% /academicresearch/filtering/ 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Exchange Rates Data - And the appropriate Rates for EMU in
% countries ( EMU12 or EMU17)
xx
disp(['Converting From : ',fromCurrencyIndex{1,1},' to : ',toCurrencyIndex{1,1}])
disp(' ')
emuCurrencySet = parameters2use.emuCurrencySet;
isEmuHome = ismember(fromCurrencyIndex,emuCurrencySet);
isEmuBase = ismember(toCurrencyIndex,emuCurrencySet);
ratesBase = load([parameters2use.dir_exchangeRates,'sterling_rates.mat']);
rates2ukDates=ratesBase.tUse;rates2uk=ratesBase.ret;rates2ukNames=ratesBase.stocksUse;rates2ukCurrency=ratesBase.currencyUse; clear ratesBase
% this loads rates2uk, rates2ukCurrency, rates2ukDates, rates2ukNames,
% rates2ukStocks variables
rates2ukDates=rates2ukDates'; % check if this is still needed to be transposed
priceVariables = parameters2use.priceVariables;
if isEmuHome
    % check what is this and define it in code
    rates=load([ratesDirectory,'uk2emuRates.mat']);uk2emuRatesDates=rates.uk2emuRatesDates;uk2emuRates=rates.uk2emuRates;uk2emuRatesNames=rates.uk2emuRatesNames;uk2emuRatesCurrency=rates.uk2emuRatesCurrency;clear rates
    uk2emuRatesDates=uk2emuRatesDates';
    % This loads uk2emuRate, uk2emuRatesCurrency, uk2emuRateDates, uk2emuRatesNames,
    % uk2emuRateStocks variables. Define currency to convert from and confrom rate on xInFromCurrency's dates 
    [isHomeCountry,targetHomeColn]=ismember(fromCurrencyIndex,uk2emuRatesCurrency);
    if isHomeCountry
        homeRate = uk2emuRates(:,targetHomeColn);
        disp(['From Currency Rate Applied : ',uk2emuRatesNames{1,targetHomeColn}])
        [homeRateConformed,homeDatesConformed] = singledateCut(homeRate,uk2emuRatesDates,t);
    else
        disp('Currency To Convert From does not Exists in Rates Files')
        disp('Removing... Stock Obs')
        homeRateConformed = nan(size(xInFromCurrency));
        homeDatesConformed = t;
    end
else
    [isHomeCountry,targetHomeColn]=ismember(fromCurrencyIndex,rates2ukCurrency);
    if isequal(fromCurrencyIndex(1,1),{char(163)})
        isHomeCountry=true;
    end
    if isHomeCountry
        if isequal(fromCurrencyIndex(1,1),{char(163)})
            homeRate = ones(size(rates2uk,1),1);
            disp(['From Currency Rate Applied : ',char(163)])
        else
            homeRate = rates2uk(:,targetHomeColn);
            disp(['From Currency Rate Applied : ',rates2ukNames{1,targetHomeColn}])
        end
        [homeRateConformed,homeDatesConformed] = singledateCut(homeRate,rates2ukDates,t);
    else
        disp('Currency To Convert From does not Exists in Rates Files')
        disp('Removing... Stock Obs')
        homeRateConformed = nan(size(xInFromCurrency));
        homeDatesConformed = t;
    end
    % handle UK pennies - this should apply only to currency data which are expressed in pennies
    if isequal(fromCurrencyIndex(1,1),{char(163)}) && ismember(nameOfX,priceVariables)
        xInFromCurrency = xInFromCurrency./100;
    end
end
% define currency to convert to and confrom rate on xInFromCurrency's dates
if isEmuBase
    rates=load([ratesDirectory,'uk2emuRates.mat']);uk2emuRatesDates=rates.uk2emuRatesDates;uk2emuRates=rates.uk2emuRates;uk2emuRatesNames=rates.uk2emuRatesNames;uk2emuRatesCurrency=rates.uk2emuRatesCurrency;clear rates
    uk2emuRatesDates=uk2emuRatesDates';
    [isBaseCountry,targetBaseColn]=ismember(toCurrencyIndex,uk2emuRatesCurrency);
    if isBaseCountry
        baseRate = uk2emuRates(:,targetBaseColn);
        disp(['To Currency Rate Applied : ',uk2emuRatesNames{1,targetBaseColn}])
        [baseRateConformed,baseDatesConformed] = singledateCut(baseRate,uk2emuRatesDates,t);
    else
        disp('Currency To Convert To does not Exists in Rates Files')
        disp('Removing... Stock Obs')
        baseRateConformed = nan(size(xInFromCurrency));
        baseDatesConformed = t;
    end
else
    [isBaseCountry,targetBaseColn]= ismember(toCurrencyIndex,rates2ukCurrency);
    if isequal(toCurrencyIndex(1,1),{char(163)})
        isBaseCountry=true;
    end
    if isBaseCountry
        if isequal(toCurrencyIndex(1,1),{char(163)})
            baseRate = ones(size(rates2uk,1),1);
            disp(['To Currency Rate Applied : ',char(163)])
        else
            baseRate = rates2uk(:,targetBaseColn);
            disp(['To Currency Rate Applied (inversely) : ',rates2ukNames{1,targetBaseColn}])
        end
        [baseRateConformed,baseDatesConformed] = singledateCut(baseRate,rates2ukDates,t);
    else
        disp('Currency To Convert To does not Exists in Rates Files')
        disp('Removing... Stock Obs')
        baseRateConformed = nan(size(xInFromCurrency));
        baseDatesConformed = t;
    end
    if isequal(toCurrencyIndex(1,1),{char(163)}) && ismember(nameOfX,priceVariables)
        xInFromCurrency = xInFromCurrency.*100; % So pricing data of another currency to be transformed to Pennies
    end
end
% check conformation success on both sides
if isequal(homeDatesConformed,baseDatesConformed,t)&& isequal(size(homeRateConformed,1),size(baseRateConformed,1),size(xInFromCurrency,1))
    disp('Conformation Successful')
else
    disp('Error in rates Conformation')
    disp('Exiting Conversion...')
    xInToCurrency =[];
    return
end
% convert to base currency - here again we have to evaluate EMU differently
if isEmuHome % emu via rates are inverted compared to other rates against 
    homeRateConformed = (1./homeRateConformed);
end
if isEmuBase 
    baseRateConformed = (1./baseRateConformed);
end
conversionRate = ((1./homeRateConformed).*baseRateConformed);
xInToCurrency = xInFromCurrency.*repmat(conversionRate,1,size(xInFromCurrency,2));
end
%% Cut to subsample dates 
function [xNew,tNew] = singledateCut(x,tInitialSeries,tLowerFrequency)
% This function conforms a daily frequency serie to a lower frequency
% serie, by providing the date vectors. x and tUsed are the initial
% variables-dates
isNewDate = ismember(tInitialSeries,tLowerFrequency);
xNew = x(isNewDate,:);
tNew =tInitialSeries(isNewDate);
% check consistency
if isequaln(tLowerFrequency,tNew)
    disp('Dates Correctly Conformed')
end
end
%% Fix for German data and XETRA 


%% Landis and Skouras (2019), "A granular approach to international equity data from Thomson Datastream". 





