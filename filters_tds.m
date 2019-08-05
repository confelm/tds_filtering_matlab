function [variableClean,mvClean,tUse,stocksUse,ixDateUse,ixStocksUse]=filters_tds(parameters,loadDirectory,variableIndex,filterIndex)
% Apply ls daily filters 

% define default values of input arguments 
if nargin<=3
    filterIndex = ones(16,1);
end
% load data and get type of filtering 
disp(['applying filters in variable : ',variableIndex])
% Load Data - prices or return indexes 
variable = load([loadDirectory,variableIndex,'.mat']);
stocksUse = variable.stocksUse; 
tUse = variable.tUse; 
eval(['variable = variable.',variableIndex,';']); 
% Define Filtering Variables
nDays = size(variable,1); 
kStocks = size(variable,2); 
isExclude = false(size(variable));
isVariableRet = isequal(variableIndex,'ret');
isVariableRI = isequal(variableIndex,'ri');

% Inform User that a non closing or non openning pricing datatype will be
% handled as a closing price. 
if ~isVariableRet && ~isVariableRI && ~isequal(variableIndex,'p') && ~isequal(variableIndex,'po')
    disp(' ')
    disp('variable will be filtered as a price... -')
    disp(' ')
end
%% apply all filters
% Fix Price Variable
if isVariableRet || isVariableRI   % if not prices are loaded 
    disp('load price variable')
    priceVariable = load([loadDirectory,'p.mat'],'p'); % use prices to track the zero or negatrive entries 
    priceVariable = priceVariable.p;
else
    disp('p variable is base...')
    priceVariable=variable;
end
% Fix RI Variable
if isVariableRI
    disp('ri variable is base')
    riVariable = variable;
else
    disp('loading ri variable')
    riVariable =  load([loadDirectory,'ri.mat'],'ri'); % Use Prices to get Negative-Zero entries.
    riVariable = riVariable.ri;
end
%% Filter 1 - Remove Invalid Negative-Zero Pricing Data -Applied it Initially.
isInvalidPricing = priceVariable<=0|riVariable<=0; % So to cover Both Possibilities of non sense values in pricing Data 
if filterIndex(1)
    disp('Applying Filter 1 - Non Sense price values -')
    isExclude(isInvalidPricing)=true;
    priceVariable(priceVariable<=0)=nan;
    priceVariable(riVariable<=0)=nan;
    if isVariableRet || isVariableRI
        variable(priceVariable<=0)=nan;
        variable(riVariable<=0)=nan;
    end
    % if Filtering is Applied to Returns, then we will need to correct t+1
    % return entry also, since a bad price will affect both t,t+1 returns.
    if isVariableRet
        isInvalidPricingForward = [false(1,size(isInvalidPricing,2));isInvalidPricing(1:end-1,:)];
        clear isInvalidPricing
        isExclude(isInvalidPricingForward)=true;
        clear isInvalidPricingForward
    else
        clear isInvalidPricing
    end
else
    disp('Skipping Filter 1 - Non Sense price values -')
end
%% Filter 2 - Re-Compute Incorectly Adjusted Data : 
if filterIndex(2)
    disp('Applying Filter 2 - Capital Adjustment Inconsistencies -')
    if isVariableRet || isVariableRI
        p=load([loadDirectory,'p.mat'],'p'); p=p.p; p(p<=0)=nan;
        up = load([loadDirectory,'up.mat'],'up'); up=up.up;up(up<=0)=nan;
        af = load([loadDirectory,'af.mat'],'af');af=af.af;af(af<=0)=nan;
        varDiff = ((abs((up.*af)-p))./(p)); clear p up af
        isExclude(varDiff>parameters.filters.adjustmentCutoff.price)=true;
        priceVariable(varDiff>parameters.filters.adjustmentCutoff.price)=nan;
        variable(varDiff>parameters.filters.adjustmentCutoff.price)=nan;
        clear varDiff
    else
        unadjVariableIndex = ['u',variableIndex];
        uVariable = load([loadDirectory,unadjVariableIndex,'.mat'],unadjVariableIndex);
        eval(['uVariable=uVariable.',unadjVariableIndex,';'])
        uVariable(uVariable<=0) = nan;
        af = load([loadDirectory,'af.mat'],'af'); af = af.af;
        af(af<=0) = nan;
        varDiff = ((abs((uVariable.*af)-priceVariable))./(priceVariable));
        isExclude(varDiff>parameters.filters.adjustmentCutoff.price)=true;
        af(varDiff>parameters.filters.adjustmentCutoff.price)=nan;
        priceVariableNew = uVariable.*af;
        clear uVariable unadjVariableIndex varDiff
        priceVariableNew(isnan(priceVariable)) = nan; % Adjust for the already applied filters.
        priceVariable = priceVariableNew; clear priceVariableNew
    end
else
    disp('Skipping Filter 2 - Capital Adjustment Inconsistencies -')
end
%% Compute Returns Variable After Initial Filters for capital Adjustments and non sense values  
if isVariableRet
    disp('Using Ret Variable to compute Returns')
    retVariable = variable; 
elseif isVariableRI
    disp('Using RI Variable to compute Returns')
    retVariable=[nan(1,size(variable,2));-1+variable(2:end,:)./variable(1:end-1,:)];
else
    disp('Using Price Variable to compute Returns')
    retVariable=[nan(1,size(priceVariable,2));-1+priceVariable(2:end,:)./priceVariable(1:end-1,:)];
end
%% Filter 3 : Remove all Prices of dead stocks
% Corrects the padded prices at the end of time series for dead-inactive
% stocks.
if filterIndex(3)
    disp('Applying Filter 3 - Dead Stocks -')
    for iStock=1:kStocks
        ixPotentialDead=find((retVariable(:,iStock)~=0) & ~isnan(retVariable(:,iStock)),1,'last'); % convoluted syntax handles nans properly
        if nDays-ixPotentialDead>parameters.filters.nDaydead % constant for more than the 2 last weeks
            isExclude(ixPotentialDead+1:end,iStock)=true;
        end
    end
    clear iStock ixPotentialDead
else
    disp('Skipping Filter 3 - Dead Stocks -')
end
%% Filter 4 :  Filter Long Term Illiquid Prices
% Essential When a zero return Price is enhanced with Dividends in a total return index. 
if filterIndex(4)
    disp('Applying Filter 4 - staleness -')
    isExclude(mvav(abs(retVariable),parameters.filters.nDayStale)==0)=true;
else
    disp('Skipping Filter 4 - staleness -')
end
%% Filter 5 : Massive Price-Return Daily Reversals
% Remove massive Price Reversals (of >100%, <-50%) within two days-
% if Rt >100% and (1+Rt)(1+Rt-1)-1< 0
if filterIndex(5)
    disp('Applying Filter 5 - Daily Reversals -')
    retLead=[retVariable(2:end,:);nan(1,size(retVariable,2))];
    isReturnReversal = (retVariable>parameters.filters.retcutoff & retLead<parameters.filters.retLeadcutoff)|...
        (retVariable<parameters.filters.retLeadcutoff & retLead>parameters.filters.retcutoff);
    isExclude(isReturnReversal)=true; % price daily reversals larger than 100%.
    % If Applied To rets Remove Next Day's Return Also, as it may be affected
    % by Reversal
    if isVariableRet
        isReturnReversalForward = [false(1,size(isReturnReversal,2));isReturnReversal(1:end-1,:)];
        clear isReturnReversal
        isExclude(isReturnReversalForward)=true;
        clear isReturnReversalForward
    else
        clear isReturnReversal
    end
    clear retLead
else
    disp('Skipping Filter 5 - Daily Reversals -')
end
%% Filter 6 : Exclude Problematic tiny stocks
if filterIndex(6)
    disp('Applying Filter 6 - Penny Stocks -')
    up = load([loadDirectory,'up.mat'],'up'); up=up.up;up(up<=0)=nan;
    for iDate=1:size(up,1)
        pennyBreakpointDay = prctile(up(iDate,:),parameters.filters.priceLevelcutoff);
        isExclude(iDate,up(iDate,:)<pennyBreakpointDay)=true;
    end
    clear up pennyBreakpointDay
else
    disp('Skipping Filter 6 - Penny Stocks -')
end
%
%% Filter 6 :  Lower Volatility Bounds 
if filterIndex(7)|| filterIndex(8)
    stdDev = nanstd(retVariable,[],1);
end
if filterIndex(7)
    disp('Applying Filter 7 - Low Volatility -')
    isExclude(:,stdDev<parameters.filters.lowerVolcutoff)=true;
else
    disp('Skipping Filter 7 - Low Volatility -')
end
%% Filter 8 : Upper Volatility Bounds % No such Filter to be applied a priori 
if filterIndex(8)
    disp('Applying Filter 8 - High Volatility -')
    isExclude(:,stdDev>parameters.filters.upperVolcutoff)=true;
else
    disp('Skipping Filter 8 - High Volatility -')
end
clear stdDev
%% Filters 8-10 Compute Stuff for all together and apply individually
% Compute Actual Days for all Stocks: 
endDate = nan(1,size(retVariable,2));
startDate = nan(1,size(retVariable,2));
for iStock=1:size(retVariable,2)
    lastIndex = find(~isnan(retVariable(:,iStock)),1,'last');
    if ~isempty(lastIndex)
         endDate(1,iStock)=lastIndex;
    end
    startIndex = find(~isnan(retVariable(:,iStock)),1,'first');
    if ~isempty(startIndex)
        startDate(1,iStock)=startIndex;
    end
end
clear iStock lastIndex startIndex 
% Define Some Useful Ratios to Use on Filters : 
positiveReturnsCounter = nan(1,size(retVariable,2));
negativeReturnsCounter = nan(1,size(retVariable,2));
illiquidityCounter = nan(1,size(retVariable,2));
totalActualDays = nan(1,size(retVariable,2));
totalTradingDays = nan(1,size(retVariable,2));
%nanDaysCounter = nan(1,size(retVariable,2));
for iStock = 1:size(retVariable,2)
    totalActualDays(1,iStock) = endDate(1,iStock)-startDate(1,iStock)+1;
    if ~isnan(endDate(1,iStock)) && ~isnan(startDate(1,iStock))
      reti = retVariable(startDate(1,iStock):endDate(1,iStock),iStock);
      totalTradingDays(1,iStock) = sum(~isnan(reti));
      %nanDaysCounter(1,iStock) = sum(isnan(reti));
      illiquidityCounter(1,iStock)= sum(reti==0);
      positiveReturnsCounter(1,iStock)=sum(reti>0);
      negativeReturnsCounter(1,iStock)=sum(reti<0);
    else 
        illiquidityCounter(1,iStock)=nan;
        positiveReturnsCounter(1,iStock)=nan;
        negativeReturnsCounter(1,iStock)=nan;
        totalTradingDays(1,iStock) = nan;
        %nanDaysCounter(1,iStock)=nan;
    end
    clear reti 
end
clear iStock endDate
%%%%%% Filter 9 %%%%%%%
if filterIndex(9)
    disp('Applying Filter 9 - Illiquidity -')
    % Filter 9a
    isLongTermIlliquid = (totalTradingDays./totalActualDays)<parameters.filters.longTermIlliquidCutoff;
    isExclude(:,isLongTermIlliquid)=true;
    clear isLongTermIlliquid
    % Filter 9b
    isShortTermIlliquid = (illiquidityCounter./totalTradingDays)>parameters.filters.shortTermIlliquidCutoff;
    isExclude(:,isShortTermIlliquid)=true;
    clear isShortTermIlliquid
else
    disp('Skipping Filter 9 - Illiquidity -')
end
clear totalActualDays
%%%%%%%% Filter 10 %%%%%%%
if filterIndex(10)
    disp('Applying Filter 10 - Implausibility -')
    isExtremePositiveOrNegative = (((positiveReturnsCounter+illiquidityCounter)./totalTradingDays)>parameters.filters.extremePerformanceCutoff | ((negativeReturnsCounter+illiquidityCounter)./totalTradingDays)>parameters.filters.extremePerformanceCutoff);
    isExclude(:,isExtremePositiveOrNegative)=true;
    clear isExtremePositiveOrNegative
else
    disp('Skipping Filter 10 - Implausibility -')
end
clear positiveReturnsCounter negativeReturnsCounter illiquidityCounter
%%%%%%%% Filter 11 %%%%%%%%
% Filter 11
if filterIndex(11)
    disp('Applying Filter 11 - Few Observations -')
    % Filter 11a
    isTooShortValid = (totalTradingDays<parameters.filters.shortPeriodCutoff & startDate<size(retVariable,1)-parameters.filters.shortPeriodCutoff);
    isExclude(:,isTooShortValid)=true;
    clear isTooShortValid
    % Filter 11b
    isExtraShort  = totalTradingDays<parameters.filters.extraShortCutoff;
    isExclude(:,isExtraShort)=true;
    clear isExtraShort
else
    disp('Skipping Filter 11 - Few Observations -')
end
clear totalTradingDays startDate
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Static Filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Common Stocks
if filterIndex(12)
    disp('Applying Filter 12 - Common Stock Filters -')
    shareType = load([loadDirectory,'shareType.mat'],'isCommon_trac','isNonCommon_abbs');
    isExclude(~shareType.isCommon_trac|shareType.isNonCommon_abbs)=true;
    clear shareType
else
    disp('Skipping Filter 12 - Common Stock Filters -')
end
% Stocks Not trading in Home Currency
if filterIndex(13)
    disp('Applying Filter 13 - Stocks Not Trading in Home Currency -')
    currencyCode = load([loadDirectory,'currencyCode.mat'],'isHomeCurrency');
    isExclude(~currencyCode.isHomeCurrency)=true;
    clear currencyCode
else
    disp('Skipping Filter 13 - Stocks Not Trading in Home Currency -')
end
% Non Domestic Stocks
if filterIndex(14)
    disp('Applying Filter 14 - Non Domestic Stocks -')
    geogn = load([loadDirectory,'geogn.mat'],'isDomestic');
    isExclude(~geogn.isDomestic)=true;
    clear geogn
else
    disp('Skipping Filter 14 - Non Domestic Stocks -')
end
% Cross Listings 
if filterIndex(15)
    disp('Applying Filter 15 - Cross Listings -')
    exchangeCode = load([loadDirectory,'exchangeCode.mat'],'isRemoved_crossList');
    isExclude(exchangeCode.isRemoved_crossList)=true;
    clear exchangeCode
else
    disp('Skipping Filter 15 - Cross Listings -')
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Apply Filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample Filter - Eliminate Non Trading Days, After Filtering
if isVariableRet || isVariableRI
    variableClean = variable; clear priceVariable variable
else
    variableClean = priceVariable; clear priceVariable variable
end
variableClean(isExclude) = nan; 

if isVariableRet
    retVariableClean = variableClean; 
else
    retVariableClean=[nan(1,size(variableClean,2));-1+variableClean(2:end,:)./variableClean(1:end-1,:)];
end
%% remove non actual trading dates
% All markets Except German Venues which need to be combined later. 
if filterIndex(16)
    disp('Applying Filter 16 - Remove Holidays After Filters -')
    variableClean(mean(~isnan(retVariableClean) & retVariableClean~=0,2)<parameters.filters.validStockDay,:)=nan;
else
    disp('Skipping Filter 16 - Market Holidays -')
end
clear retvariableClean
%% Apply Filters to market cap
mv = load([loadDirectory,'mv.mat'],'mv'); 
mvClean = mv.mv; clear mv
mvClean(isnan(variableClean))=nan;
%{
% Do This outside this function 
%% Remove Stocks and dates but Keep Indexes for rest of the data. 
isAllNanStocks=mean(isnan(variableClean))==1;
isAllNanDays=mean(isnan(variableClean),2)==1;   
% Notice the Complete Filtering of Stocks or Days : 
if sum(isAllNanStocks)==size(variableClean,2) || sum(isAllNanDays)==size(variableClean,1)
    if sum(isAllNanStocks)==size(variableClean,2)
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Stocks Completely Removed')
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    end
    if sum(isAllNanDays)==size(variableClean,1)
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Days Completely Removed')
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    end
end
%% Remove Stocks for all variables and keep the variables to be used. 
variableClean=variableClean(~isAllNanDays,~isAllNanStocks);
stocksUse=stocksUse(~isAllNanStocks);
tUse=tUse(~isAllNanDays);
mvClean=mvClean(~isAllNanDays,~isAllNanStocks);
ixStocksUse = ~isAllNanStocks; 
ixDateUse = ~isAllNanDays; 
%}
end
%%%%%%%%%%%%%%%%%%%%%%%
%% Embeeded Functions %
%%%%%%%%%%%%%%%%%%%%%%%
function xMvav=mvav(X,windowLength,Mask)
% mvav(X,windowLength,Mask)
% Mask  - (Optional input) binary mask for X showing observations that need to be used when
% taking lags.

if nargin==2
    xMvav=filter(ones(1,windowLength)/windowLength,1,X,NaN(windowLength-1,1)); 
else
    xMvav=NaN(size(X));
    xMvav(Mask)=filter(ones(1,windowLength)/windowLength,1,X(Mask),NaN(windowLength-1,1)); 
end
end