% load data and convert monthly observations into quarterly
cd('Data')
picPath = [pwd, '/data_pics/'];


%dates for datasets that are not in differences
t1 = datetime('1984-Jan-01');
t2 = datetime('2007-Dec-01');
t = t1:calmonths(3):t2;

%dates for datasets that are in differences
t1 = datetime('1947-Jan-01');
t2 = datetime('2014-Oct-01');
t_full = t1:calmonths(3):t2;

t1_plot = datetime('1965-Jan-01');
t2_plot = t2;
t_plot = t1_plot:calmonths(3):t2_plot;
dates_skip = length(t_full) - length(t_plot);



startDate = datenum(0084,01,01);
startDateDiff = datenum(0083,08,01);
endDate = datenum(0108,01,01);

GDP = xlsread('GDP.xls');           %nominal GDP - used for labor share

GDPdat = xlsread('GDPC96.xls');
    %real GDP, seasonally adjusted, annual rate, quarterly from 1947Q1

%note: that the files below needed to be cleaned, where the notes section 
%of each was deleted, so matlab could read the file

FedFunds = xlsread('FEDFUNDS.xls');
    %Effective Fed Funds rate, not adjusted, monthly from July 1954
GDPDEF = xlsread('GDPDEF');
    %Implicit Price Deflator, seasonally adjusted, from 1947Q1
LaborIncome = xlsread('COE.xls');
    %compensation of employees, seasonally adjusted
Hours = xlsread('HOHWMN02USQ065S.xls');
    %weekly hours, quarterly, from Jan 1955

%plot the labor share
p = plot(t_plot,log(LaborIncome(dates_skip+1:end,2)./GDP(dates_skip+1:end,2)), 'linewidth',4,...
    'DatetimeTickFormat','yyyy');
set(gca,'fontsize',20,'fontweight','demi')
ylim([-0.7 -0.5])
xlim([datenum(1965,1,1) datenum(2016,1,1)])
set(gca, 'Xtick', datenum(1965,1,1):(datenum(1975,2,1)-datenum(1965,1,1)): datenum(2016,4,4) )
box off
print('-dpng',[picPath, 'lab_share_full'])

PCE_S = xlsread('PCESV.xls');
%personal consumption expeditures, services, seasonally adjusted annual
%rate
PCE_ND = xlsread('PCND.xls');
%personal consumption expenditures, nondurable goods, seasonally adjusted,
%quarterly at annual rate


%plot consumption share based on nominal GDP
plot(t_plot,log((PCE_S(dates_skip+1:end,2) + PCE_ND(dates_skip+1:end,2))./GDP(dates_skip+1:end,2)),...
    'linewidth',4,...
    'DatetimeTickFormat','yyyy');
set(gca,'fontsize',20,'fontweight','demi')
ylim([-1 -0.2])
xlim([datenum(1965,1,1) datenum(2016,1,1)])
set(gca, 'Xtick', datenum(1965,1,1):(datenum(1975,2,1)-datenum(1965,1,1)): datenum(2016,4,4) )
box off
print('-dpng',[picPath, 'cons_share_S_ND'])

GDP = GDP(datenum(datestr(GDP(:,1)))>=startDate & datenum(datestr(GDP(:,1)))<endDate,2) ;

    
%Delete all antries from before 1955 Q1, since they are not common    
Hours = log( Hours(datenum(datestr(Hours(:,1)))>=startDate & datenum(datestr(Hours(:,1)))<endDate,2) );
LaborIncome = LaborIncome(datenum(datestr(LaborIncome(:,1)))>=startDate & datenum(datestr(LaborIncome(:,1)))<endDate,2) ;

LabShare = log(LaborIncome./GDP);

%this is what we use for now, but the 50s are too low, and the entire
%series is too low

GDPDEF = GDPDEF(datenum(datestr(GDPDEF(:,1)))>=startDateDiff & datenum(datestr(GDPDEF(:,1)))<endDate,2);
GDPDEF = diff(log(GDPDEF));

FedFunds =  FedFunds(datenum(datestr(FedFunds(:,1)))>=startDate & datenum(datestr(FedFunds(:,1)))<endDate,2) ;
GDPdat = GDPdat(datenum(datestr(GDPdat(:,1)))>=startDateDiff & datenum(datestr(GDPdat(:,1)))<endDate,2);

GDPdat = diff(log(GDPdat));

%Common start date is 1955 Q1

%convert federal funds from monthly to quarterly by taking the mean
nQuarters = length(GDPdat);
FedFundsQuarterly = zeros(nQuarters,1);
j = 1;
for i = 1:nQuarters
    FedFundsQuarterly(i) = mean( FedFunds(j:j+2));
    j = j+3;
end

FedFundsQuarterly = FedFundsQuarterly/400;

%save all datasets
save('Hours.mat', 'Hours')
save('GDP.mat', 'GDPdat')
save('FedFunds.mat', 'FedFundsQuarterly')
save('Labor.mat', 'LabShare')
save('GDPdef.mat', 'GDPDEF')


