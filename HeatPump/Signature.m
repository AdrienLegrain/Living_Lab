close all
clear all
ImportELdata

%% 
figure('Name','Data','Units','inches','Position',[0 0 7.5 5.5],'PaperPositionMode','auto')
hold on

subplot(2,2,1)
hold on
plot(EL_Heat{1:12,2},EL_Heat{1:12,3},'-*')
plot(EL_Heat{1:12,2},EL_Heat{1:12,4},'-*')
plot(EL_Heat{1:12,2},EL_Heat{1:12,5},'-*')
plot(EL_Heat{1:12,2},EL_Heat{1:12,6},'-*')
plot(EL_Heat{1:12,2},EL_Heat{1:12,7},'-*')
plot(EL_Heat{1:12,2},EL_Heat{1:12,8},'-*')
plot(EL_Heat{1:12,2},EL_Heat{1:12,9},'-*')
plot(EL_Heat{1:12,2},EL_Heat{1:12,10},'-*')
plot(EL_Heat{1:12,2},EL_Heat{1:12,11},'-*')
title('Data from 2017')
xlabel('Outside temperature [°C/Month]')
ylabel('Energy for heating [kWh]')
legend('ELL','ELA','ELB','ELD','ELE','DIA','ELG','ELH','UOR')
% saveas(gcf,'Signature2017.png')

subplot(2,2,2)
hold on
plot(EL_Heat{13:24,2},EL_Heat{13:24,3},'-*')
plot(EL_Heat{13:24,2},EL_Heat{13:24,4},'-*')
plot(EL_Heat{13:24,2},EL_Heat{13:24,5},'-*')
plot(EL_Heat{13:24,2},EL_Heat{13:24,6},'-*')
plot(EL_Heat{13:24,2},EL_Heat{13:24,7},'-*')
plot(EL_Heat{13:24,2},EL_Heat{13:24,8},'-*')
plot(EL_Heat{13:24,2},EL_Heat{13:24,9},'-*')
plot(EL_Heat{13:24,2},EL_Heat{13:24,10},'-*')
plot(EL_Heat{13:24,2},EL_Heat{13:24,11},'-*')
title('Data from 2018')
xlabel('Outside temperature [°C/Month]')
ylabel('Energy for heating [kWh]')
legend('ELL','ELA','ELB','ELD','ELE','DIA','ELG','ELH','UOR')
% saveas(gcf,'Signature2018.png')


subplot(2,2,3)
hold on
plot(EL_Heat{25:36,2},EL_Heat{25:36,3},'-*')
plot(EL_Heat{25:36,2},EL_Heat{25:36,4},'-*')
plot(EL_Heat{25:36,2},EL_Heat{25:36,5},'-*')
plot(EL_Heat{25:36,2},EL_Heat{25:36,6},'-*')
plot(EL_Heat{25:36,2},EL_Heat{25:36,7},'-*')
plot(EL_Heat{25:36,2},EL_Heat{25:36,8},'-*')
plot(EL_Heat{25:36,2},EL_Heat{25:36,9},'-*')
plot(EL_Heat{25:36,2},EL_Heat{25:36,10},'-*')
plot(EL_Heat{25:36,2},EL_Heat{25:36,11},'-*')
title('Data from 2019')
xlabel('Outside temperature [°C/Month]')
ylabel('Energy for heating [kWh]')
legend('ELL','ELA','ELB','ELD','ELE','DIA','ELG','ELH','UOR')
% saveas(gcf,'Signature2019.png')

% Outside temperatures
x = linspace(1,height(EL_Heat),36);
x2 = datetime(EL_Heat{:,1},'InputFormat','dd.MM.yyyy');

subplot(2,2,4)
plot(x2,EL_Heat{:,2},'-*')
xlabel('Month')
ylabel('Temperature [°C]')
title('Outside temperatures')
% saveas(gcf,'OutsideTemperatures.png')

print -depsc2 dataHeating.eps


%% Fit of the signatures

x = EL_Heat(:,2);
y = EL_Heat(:,3);
FitToSave = [1,2,2,1,2,1,2,1,1];
FitsValue = [];
figure('Name','Fits','Units','inches','Position',[0 0 7.5 5.5],'PaperPositionMode','auto')
hold on
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
for building = 3:size(EL_Heat,2)
    subplot(3,3,building-2)
    title (EL_Heat.Properties.VariableNames{building}(1:3))
    hold on
    for year = 1:3
        criteria = EL_Heat{(year-1)*12+1:(year-1)*12+12,building}>= 100;
        y = EL_Heat{(year-1)*12+1:(year-1)*12+12,building};
        x = EL_Heat{(year-1)*12+1:(year-1)*12+12,2};
        x2 = x(criteria);
        y2 = y(criteria);
        [fitresult, gof] = fit( x2, y2, fittype( 'poly1'));
        if FitToSave(building-2)==year
            FitsValue = [FitsValue; fitresult.p1, fitresult.p2];
        end
        
        if year==1
            h = plot( fitresult, x2, y2);
            h(2).Color=[1 0 0.2];
        elseif year == 2
            h = plot( fitresult, x2, y2);
            h(2).Color=[0.2 0.8 0.2];
        else
            h = plot( fitresult, x2, y2);
            h(2).Color=[1 0.8 0];
        end
    end
    xlabel('T [°C]')
    ylabel('Energy [kWh]')
    legend('off')
end
print -depsc2 Fits.eps

