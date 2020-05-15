close all
clear all
ImportELdata

%% 
figure (1)
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
saveas(gcf,'Signature2017.png')


figure (2)
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
saveas(gcf,'Signature2018.png')


figure (3)
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
saveas(gcf,'Signature2019.png')
%% Outside temperatures

x = linspace(1,height(EL_Heat),36)
x2 = datetime(EL_Heat{:,1},'InputFormat','dd.MM.yyyy')

figure (10)
plot(x2,EL_Heat{:,2},'-*')
xlabel('Month')
ylabel('Temperature [°C]')
title('Outside temperatures')
saveas(gcf,'OutsideTemperatures.png')

figure (20)
plot(x2,EL_Heat{:,3},'*')

%% Fit of the signatures

x = EL_Heat(:,2);
y = EL_Heat(:,3);

for building = 3:size(EL_Heat,2)
    FitBuilding1 = [];
    FitBuilding2 = [];
%     building = 3
    figure (building+30)
    hold on
    title (['Building ', EL_Heat.Properties.VariableNames{building}(1:3)])
    for year = 1:3
%         building = 3;
%         year = 1;
        criteria = EL_Heat{(year-1)*12+1:(year-1)*12+12,building}>= 100;
        y = EL_Heat{(year-1)*12+1:(year-1)*12+12,building};
        x = EL_Heat{(year-1)*12+1:(year-1)*12+12,2};
        
        x2 = x(criteria);
        y2 = y(criteria);
%         [t1,t2] = createFit(x2,y2)
        % Fit model to data.
        [fitresult, gof] = fit( x2, y2, fittype( 'poly1'));
       
        h = plot( fitresult, x2, y2 );

    end
    
end

