%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /home/thomas/ownCloud2/MES3/Semester project/EL - Consommations Electricité Chaleur 2017-2018-2019.xlsx
%    Worksheet: Feuil1
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Auto-generated by MATLAB on 2020/04/30 14:59:05

%% Import the data
[~, ~, raw] = xlsread('/home/thomas/ownCloud2/MES3/Semester project/EL - Consommations Electricité Chaleur 2017-2018-2019.xlsx','Feuil1');
raw = raw(4:end,:);
stringVectors = string(raw(:,[2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95,98,101,104,106,107]));
stringVectors(ismissing(stringVectors)) = '';
raw = raw(:,[1,3,4,6,7,9,10,12,13,15,16,18,19,21,22,24,25,27,28,30,31,33,34,36,37,39,40,42,43,45,46,48,49,51,52,54,55,57,58,60,61,63,64,66,67,69,70,72,73,75,76,78,79,81,82,84,85,87,88,90,91,93,94,96,97,99,100,102,103,105,106,108]);

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Create table
ELbuildings = table;

%% Allocate imported array to column variable names
ELbuildings.EndOfMonth = stringVectors(:,37);
ELbuildings.MeanTemp = data(:,72);
ELbuildings.ELL_ElecForce = data(:,2);
ELbuildings.ELL_ElecService = data(:,4);
ELbuildings.ELL_ElecLight = data(:,6);
ELbuildings.ELL_ElecLightOutside = data(:,8);
ELbuildings.ELL_Heat = data(:,10);
ELbuildings.ELAELB_ElecForce = data(:,12);
ELbuildings.ELAELB_ElecService = data(:,14);
ELbuildings.ELAELB_ElecLight = data(:,16);
ELbuildings.ELAELB_ElecCafeteriaForce = data(:,18);
ELbuildings.ELAELB_ElecCafeteriaECS= data(:,20);
ELbuildings.ELA_Heat = data(:,22);
ELbuildings.ELB_Heat = data(:,24);
ELbuildings.ELDELEDIA_ElecForce = data(:,26);
ELbuildings.ELDELEDIA_ElecService = data(:,28);
ELbuildings.ELDELEDIA_ElecLight = data(:,30);
ELbuildings.ELDELEDIA_ElecLightOutside = data(:,32);
ELbuildings.ELD_Heat = data(:,34);
ELbuildings.ELE_Heat = data(:,36);
ELbuildings.DIA_Heat = data(:,38);
ELbuildings.ELGELH_ElecForce = data(:,40);
ELbuildings.ELGELH_ElecService = data(:,42);
ELbuildings.ELGELH_ElecLight = data(:,44);
ELbuildings.ELGELH_ElecSM2Lab = data(:,46);
ELbuildings.ELG_Heat = data(:,48);
ELbuildings.ELH_Heat = data(:,50);
ELbuildings.UOR_ElecLight = data(:,52);
ELbuildings.UOR_ElecForce = data(:,54);
ELbuildings.UOR_ElecLight = data(:,56);
ELbuildings.UORMJC_ElecForce = data(:,58);
ELbuildings.UORMJC_ElecLight = data(:,60);
ELbuildings.UORMJC_ElecForce = data(:,62);
ELbuildings.UOR_ElecPrincipalLight = data(:,64);
ELbuildings.UOR_ElecWelcomeForce = data(:,66);
ELbuildings.UOR_ElecWelcomeLight = data(:,68);
ELbuildings.UOR_Heat = data(:,70);


%% Clear temporary variables
clearvars data raw stringVectors;

ELL_Elec = ELbuildings(:,[1:2,3:6])
ELA_ELB_Elec = ELbuildings(:,[1:2,8:12])
ELD_ELE_DIA_Elec= ELbuildings(:,[1:2,15:18])
ELG_ELH_Elec= ELbuildings(:,[1:2,22:25])
UOR_Elec = ELbuildings(:,[1:2,28:34])
EL_Heat = ELbuildings(:,[1,2,7, 13,14,19,20,21,26,27,35])