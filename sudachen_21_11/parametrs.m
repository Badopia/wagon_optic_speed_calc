clear
%% Датчики и их расположение
delta = 50;         % База двух датчиков вдоль пути, см; 
T = 7;              % Длительность моделирования, с
d0 = 0.40;          % Боковой зазор от датчика до плоскости борта, м
L = delta / 100;    % База датчиков.м 

%% Движение вагона
speed = 2;          % Начальная скорость, м/с
acceleration = 0.5; % Постоянное ускорение, м/с^2


%% Параметры дальномеров DE-CM2THIU
Ts = 0.010;          % Шаг записи и моделирования, с (100 Гц)
FoV = 2;             % Дальность, м; используется при расчёте maxRange
std = 0.003;         % СКО шума расстояния, м
rangeMargin = 0.995;         % Множитель запаса дальности
maxRange = FoV * rangeMargin; % Верхний предел измерения, м

%% Файл одного профиля вагона
dataFolder = fullfile(fileparts(mfilename('fullpath')), 'data');
profileFile = fullfile(dataFolder, 'wagon_12_6995_profile.mat'); % Загружается в main.m

%% Фон датчиков
bg = 2;                    % Значение при отсутствии возврата, м

%% Корреляционная оценка скорости и ускорения — solve.m
vMin = 2;                     % Нижняя граница поиска скорости, м/с
vMax = 10;                    % Верхняя граница поиска скорости, м/с
minCorr = 0.80;              % Порог корреляции из условия peak(k) >= 0.80

% Потоковый расчёт по скользящему окну с компенсацией ускорения.
solveWindow = 2;             % Длина окна, с; первая оценка после накопления окна
solveHop = L/vMin;           % Период обновления оценки, с
solveQLimit = 0.8;           % Ограничение параметра деформации |q| < 1
solveQStep = 0.01;           % Шаг поиска q

%% Обработка двух потоков со временем — solve.m
analysis=struct('Ts',Ts,'L',L,'vMin',vMin,'vMax',vMax, ...
    'windowSeconds',solveWindow,'overlapFraction',1-solveHop/solveWindow, ...
    'bg',bg,'minValidFraction',0.50,'minValidSamples',50, ...
    'minRange',0.03,'minStd',0.01,'minChanges',2,'changeThreshold',0.02, ...
    'smoothSamples',3,'maxGapSeconds',2.5*Ts, ...
    'minCorrelation',0.90,'minLagPeakMargin',0.001,'minQContrast',0.0001, ...
    'rejectBoundary',true, ...
    'printWindows',true, ...
    'motion',struct('minCorr',minCorr,'qLimit',solveQLimit,'qStep',solveQStep));
% Второй участок всегда продлевается на L/vMin секунд.
% Ускорение удерживается до следующего принятого окна, скорость меняется по v=vm+a*(t-tMid).
% После изменения параметров запустите main.m и начните новый поток с state=[].
