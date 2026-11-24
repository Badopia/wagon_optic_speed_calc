%% demo_acceleration_comp.m
% Демонстрация и поверка методики Strauß 2009 на данных optic_vagon.
%
% Скрипт запускает функцию acceleration_compensation в обоих режимах
% выделения окон ('feature' и 'fixed'), сравнивает результат с исходным
% методом main.m (без warping) и строит графики.
%
% Запуск: из папки optic_vagon\acceleration_comp, командой demo_acceleration_comp
%   (папка определяется по mfilename, поэтому запускать можно и из Matlab,
%   если текущая папка другая).

clear; close all; clc;

thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);                       % папка optic_vagon
addpath(thisDir);                                   % путь к функции

% --- Профили для проверки ------------------------------------------------
% Профиль выбирается ОДИН раз. Раньше переменная profileFile
% переопределялась подряд три раза, и в работе оставался только последний
% (step_prof.mat) - именно этого хотели избежать, поэтому теперь список
% явный и цикл идёт по всем имеющимся профилям.
profiles = struct( ...
    'name',  {'wagon_12_6995', 'sine'}, ...
    'file',  {fullfile(rootDir, 'data', 'wagon_12_6995_profile.mat'), ...
              fullfile(rootDir, 'data', 'sine_prof.mat')}, ...
    'note',  {'синтетический профиль полувагона 12-6995', ...
              'синус по x: периодический рельеф, для проверки ложных пиков'});
if isfile(fullfile(rootDir, 'data', 'step_prof.mat'))
    profiles(end+1) = struct('name', 'step', ...
        'file', fullfile(rootDir, 'data', 'step_prof.mat'), ...
        'note', 'прямоугольные ступеньки: вырожденный рельеф');
end

fprintf('Скорость 4.0 м/с, ускорение -0.5 м/с^2, база 0.5 м, шум 0.033 м\n\n');

%% --- 0. Сквозная поверка математики (проверочный скрипт) ---------------
fprintf('--- 0. Математическая поверка (verify_acceleration_comp) ---\n');
verify_acceleration_comp;

%% --- 1. Режим 'fixed': регулярные окна с перекрытием ---------------------
% Окно 3.5 с выбрано осознанно: ускорение оценивается как a = 2*q*v/D,
% поэтому чем короче окно, тем сильнее ошибка q раздувает разброс a.
optsFixed = struct('windowMode', 'fixed', 'winSeconds', 3.5, ...
                   'hopSeconds', 0.5, 'randomSeed', 42, 'makePlot', false);

%% --- 2. Режим 'feature': окна вокруг деталей вагона ----------------------
% Детали дают мало времени, поэтому окно расширяется до минимума,
% достаточного для корреляции; такие окна помечаются 'ledge'.
optsFeat = struct('windowMode', 'feature', 'randomSeed', 42, 'makePlot', false);

results = struct('name', {}, 'rFixed', {}, 'rFeat', {});
for p = 1:numel(profiles)
    if ~isfile(profiles(p).file)
        fprintf('Профиль не найден, пропуск: %s\n', profiles(p).file);
        continue
    end
    fprintf('\n================= %s =================\n', profiles(p).name);
    fprintf('%s\n\n', profiles(p).note);

    rFixed = acceleration_compensation(profiles(p).file, optsFixed);
    rFeat  = acceleration_compensation(profiles(p).file, optsFeat);

    fprintf('--- %s: fixed (окно %.1f с, шаг %.1f с) ---\n', ...
        profiles(p).name, optsFixed.winSeconds, optsFixed.hopSeconds);
    printTable(rFixed.table);
    printQuality('fixed', rFixed);

    fprintf('\n--- %s: feature (окна вокруг деталей) ---\n', profiles(p).name);
    printTable(rFeat.table);
    printQuality('feature', rFeat);

    results(end+1) = struct('name', profiles(p).name, ...
        'rFixed', rFixed, 'rFeat', rFeat); %#ok<SAGROW>
end

%% --- 3. Сводка ----------------------------------------------------------
fprintf('\n=== ИТОГ ===\n');
fprintf('%-16s %-10s %8s %8s %10s %10s\n', ...
    'профиль', 'режим', 'окон', 'принято', 'RMSE v', 'mean a');
for k = 1:numel(results)
    q1 = results(k).rFixed;
    fprintf('%-16s %-10s %8d %8d %10.4f %10.4f\n', results(k).name, 'fixed', ...
        q1.quality.count, q1.quality.accepted, q1.quality.rmseSpeed, ...
        q1.quality.meanAcc);
    q2 = results(k).rFeat;
    fprintf('%-16s %-10s %8d %8d %10.4f %10.4f\n', results(k).name, 'feature', ...
        q2.quality.count, q2.quality.accepted, q2.quality.rmseSpeed, ...
        q2.quality.meanAcc);
end

%% --- 4. Графики ---------------------------------------------------------
% Берём профиль вагона - на нём видно физику задачи.
idxWagon = find(strcmp({results.name}, 'wagon_12_6995'), 1);
if ~isempty(idxWagon)
    r2 = results(idxWagon).rFixed;
    r1 = results(idxWagon).rFeat;

    fig = figure('Name', 'Скорость и ускорение по методике Strauß 2009', ...
        'NumberTitle', 'off', 'WindowState', 'normal');
    t = r2.signals.t;

    subplot(3, 1, 1);
    plot(t, r2.signals.y1, 'k', 'LineWidth', 1.0); hold on;
    plot(t, r2.signals.y2, 'r', 'LineWidth', 0.8);
    for k = 1:height(r2.table)
        if r2.table.Accepted(k)
            xline(r2.table.Start_s(k), '-b');
        elseif r2.table.Region(k) == "ledge"
            xline(r2.table.Start_s(k), ':m');
        else
            xline(r2.table.Start_s(k), ':k');
        end
    end
    grid on; xlabel('Время, с'); ylabel('Расстояние, м');
    title('Сигналы дальномеров; сплошная - принято, пунктир - отказ, magenta - кромка');
    legend('Датчик 1', 'Датчик 2', 'Location', 'best');

    subplot(3, 1, 2);
    plot(r2.table.tMid_s, r2.table.SpeedTrue_mps, 'k-', 'LineWidth', 1.5); hold on;
    ok = r2.table.Accepted;
    plot(r2.table.tMid_s(ok), r2.table.Speed_mps(ok), 'bo', ...
        'MarkerFaceColor', 'b');
    plot(r2.table.tMid_s(~ok), r2.table.Speed_mps(~ok), 'ko', 'MarkerSize', 9);
    grid on; xlabel('Время, с'); ylabel('Скорость, м/с');
    legend('Истинная', 'Принята', 'Отброшена', 'Location', 'northeast');
    title('Оценка скорости: v = L/(lag*Ts)');

    subplot(3, 1, 3);
    yline(r2.settings.accel0, 'k-', 'LineWidth', 1.5); hold on;
    plot(r2.table.tMid_s(ok), r2.table.Acc_mps2(ok), 'rs', ...
        'MarkerFaceColor', 'r');
    grid on; xlabel('Время, с'); ylabel('Ускорение, м/с^2');
    legend('Истинное', 'Оценка', 'Location', 'best');
    title('Оценка ускорения из параметра q: a = 2*q*v/D');
    sgtitle(sprintf('Профиль %s: база %.2f м, окно %.1f с, шаг %.1f с', ...
        results(idxWagon).name, r2.settings.sensorGap, ...
        r2.settings.winSeconds, r2.settings.hopSeconds));

    outDir = fullfile(thisDir, 'results');
    if ~isfolder(outDir)
        mkdir(outDir);
    end
    png = fullfile(outDir, 'demo_modes.png');
    exportgraphics(fig, png, 'Resolution', 130);
    fprintf('\nГрафик: %s\n', png);
end

%% ==================== вспомогательные функции ============================
function printTable(tbl)
% Печать таблицы результатов в компактном виде.
% Столбец "итог": OK - надежное окно, OK* - принято, но окно на
% передней кромке вагона (в нём каналы видят разные части вагона).
    fprintf('%-7s %-7s %-9s %-18s %9s %9s %8s %8s %7s %-6s %s\n', ...
        'Start', 'End', 'Region', 'Status', 'v, м/с', 'v ист.', ...
        'a, м/с^2', 'rho', 'q', 'итог', 'причина');
    for k = 1:height(tbl)
        verdict = "нет";
        if tbl.Accepted(k) && tbl.Region(k) == "run"
            verdict = "OK";
        elseif tbl.Accepted(k)
            verdict = "OK*";
        end
        fprintf('%7.2f %7.2f %-9s %-18s %9.4f %9.4f %8.3f %8.4f %7.3f %-6s %s\n', ...
            tbl.Start_s(k), tbl.End_s(k), tbl.Region(k), ...
            string(tbl.Status(k)), tbl.Speed_mps(k), ...
            tbl.SpeedTrue_mps(k), tbl.Acc_mps2(k), ...
            tbl.Correlation(k), tbl.WarpFactor(k), verdict, ...
            strjoin(cellstr(tbl.Note(k)), ''));
    end
end

function printQuality(name, r)
% Печать сводных метрик качества с раздельной статистикой по регионам.
    q = r.quality;
    if q.accepted == 0
        fprintf('  %s: ни одного принятого окна из %d\n', name, q.count);
        return
    end
    fprintf('  %s: принято %d из %d окон\n', name, q.accepted, q.count);
    fprintf('    RMSE скорости = %.4f м/с, max |ошибка| = %.4f м/с\n', ...
        q.rmseSpeed, q.maxSpeedErr);
    fprintf('    среднее ускорение = %+.4f м/с^2 (истина %+.3f)\n', ...
        q.meanAcc, r.settings.accel0);
    fprintf('    регион run: %d принято, RMSE %.4f | ledge: %d принято, RMSE %.4f\n', ...
        q.runCount, q.runRmseSpeed, q.ledgeCount, q.ledgeRmseSpeed);
end



