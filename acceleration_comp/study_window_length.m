%% study_window_length.m
% Исследование зависимости точности оценки скорости и ускорения от длины
% окна и его перекрытия (шага скольжения).
%
% ПЛАН ЭКСПЕРИМЕНТА
%   1. Основной перебор: длина окна x перекрытие x прогоны шума.
%   2. Повтор по начальным условиям (v0, a0): даёт статистику, не зависящую
%      от одной "удачной" конфигурации профиля.
%   3. Разделение точности и разрешения: показывается, что перекрытие
%      управляет ЧАСТОТОЙ оценок, а не их точностью.
%
% Метрики: RMSE (полная ошибка = bias + разброс), bias (систематика),
%          std (разброс), покрытие записи.
%
% Запуск: из папки optic_vagon\acceleration_comp, команда study_window_length

clear; close all; clc;

thisDir  = fileparts(mfilename('fullpath'));
rootDir  = fileparts(thisDir);
profileFile = fullfile(rootDir, 'data', 'wagon_12_6995_profile.mat');
addpath(thisDir);

outDir = fullfile(thisDir, 'results');
if ~isfolder(outDir)
    mkdir(outDir);
end

fprintf('Профиль: %s\n', profileFile);
fprintf('Результаты: %s\n\n', outDir);

%% --- Параметры плана -----------------------------------------------------
winList     = [2.0 2.5 3.0 3.5 4.0];                     % длины окон, с
hopFracList = [1.00 0.50 0.25];                         % шаг как доля окна
v0List      = [2.0 2.5];                                % начальные скорости, м/с
a0List      = 0.5;                                     % ускорения, м/с^2
nSeed       = 3;                                        % прогонов шума

% Обоснование плана.
% 1. Скорости и база. По умолчанию speed0 = 4 м/с, sensorGap = 0.5 м,
%    speedMin = 0.3, speedMax = 3.0. Окна длиной до 1.75 с при этих
%    настройках физически непригодны (см. ниже), а на 4 м/с поиск
%    упирается в верхнюю границу speedMax = 3 м/с. Поэтому перебор идёт
%    на скоростях 2.0...2.5 м/с: они лежат в (speedMin, speedMax) и
%    ускорение 0.5 м/с^2 за 4 с выводит вагон к границе диапазона только
%    к концу записи, а не с первых окон.
% 2. Ускорение одно (0.5 м/с^2): иначе перебор не укладывается в
%    разумное время. Каждый прогон занимает 5...15 с.
%
% Физическая нижняя граница длины окна:
%   окно должно быть длиннее максимальной задержки + требуемое перекрытие:
%   lagMax = L/(speedMin*Ts) = 0.5/(0.3*0.01) = 166 отсч., minOverlap = 30,
%   итого >= 197 отсч. = 1.97 с. Окна короче функция отвергает ошибкой
%   accComp:windowTooShort (раньше она молча выдавала ноль оценок).
minWin = (floor(0.5/(0.3*0.010)) + 30 + 1) * 0.010;
fprintf('Физическая нижняя граница окна: %.2f с\n', minWin);
winList = winList(winList >= minWin);
fprintf('Длины окон: %s с\n', mat2str(winList));
fprintf('Перекрытия (шаг/окно): %s\n', mat2str(hopFracList));
fprintf('Повторов на точку: %d x (v0,a0) = %d прогонов\n\n', ...
    nSeed, numel(v0List)*numel(a0List)*nSeed);

%% --- Основной перебор ----------------------------------------------------
nW = numel(winList);
nF = numel(hopFracList);

rmseV    = nan(nW, nF);
biasV    = nan(nW, nF);
stdV     = nan(nW, nF);
rmseA    = nan(nW, nF);
stdA     = nan(nW, nF);
coverage = nan(nW, nF);
nWinAvg  = nan(nW, nF);
resSec   = nan(nW, nF);   % временное разрешение = шаг скольжения

tStart = tic;
for iW = 1:nW
    W = winList(iW);
    for iF = 1:nF
        hop = W * hopFracList(iF);
        resSec(iW, iF) = hop;

        eV = []; eA = []; cv = []; nw = [];
        nFailed = 0;   % прогонов, где ни одно окно не было принято
        for v0 = v0List
            for a0 = a0List
                for sd = 1:nSeed
                    o = struct('windowMode', 'fixed', 'winSeconds', W, ...
                               'hopSeconds', hop, 'speed0', v0, ...
                               'accel0', a0, 'randomSeed', sd, ...
                               'makePlot', false);
                    % Окно короче физического минимума отвергается
                    % функцией с ошибкой accComp:windowTooShort - в
                    % переборе такие точки пропускаются, а не роняют
                    % весь прогон.
                    try
                        r = acceleration_compensation(profileFile, o);
                    catch ME
                        if strcmp(ME.identifier, 'accComp:windowTooShort')
                            nFailed = nFailed + 1;
                            continue
                        end
                        rethrow(ME);
                    end
                    ok = r.table.Accepted;
                    if ~any(ok)
                        nFailed = nFailed + 1;
                        continue
                    end

                    eV = [eV; r.table.SpeedErr_mps(ok)];              %#ok<AGROW>
                    fin = ok & isfinite(r.table.Acc_mps2);
                    if any(fin)
                        eA = [eA; r.table.Acc_mps2(fin) - ...
                                  r.table.AccTrue_mps2(fin)];         %#ok<AGROW>
                    end

                    span = sum(r.table.End_s(ok) - r.table.Start_s(ok));
                    cv = [cv; min(span / r.settings.T, 1)];           %#ok<AGROW>
                    nw = [nw; nnz(ok)];                               %#ok<AGROW>
                end
            end
        end

        if ~isempty(eV)
            rmseV(iW,iF) = sqrt(mean(eV.^2));
            biasV(iW,iF) = mean(eV);
            stdV(iW,iF)  = std(eV);
        end
        if ~isempty(eA)
            rmseA(iW,iF) = sqrt(mean(eA.^2));
            stdA(iW,iF)  = std(eA);
        end
        coverage(iW,iF) = mean(cv);
        nWinAvg(iW,iF)  = mean(nw);

        fprintf('W=%4.2f с hop=%4.2f с: RMSE v=%.4f bias=%+.4f | RMSE a=%.3f std a=%.3f | покрытие %3.0f%% | окон %4.1f\n', ...
            W, hop, rmseV(iW,iF), biasV(iW,iF), rmseA(iW,iF), stdA(iW,iF), ...
            100*coverage(iW,iF), nWinAvg(iW,iF));
    end
end
fprintf('\nПеребор завершён за %.1f с\n\n', toc(tStart));

%% --- Сводный критерий и оптимум -----------------------------------------
bestV = min(rmseV(:));
bestA = min(rmseA(:));
score = (rmseV ./ bestV) + (rmseA ./ bestA);
score(coverage < 0.90) = NaN;            % нужно почти полное покрытие

[sBest, iBest] = min(score(:));
[iWb, iFb] = ind2sub(size(score), iBest);

fprintf('=== Оптимум (покрытие >= 90%%) ===\n');
fprintf('Длина окна = %.2f с, шаг %.2f с (перекрытие %.0f%%)\n', ...
    winList(iWb), resSec(iWb,iFb), 100*(1-hopFracList(iFb)));
fprintf('RMSE скорости = %.4f м/с, RMSE ускорения = %.3f м/с^2\n\n', ...
    rmseV(iWb,iFb), rmseA(iWb,iFb));

%% --- Рекомендуемый интервал ---------------------------------------------
% "Хорошие" окна: балл не хуже 1.5 от лучшего при покрытии >= 90%.
threshold = 1.5;
good = score <= sBest * threshold;
goodWin = winList(any(good, 2));

if isempty(goodWin)
    recWin = [winList(iWb) winList(iWb)];
else
    recWin = [min(goodWin) max(goodWin)];
end

fprintf('=== РЕКОМЕНДУЕМЫЙ ИНТЕРВАЛ ДЛИНЫ ОКНА ===\n');
fprintf('%.2f .. %.2f с\n\n', recWin(1), recWin(2));

%% --- Печать таблиц -------------------------------------------------------
fprintf('%-9s %-9s %10s %10s %10s %10s %9s\n', ...
    'W, с', 'перекр.', 'RMSE v', 'bias v', 'RMSE a', 'std a', 'покрытие');
for iW = 1:nW
    for iF = 1:nF
        if isnan(rmseV(iW,iF))
            continue
        end
        fprintf('%-9.2f %8.0f%% %10.4f %+10.4f %10.4f %10.4f %8.0f%%\n', ...
            winList(iW), 100*(1-hopFracList(iF)), rmseV(iW,iF), ...
            biasV(iW,iF), rmseA(iW,iF), stdA(iW,iF), 100*coverage(iW,iF));
    end
end

%% --- Сохранение ----------------------------------------------------------
save(fullfile(outDir, 'window_study.mat'), 'winList', 'hopFracList', ...
    'rmseV', 'biasV', 'stdV', 'rmseA', 'stdA', 'coverage', 'nWinAvg', ...
    'resSec', 'score', 'iWb', 'iFb', 'recWin', 'v0List', 'a0List', 'nSeed');

% CSV: одна строка на каждую комбинацию (длина окна, перекрытие).
% Матрицы метрик заполнялись как rmseV(iW, iF), поэтому линеаризация (:) идёт
% по столбцам: сначала все длины окон при первом перекрытии, затем при втором
% и т.д. Индексы через ndgrid восстанавливают ту же последовательность,
% что гарантирует согласованность строк.
[iWc, iFc] = ndgrid(1:nW, 1:nF);
Wcol = winList(iWc);  Wcol = Wcol(:);
Ocol = 100 * (1 - hopFracList(iFc));  Ocol = Ocol(:);

Tcsv = table(Wcol, Ocol, rmseV(:), biasV(:), stdV(:), rmseA(:), stdA(:), ...
             coverage(:), ...
             'VariableNames', {'Window_s', 'Overlap_pct', 'RMSE_v', ...
                               'Bias_v', 'Std_v', 'RMSE_a', 'Std_a', 'Coverage'});
Tcsv = Tcsv(~isnan(Tcsv.RMSE_v), :);
writetable(Tcsv, fullfile(outDir, 'window_study.csv'));

%% --- Графики -------------------------------------------------------------
fig = figure('Name', 'Точность vs длина окна и перекрытие', ...
    'NumberTitle', 'off', 'WindowStyle', 'normal', ...
    'Units', 'normalized', 'Position', [0.03 0.05 0.95 0.85]);

styles = {'-o', '-s', '-^'};
colors = {[0.00 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19]};

% (а) RMSE скорости
subplot(2, 2, 1); hold on; grid on;
for iF = 1:nF
    plot(winList, rmseV(:,iF), styles{iF}, 'Color', colors{iF}, ...
        'LineWidth', 1.8, 'MarkerSize', 7, ...
        'DisplayName', sprintf('перекрытие %.0f%%', 100*(1-hopFracList(iF))));
end
xlabel('Длина окна, с'); ylabel('RMSE скорости, м/с');
title('(а) Точность скорости');
legend('Location', 'northeast');

% (б) RMSE ускорения
subplot(2, 2, 2); hold on; grid on;
for iF = 1:nF
    plot(winList, rmseA(:,iF), styles{iF}, 'Color', colors{iF}, ...
        'LineWidth', 1.8, 'MarkerSize', 7, ...
        'DisplayName', sprintf('перекрытие %.0f%%', 100*(1-hopFracList(iF))));
end
xlabel('Длина окна, с'); ylabel('RMSE ускорения, м/с^2');
title('(б) Точность ускорения');
legend('Location', 'northeast');

% (в) Сводный балл с зоной оптимума
subplot(2, 2, 3); hold on; grid on;
for iF = 1:nF
    plot(winList, score(:,iF), styles{iF}, 'Color', colors{iF}, ...
        'LineWidth', 1.8, 'MarkerSize', 7, ...
        'DisplayName', sprintf('перекрытие %.0f%%', 100*(1-hopFracList(iF))));
end
patch([recWin(1) recWin(2) recWin(2) recWin(1)], ...
      [0 0 100 100], [0.85 0.95 0.85], 'EdgeColor', 'none', ...
      'FaceAlpha', 0.5, 'DisplayName', 'оптимальный интервал');
yline(threshold, 'k:', 'LineWidth', 1.0, 'DisplayName', 'порог');
uistack(findobj(gca, 'Type', 'line'), 'top');
xlabel('Длина окна, с'); ylabel('Сводный балл (меньше = лучше)');
title(sprintf('(в) Сводный критерий; оптимум %.2f с', winList(iWb)));
legend('Location', 'northeast');
ylim([0 min(max(score(:)), 10)]);

% (г) Разрешение по времени
subplot(2, 2, 4); hold on; grid on;
for iF = 1:nF
    plot(winList, nWinAvg(:,iF), styles{iF}, 'Color', colors{iF}, ...
        'LineWidth', 1.8, 'MarkerSize', 7, ...
        'DisplayName', sprintf('перекрытие %.0f%%', 100*(1-hopFracList(iF))));
end
xlabel('Длина окна, с'); ylabel('Оценок за проезд, шт.');
title('(г) Число оценок на 5-секундный проезд');
legend('Location', 'northeast');

sgtitle(sprintf(['Точность оценки от длины окна и перекрытия ' ...
    '(L = 2.0 м, %d сценариев x %d seed)'], ...
    numel(v0List)*numel(a0List), nSeed));

pngName = fullfile(outDir, 'accuracy_vs_window.png');
exportgraphics(fig, pngName, 'Resolution', 150);
fprintf('\nГрафик: %s\n', pngName);

%% --- 3D-визуализация: поверхность точности -------------------------------
% 2D-графики выше отвечают на вопрос "как зависит от длины окна" для каждого
% перекрытия отдельно. 3D-поверхность показывает обе переменные ОДНОВРЕМЕННО
% и позволяет сразу увидеть общую форму зависимости: есть ли оптимум внутри
% области, насколько плоска "плато" около оптимума, где начинается обрыв.
%
% Оси построены так, чтобы читались физически:
%   X - длина окна, с
%   Y - перекрытие, % (0 = окна встык, 75 = сильное перекрытие)
%   Z - RMSE (скорости или ускорения)

% Сетка для поверхности. Матрицы rmseV/rmseA уже имеют размер nW x nF,
% поэтому для surf годятся как есть: X - длина окна (столбцы), Y - перекрытие.
[Wgrid, Fgrid] = meshgrid(winList, 100*(1-hopFracList));

fig3 = figure('Name', '3D: поверхность точности', 'NumberTitle', 'off', ...
    'WindowStyle', 'normal', 'Units', 'normalized', ...
    'Position', [0.02 0.04 0.96 0.88]);

% ---- (a) Поверхность RMSE скорости с наложением экспериментальных точек
ax1 = subplot(2, 2, 1);
s1 = surf(ax1, Wgrid, Fgrid, rmseV.', 'EdgeColor', [0.35 0.35 0.35], ...
    'FaceAlpha', 0.85);
s1.FaceColor = 'interp';
shading(ax1, 'interp');
colormap(ax1, 'parula');
hold(ax1, 'on');
% Точки данных - чтобы видеть фактическое положение измерений
scatter3(ax1, Wgrid(:), Fgrid(:), rmseV(:), 45, rmseV(:), 'filled', ...
    'MarkerEdgeColor', 'k');
% Оптимум
if iFb > 0
    plot3(ax1, winList(iWb), 100*(1-hopFracList(iFb)), rmseV(iWb,iFb), ...
        'p', 'MarkerSize', 18, 'MarkerFaceColor', [1 0.85 0.1], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
end
grid(ax1, 'on'); box(ax1, 'on');
xlabel(ax1, 'Длина окна, с');
ylabel(ax1, 'Перекрытие, %');
zlabel(ax1, 'RMSE скорости, м/с');
title(ax1, '(а) Точность скорости');
view(ax1, -37.5, 30);
cb1 = colorbar(ax1); cb1.Label.String = 'RMSE, м/с';

% ---- (б) Поверхность RMSE ускорения (логарифмическая шкала по Z)
% Диапазон 0.017..1.5 м/с^2 - это почти два порядка, линейная шкала
% делает лучшие зоны плоскими. Lg-шкала показывает относительные различия.
ax2 = subplot(2, 2, 2);
Zlog = log10(max(rmseA, 1e-3));
s2 = surf(ax2, Wgrid, Fgrid, Zlog.', 'EdgeColor', [0.35 0.35 0.35], ...
    'FaceAlpha', 0.85);
s2.FaceColor = 'interp';
shading(ax2, 'interp');
colormap(ax2, 'parula');
hold(ax2, 'on');
scatter3(ax2, Wgrid(:), Fgrid(:), Zlog(:), 45, Zlog(:), 'filled', ...
    'MarkerEdgeColor', 'k');
plot3(ax2, winList(iWb), 100*(1-hopFracList(iFb)), Zlog(iWb,iFb), ...
    'p', 'MarkerSize', 18, 'MarkerFaceColor', [1 0.85 0.1], ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
grid(ax2, 'on'); box(ax2, 'on');
xlabel(ax2, 'Длина окна, с');
ylabel(ax2, 'Перекрытие, %');
zlabel(ax2, 'log_{10}(RMSE ускорения)');
title(ax2, '(б) Точность ускорения (лог. шкала)');
view(ax2, -37.5, 30);

% ---- (в) Совмещённые поверхности: компромисс скорость/ускорение
% Нормированные RMSE на одних осях: видно, что ускорение - узкое место.
ax3 = subplot(2, 2, 3);
sn = rmseV ./ min(rmseV(:));
sa = rmseA ./ min(rmseA(:));
surf(ax3, Wgrid, Fgrid, sn.', 'EdgeColor', 'none', 'FaceAlpha', 0.45);
hold(ax3, 'on');
surf(ax3, Wgrid, Fgrid, sa.' + 0.02, 'EdgeColor', 'none', 'FaceAlpha', 0.45);
colormap(ax3, 'parula');
grid(ax3, 'on'); box(ax3, 'on');
xlabel(ax3, 'Длина окна, с');
ylabel(ax3, 'Перекрытие, %');
zlabel(ax3, 'Нормированная ошибка (к лучшей)');
title(ax3, '(в) Скорость vs ускорение');
view(ax3, -37.5, 25);
legend(ax3, {'скорость', 'ускорение'}, 'Location', 'northeast');

% ---- (г) Контурная карта (contour3): изолинии точности как рельеф
% Удобно для чтения "где плато": линии уровня позволяют без поворота
% определить границы оптимальной области.
ax4 = subplot(2, 2, 4);
[C, h] = contour3(ax4, Wgrid, Fgrid, rmseV.', 10, 'LineWidth', 1.6);
clabel(C, h, 'FontSize', 8, 'Color', [0.1 0.1 0.1]);
hold(ax4, 'on');
% Зона рекомендуемого интервала на плоскости основания
patch(ax4, [recWin(1) recWin(2) recWin(2) recWin(1)], [0 0 75 75], ...
    zeros(1,4), [0.2 0.75 0.3], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot3(ax4, winList(iWb), 100*(1-hopFracList(iFb)), rmseV(iWb,iFb), ...
    'p', 'MarkerSize', 18, 'MarkerFaceColor', [1 0.85 0.1], ...
    'MarkerEdgeColor', 'k');
grid(ax4, 'on'); box(ax4, 'on');
xlabel(ax4, 'Длина окна, с');
ylabel(ax4, 'Перекрытие, %');
zlabel(ax4, 'RMSE скорости, м/с');
title(ax4, '(г) Изолинии точности скорости');
view(ax4, -37.5, 35);

sgtitle(fig3, sprintf(['Точность как функция (длина окна, перекрытие). ' ...
    'Звезда - оптимум: окно %.2f с, перекрытие %.0f%%'], ...
    winList(iWb), 100*(1-hopFracList(iFb))));

png3 = fullfile(outDir, 'accuracy_vs_window_3d.png');
exportgraphics(fig3, png3, 'Resolution', 150);
fprintf('3D-график: %s\n', png3);

%% --- 3D-график: траектория "длина окна -> точность" в пространстве ------
% Отдельный вид: каждая точка эксперимента - в 3D пространстве
% (длина окна, перекрытие, точность), соединённая траекториями по перекрытию.
% Показывает, что кривые почти не расходятся - перекрытие вторично.

fig3b = figure('Name', '3D: траектории точности', 'NumberTitle', 'off', ...
    'WindowStyle', 'normal', 'Units', 'normalized', ...
    'Position', [0.08 0.12 0.72 0.68]);
ax = axes(fig3b);
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

for iF = 1:nF
    Y = repmat(100*(1-hopFracList(iF)), nW, 1);
    plot3(ax, winList.', Y, rmseV(:,iF), styles{iF}, ...
        'Color', colors{iF}, 'LineWidth', 2.0, 'MarkerSize', 8, ...
        'DisplayName', sprintf('перекрытие %.0f%%', 100*(1-hopFracList(iF))));
end

% Оптимум
plot3(ax, winList(iWb), 100*(1-hopFracList(iFb)), rmseV(iWb,iFb), ...
    'p', 'MarkerSize', 20, 'MarkerFaceColor', [1 0.85 0.1], ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'оптимум');

% Вертикальные "нити" до плоскости основания - помогают читать высоту
for iW = 1:nW
    for iF = 1:nF
        if isnan(rmseV(iW,iF))
            continue
        end
        plot3(ax, [winList(iW) winList(iW)], ...
            [100*(1-hopFracList(iF)) 100*(1-hopFracList(iF))], ...
            [0 rmseV(iW,iF)], '-', 'Color', [0.75 0.75 0.75], ...
            'LineWidth', 0.6, 'HandleVisibility', 'off');
    end
end

xlabel(ax, 'Длина окна, с');
ylabel(ax, 'Перекрытие, %');
zlabel(ax, 'RMSE скорости, м/с');
title(ax, 'Траектории точности в пространстве (окно, перекрытие, ошибка)');
legend(ax, 'Location', 'northeast');
view(ax, -42, 22);

png3b = fullfile(outDir, 'accuracy_trajectories_3d.png');
exportgraphics(fig3b, png3b, 'Resolution', 150);
fprintf('3D-траектории: %s\n', png3b);

fprintf('\nТаблица: %s\n', fullfile(outDir, 'window_study.csv'));
fprintf('Готово.\n');


