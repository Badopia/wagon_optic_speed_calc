%% demo_acceleration_comp.m
% Демонстрация и поверка методики Strauß 2009 на данных optic_vagon.
%
% Скрипт запускает функцию acceleration_compensation в обоих режимах
% выделения сегментов ('feature' и 'fixed'), сравнивает результат с
% исходным методом main.m (без warping) и строит графики в отдельном окне.
%
% Запуск: из папки optic_vagon, командой  demo_acceleration_comp

clear; close all; clc;

thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);                       % папка optic_vagon
profileFile = fullfile(rootDir, 'data', 'wagon_12_6995_profile.mat');

if ~isfile(profileFile)
    error('Файл профиля не найден: %s', profileFile);
end

addpath(thisDir);                                   % путь к функции

fprintf('Профиль: %s\n\n', profileFile);

%% --- 1. Режим 'feature': окна вокруг деталей вагона ---------------------
% Детали дают мало времени (покрытие ~35%), поэтому для ускорения окно
% искусственно расширяется до минимума, достаточного для корреляции.
opts = struct('windowMode', 'feature', 'randomSeed', 42, 'makePlot', false);
r1 = acceleration_compensation(profileFile, opts);

fprintf('=== Режим feature: окна вокруг деталей ===\n');
disp(summarize(r1.table));
printQuality('feature', r1.table);

%% --- 2. Режим 'fixed': фиксированные окна с перекрытием -----------------
% Окно 2 с выбрано осознанно: ускорение оценивается как a = 2*q*v/(N*Ts),
% поэтому чем короче окно, тем сильнее шум q раздувает разброс a.
% Окно 1 с даёт std(a) ~ 0.5 м/с^2, окно 2 с - ~0.24, окно 3 с - ~0.06.
opts = struct('windowMode', 'fixed', 'winSeconds', 2.0, ...
              'hopSeconds', 0.5, 'randomSeed', 42, 'makePlot', false);
r2 = acceleration_compensation(profileFile, opts);

fprintf('\n=== Режим fixed: окна 2 с, шаг 0.5 с ===\n');
disp(summarize(r2.table));
printQuality('fixed', r2.table);

%% --- 3. Проверка на равномерном движении --------------------------------
opts = struct('windowMode', 'fixed', 'winSeconds', 2.0, 'hopSeconds', 0.5, ...
              'accel0', 0.0, 'randomSeed', 42, 'makePlot', false);
r3 = acceleration_compensation(profileFile, opts);

fprintf('\n=== Контроль: a = 0 (равномерное движение) ===\n');
printQuality('a=0', r3.table);
fprintf('  Ожидание: ошибка скорости ~0, найденный q ~0, a ~0\n');

%% --- 4. Сводка и графики ------------------------------------------------
fprintf('\n=== Итог ===\n');
fprintf('feature: окно %2d..%2d, средняя |ошибка| v = %.3f м/с, a ср. = %.3f м/с^2\n', ...
    min(r1.table.Start_s), max(r1.table.End_s), ...
    mean(abs(r1.table.SpeedErr_mps(r1.table.Accepted))), ...
    mean(r1.table.Acc_mps2(r1.table.Accepted)));
fprintf('fixed:   окно %2d..%2d, средняя |ошибка| v = %.3f м/с, a ср. = %.3f м/с^2\n', ...
    min(r2.table.Start_s), max(r2.table.End_s), ...
    mean(abs(r2.table.SpeedErr_mps(r2.table.Accepted))), ...
    mean(r2.table.Acc_mps2(r2.table.Accepted)));

if any(r2.table.Accepted)
    fh = figure('Name', 'Скорость и ускорение по методике Strauß 2009', ...
        'NumberTitle', 'off', 'WindowStyle', 'normal');
    if ~matlab.ui.internal.isFigureShowEnabled()
        fh.Visible = 'on';
    end

    subplot(2, 1, 1);
    plot(r2.table.tMid_s, r2.table.SpeedTrue_mps, 'k-', 'LineWidth', 1.5); hold on;
    plot(r2.table.tMid_s(r2.table.Accepted), r2.table.Speed_mps(r2.table.Accepted), ...
        'bo', 'MarkerFaceColor', 'b');
    grid on; xlabel('Время, с'); ylabel('Скорость, м/с');
    legend('Истинная', 'Оценка (warp)', 'Location', 'northwest');
    title('Оценка скорости (fixed)');

    subplot(2, 1, 2);
    yline(r2.settings.accel0, 'k-', 'LineWidth', 1.5); hold on;
    plot(r2.table.tMid_s(r2.table.Accepted), r2.table.Acc_mps2(r2.table.Accepted), ...
        'rs', 'MarkerFaceColor', 'r');
    grid on; xlabel('Время, с'); ylabel('Ускорение, м/с^2');
    legend('Истинное', 'Оценка', 'Location', 'best');
    title('Оценка ускорения из параметра q (уравнение 18)');
end

%% --- 4. 3D-визуализация: сравнение режимов и пространство ошибок ---------
if any(r1.table.Accepted) && any(r2.table.Accepted)

    fig3 = figure('Name', '3D: сравнение режимов и точность', ...
        'NumberTitle', 'off', 'WindowStyle', 'normal', ...
        'Units', 'normalized', 'Position', [0.05 0.08 0.9 0.75]);

    % (а) Пространство (t, v, a): два режима на одних осях.
    % Позволяет увидеть, что режим feature даёт редкие, но точные оценки,
    % а fixed - частые; по оси Z сразу виден разброс ускорения.
    ax1 = subplot(1, 2, 1);
    hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');

    tModel = linspace(0, r2.settings.T, 200).';
    vModel = r2.settings.speed0 + r2.settings.accel0 * tModel;
    aModel = repmat(r2.settings.accel0, size(tModel));
    plot3(ax1, tModel, vModel, aModel, 'k-', 'LineWidth', 1.8, ...
        'DisplayName', 'истинная модель');

    ok1 = r1.table.Accepted;
    ok2 = r2.table.Accepted;
    scatter3(ax1, r1.table.tMid_s(ok1), r1.table.Speed_mps(ok1), ...
        r1.table.Acc_mps2(ok1), 100, 'o', 'filled', ...
        'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'k', ...
        'DisplayName', 'feature');
    scatter3(ax1, r2.table.tMid_s(ok2), r2.table.Speed_mps(ok2), ...
        r2.table.Acc_mps2(ok2), 80, 's', 'filled', ...
        'MarkerFaceColor', [0.00 0.45 0.74], 'MarkerEdgeColor', 'k', ...
        'DisplayName', 'fixed');

    xlabel(ax1, 'Время, с'); ylabel(ax1, 'Скорость, м/с');
    zlabel(ax1, 'Ускорение, м/с^2');
    title(ax1, '(а) Оценки в пространстве (t, v, a)');
    legend(ax1, 'Location', 'northeast');
    view(ax1, -40, 22);

    % (б) Столбчатая 3D-диаграмма ошибок по окнам.
    % Три измерения: время окна, тип ошибки, величина. Удобно сравнивать
    % окна между собой и видеть, какая из величин даёт больший вклад.
    ax2 = subplot(1, 2, 2);
    hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');

    idx = find(ok2);
    for j = 1:numel(idx)
        k = idx(j);
        % Столбик ошибки скорости (ось Y = 1) и ускорения (ось Y = 2).
        bar3Offset(ax2, r2.table.tMid_s(k), 1, r2.table.SpeedErr_mps(k), ...
            [0.00 0.45 0.74]);
        bar3Offset(ax2, r2.table.tMid_s(k), 2, ...
            r2.table.Acc_mps2(k) - r2.table.AccTrue_mps2(k), ...
            [0.85 0.33 0.10]);
    end
    plot3(ax2, ax2.XLim, [1 1], [0 0], 'k-', 'LineWidth', 1.0);
    plot3(ax2, ax2.XLim, [2 2], [0 0], 'k-', 'LineWidth', 1.0);

    ax2.YTick = [1 2];
    ax2.YTickLabel = {'ошибка v, м/с', 'ошибка a, м/с^2'};
    xlabel(ax2, 'Время окна, с'); ylabel(ax2, 'Величина');
    zlabel(ax2, 'Ошибка (оценка - истина)');
    title(ax2, '(б) Ошибки по окнам');
    view(ax2, -35, 18);
end

%% ==================== вспомогательные функции ============================
function bar3Offset(ax, x0, y0, h, col)
% Одиночный 3D-столбик в точке (x0, y0) высотой h.
% Используется вместо bar3, потому что bar3 строит целую матрицу столбиков
% и не позволяет свободно чередовать цвета и позиции по осям.
    if ~isfinite(h) || h == 0
        return
    end
    dz = max(abs(h), 1e-6) * 0.02;   % запас, чтобы при нуле не вырождалось
    w = 0.10;                        % полуширина столбика по X
    yb = y0 - 0.25;
    yt = y0 + 0.25;

    % Вершины призмы: нижнее и верхнее основания.
    X = [x0-w, x0+w, x0+w, x0-w, x0-w, x0+w, x0+w, x0-w];
    Y = [yb, yb, yt, yt, yb, yb, yt, yt];
    if h >= 0
        Z = [0 0 0 0, h+dz h+dz h+dz h+dz];
    else
        Z = [h-dz h-dz h-dz h-dz, 0 0 0 0];
    end

    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    patch(ax, 'Vertices', [X(:) Y(:) Z(:)], 'Faces', faces, ...
        'FaceColor', col, 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.5, ...
        'FaceAlpha', 0.9, 'HandleVisibility', 'off');
end

function t = summarize(tbl)
% Краткое представление таблицы результатов.
    t = tbl(:, {'Start_s', 'End_s', 'Speed_mps', 'SpeedTrue_mps', ...
                'SpeedErr_mps', 'Acc_mps2', 'Correlation', ...
                'WarpFactor', 'Accepted'});
end

function printQuality(name, tbl)
% Печать сводных метрик качества.
    ok = tbl.Accepted;
    if ~any(ok)
        fprintf('  %s: ни одного принятого окна\n', name);
        return
    end
    fprintf('  %s: принято %d из %d окон\n', name, nnz(ok), height(tbl));
    fprintf('  средняя |ошибка| скорости = %.3f м/с, макс = %.3f м/с\n', ...
        mean(abs(tbl.SpeedErr_mps(ok))), max(abs(tbl.SpeedErr_mps(ok))));
    if any(isfinite(tbl.Acc_mps2(ok)))
        fprintf('  среднее ускорение = %.3f м/с^2 (истина %.3f)\n', ...
            mean(tbl.Acc_mps2(ok)), tbl.AccTrue_mps2(1));
    end
end
