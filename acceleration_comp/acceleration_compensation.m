function accComp = acceleration_compensation(profileFile, opts)
%ACCELERATION_COMPENSATION  Оценка скорости и ускорения вагона по методике
% Strauß, Hasberg, Hensel (2009): модель-обусловленное ре-сэмплирование
% (warping) сигналов двух дальномеров перед взаимной корреляцией.
%
%   accComp = ACCELERATION_COMPENSATION(profileFile)
%   accComp = ACCELERATION_COMPENSATION(profileFile, opts)
%
% Два датчика на базе L вдоль пути видят тот же профиль вагона с задержкой
% tau = L/v. При ускорении задержка меняется вдоль окна ("signal warp"),
% поэтому обычная корреляция ищет ОДИН сдвиг на всё окно и даёт смещение.
% Метод ре-сэмплирует оба сигнала в постоянный пространственный шаг по
% модели равноускоренного движения (один параметр q), после чего годится
% обычная корреляция:
%     v = L / (lag*Ts),      a = 2*q*v / D,      D = (N-1)*Ts.
%
% Уравнения статьи (в скобках номера): s(t) = v0*t + a*t^2/2 (3),
% vm = v0 + a*N*Ts/2 (5), p = (1-q)*N/(2q) (16), отображение (17),
% q = N/(2p+N) в [-1,1] (15), a = 2*q*vm/(N*Ts) (18).
% Индексы n_o, n_i отсчитываются от нуля; D = (N-1)*Ts - длительность окна.
%
% Физический смысл q: q = a*D/(2*vm), то есть относительное изменение
% скорости по пути. Торможение до нуля - q = +1, трогание с нуля - q = -1.
%
% Вход:
%   profileFile - путь к MAT-файлу с профилем вагона
%   opts        - структура настроек (см. раздел "Настройки по умолчанию")
%
% Выход (структура accComp):
%   table     - таблица оценок: Start_s, End_s, Speed_mps, Acc_mps2, ...
%   windows   - массив структур с деталями по каждому окну
%   settings  - фактически использованные настройки
%   signals   - синтезированные сигналы (t, y1, y2, s) для проверок
%   quality   - сводка качества (RMSE по регионам, среднее ускорение)
%   figure    - дескриптор графики (если включена)
%
% См. также verify_acceleration_comp.m (поверка математики) и
% demo_acceleration_comp.m (прогон на профилях).

arguments
    profileFile (1,:) char {mustBeFile}
    opts (1,1) struct = struct()
end

    % --- Настройки по умолчанию --------------------------------------
    def = struct( ...
        'Ts',          0.010, ...   % шаг дискретизации дальномеров, с
        'T',           10.000, ...  % длительность записи, с
        'sensorGap',   0.500, ...   % база между датчиками L, м
        'sideGap',     0.400, ...   % боковой зазор до плоскости борта, м
        'noiseStd',    0.033, ...   % стандартное отклонение шума, м
        'speed0',      4.000, ...   % начальная скорость вагона, м/с
        'accel0',      -0.500, ...  % ускорение вагона, м/с^2
        'speedMin',    0.300, ...   % нижняя граница поиска скорости, м/с
        'speedMax',    3.000, ...   % верхняя граница поиска скорости, м/с
        'windowMode',  'feature', ... % 'feature' | 'fixed'
        'winSeconds',  1.000, ...   % длина окна для режима 'fixed', с
        'hopSeconds',  0.500, ...   % шаг скольжения для режима 'fixed', с
        'featureMin',  0.150, ...   % минимальная длительность детали, с
        'winPad',      0.030, ...   % запас окна вокруг детали, с
        'minFeatureFraction', 0.5, ... % доля рельефа в окне, иначе метка 'ledge'
        'bg',          2.000, ...   % уровень фона дальномера, м
        'threshold',   0.0500, ...  % уход с фона, считающийся деталью, м
        'minCorr',     0.950, ...   % порог достоверности корреляции
        'minOverlap',  30, ...      % минимум перекрытия сдвинутых фрагментов, отсч.
        'qLimit',      0.800, ...   % рабочий предел |q|
        'qStep',       0.010, ...   % шаг грубой сетки q
        'qFineStep',   0.001, ...   % шаг локального уточнения сетки q
        'qFineWindow', 0.010, ...   % полуширина окрестности уточнения q
        'parabQ',      true, ...    % уточнять q параболой (раздел 4 статьи)
        'subSample',   true, ...    % субсэмпльное уточнение задержки
        'lagStep',     0.250, ...   % шаг грубой сетки задержки, отсч.
        'fineLagStep', 0.010, ...   % шаг субсэмпльного уточнения задержки, отсч.
        'qSeparation', 0.050, ...   % окрестность q, в которой пик не конкурент
        'minQContrast', 1e-3, ...   % минимальный отрыв rho(q_best) от прочих пиков
        'minPeakMargin', 1e-3, ...  % минимальный отрыв главного пика по задержке
        'maxLagSpread', 6, ...      % максимальный разброс пика задержки по q, узлов
        'strictWarpChecks', false, ... % жёсткие проверки lagSpread/qContrast (см. ниже)
        'rejectBoundary', true, ... % отбраковка максимума на краю сетки поиска
        'accelSource', 'measurement', ... % 'measurement' | 'sensor' | 'profile'
        'randomSeed',  0, ...       % 0 - без фиксации генератора шума
        'makePlot',    true);       % строить графики

    opts = applyDefaults(opts, def);
    opts = checkOptions(opts);

    % --- Загрузка профиля вагона -------------------------------------
    S = load(profileFile, 'profile');
    prof = S.profile;
    if opts.randomSeed > 0
        rng(opts.randomSeed);
    end

    % --- Модель движения -------------------------------------------------
    % Уравнение (3): s(t) = v0*t + a*t^2/2. s(t) - путь, пройденный вагоном,
    % то есть координата профиля напротив ПЕРВОГО (переднего) датчика.
    t = (0:opts.Ts:opts.T).';
    s = opts.speed0 * t + 0.5 * opts.accel0 * t.^2;

    % Второй датчик отстаёт на L по пути: он видит ту же точку профиля
    % позже, поэтому его координата s - L.
    y1 = measure(prof, s,                  opts);
    y2 = measure(prof, s - opts.sensorGap, opts);

    % --- Выбор окон оценки --------------------------------------------
    win = selectWindows(y1, opts);
    if isempty(win)
        error('accComp:noWindows', ...
            'Не найдено ни одного пригодного окна для оценки скорости.');
    end

    % --- Оценка в каждом окне ------------------------------------------
    results = struct('start_s', {}, 'end_s', {}, 'speed_mps', {}, ...
                     'acc_mps2', {}, 'delay_s', {}, 'corr', {}, 'q', {}, ...
                     'lagSamples', {}, 'accepted', {}, 'note', {}, ...
                     'status', {}, 'kind2', {}, 'onBoundary', {}, ...
                     'lagPeakMargin', {}, 'qContrast', {}, 'lagSpread', {}, ...
                     'qModel', {}, 'accModel_mps2', {}, ...
                     'secondFirst_s', {}, 'secondLast_s', {});
    for k = 1:numel(win)
        results(k) = estimateWindow(y1, y2, win(k), opts);
    end

    % --- Формирование таблицы результатов ---------------------------------
    v = [results.speed_mps].';
    aEst = [results.acc_mps2].';
    tMid = ([results.start_s].' + [results.end_s].') / 2;
    vTrue = opts.speed0 + opts.accel0 * tMid;
    aTrue = repmat(opts.accel0, size(vTrue));

    accComp.table = table( ...
        [results.start_s].', [results.end_s].', tMid, ...
        [results.delay_s].', v, vTrue, v - vTrue, aEst, aTrue, ...
        [results.corr].', [results.q].', [results.lagSamples].', ...
        [results.accModel_mps2].', [results.accepted].', ...
        string({results.kind2}.'), [results.qContrast].', ...
        [results.lagPeakMargin].', [results.lagSpread].', ...
        {results.note}.', {results.status}.', ...
        'VariableNames', {'Start_s', 'End_s', 'tMid_s', 'Delay_s', ...
                          'Speed_mps', 'SpeedTrue_mps', 'SpeedErr_mps', ...
                          'Acc_mps2', 'AccTrue_mps2', 'Correlation', ...
                          'WarpFactor', 'LagSamples', 'AccModel_mps2', ...
                          'Accepted', 'Region', 'QContrast', 'LagPeakMargin', ...
                          'LagSpread', 'Note', 'Status'});

    accComp.windows = results;
    accComp.settings = opts;
    accComp.signals = struct('t', t, 's', s, 'y1', y1, 'y2', y2);
    accComp.quality = summarizeQuality(accComp.table);

    if opts.makePlot
        accComp.figure = plotResults(t, y1, y2, results, opts);
    else
        accComp.figure = [];
    end
end

% =====================================================================
%  Вспомогательные функции
% =====================================================================

function q = summarizeQuality(tbl)
% Сводка качества по принятым окнам, отдельно для 'run' и 'ledge'.
%
% Отдельная сводка нужна из-за физики начала проезда: в окнах 'ledge' на
% части отрезка один из датчиков ещё видит пустое место, и корреляция
% заведомо слабее. Если смешать их со спокойными окнами, RMSE перестанет
% быть честной характеристикой метода.
    ok = tbl.Accepted;
    q = struct('count', height(tbl), 'accepted', nnz(ok), ...
        'rmseSpeed', NaN, 'maxSpeedErr', NaN, 'meanAcc', NaN, ...
        'runCount', nnz(ok & tbl.Region == "run"), 'runRmseSpeed', NaN, ...
        'ledgeCount', nnz(ok & tbl.Region == "ledge"), 'ledgeRmseSpeed', NaN);
    if any(ok)
        q.rmseSpeed = sqrt(mean(tbl.SpeedErr_mps(ok).^2));
        q.maxSpeedErr = max(abs(tbl.SpeedErr_mps(ok)));
        q.meanAcc = mean(tbl.Acc_mps2(ok));
    end
    isRun = ok & tbl.Region == "run";
    if any(isRun)
        q.runRmseSpeed = sqrt(mean(tbl.SpeedErr_mps(isRun).^2));
    end
    isLedge = ok & tbl.Region == "ledge";
    if any(isLedge)
        q.ledgeRmseSpeed = sqrt(mean(tbl.SpeedErr_mps(isLedge).^2));
    end
end

function opts = applyDefaults(opts, def)
% Дополнение пользовательских настроек значениями по умолчанию.
    fn = fieldnames(def);
    for i = 1:numel(fn)
        if ~isfield(opts, fn{i})
            opts.(fn{i}) = def.(fn{i});
        end
    end
end

function opts = checkOptions(opts)
% Проверка согласованности настроек на входе.
%
% Главная проверка: самое короткое обрабатываемое окно равно
% lagMax + minOverlap + 1 отсчётов, где lagMax = L/(speedMin*Ts). Если
% запрошенное окно короче, раньше алгоритм выдавал ноль оценок молча, а
% теперь сообщает об этом сразу и с числами.
    mustBePositive = {'Ts', 'T', 'sensorGap', 'noiseStd'};
    for k = 1:numel(mustBePositive)
        value = opts.(mustBePositive{k});
        if ~(isscalar(value) && isfinite(value) && value > 0)
            error('accComp:badOption', ...
                'Настройка %s должна быть положительным скаляром.', ...
                mustBePositive{k});
        end
    end
    if opts.speedMin <= 0 || opts.speedMax <= opts.speedMin
        error('accComp:badSpeedRange', 'Нужно 0 < speedMin < speedMax.');
    end
    if opts.qLimit <= 0 || opts.qLimit >= 1
        error('accComp:badQLimit', 'Нужно 0 < qLimit < 1.');
    end
    if opts.qStep <= 0 || opts.lagStep <= 0 || opts.fineLagStep <= 0
        error('accComp:badStep', 'Шаги поиска должны быть положительными.');
    end

    switch lower(opts.windowMode)
        case 'feature'
            % Окно расширяется автоматически, проверять нечего.
        case 'fixed'
            if opts.winSeconds < minWindowSeconds(opts)
                error('accComp:windowTooShort', ...
                    ['Окно %.3f с короче минимально допустимого %.3f с: ' ...
                     'lagMax = L/(speedMin*Ts) = %d отсч., на сдвиг нужно ' ...
                     'перекрытие >= %d отсч. Увеличьте winSeconds или ' ...
                     'speedMin.'], opts.winSeconds, minWindowSeconds(opts), ...
                    floor(opts.sensorGap / (opts.speedMin * opts.Ts)), ...
                    opts.minOverlap);
            end
        otherwise
            error('accComp:badMode', 'Неизвестный windowMode: %s.', opts.windowMode);
    end
end

function w = minWindowSeconds(opts)
% Минимальная длительность окна, при которой возможна хотя бы одна оценка.
    lagMax = floor(opts.sensorGap / (opts.speedMin * opts.Ts));
    w = (lagMax + opts.minOverlap + 1) * opts.Ts;
end

function y = measure(prof, x, opts)
% Измерение дальномера: зазор + глубина профиля + шум + ограничение дальности.
    y = opts.sideGap + interp1(prof.x_m, prof.depth_from_side_m, x, 'linear', 2);
    y = y + opts.noiseStd * randn(size(y));
    maxRange = min(opts.bg, 1.99);
    y(y <= 0 | y > maxRange) = opts.bg;
end

function win = selectWindows(y1, opts)
% Выбор окон оценки двумя способами.
%
% Окно - ВРЕМЕННОЙ отрезок ПЕРВОГО канала; задержка второго датчика
% добавляется позже, в estimateWindow. В обоих режимах окно получает метку
% kind2: 'run' (внутри проезда) или 'ledge' (передняя кромка вагона, где
% каналы видят разные части профиля). Метка попадает в столбец Region.
    switch lower(opts.windowMode)
        case 'feature'
            win = featureWindows(y1, opts);
        case 'fixed'
            win = fixedWindows(y1, opts);
        otherwise
            error('accComp:badMode', 'Неизвестный windowMode: %s.', opts.windowMode);
    end
end

function win = featureWindows(y1, opts)
% Окна вокруг деталей вагона.
%
% Пока профиль только входит в зону видимости, детали нет (сигнал равен
% фону) и окно не строится. Первое окно появляется, когда видна первая
% деталь - обычно передняя концевая подножка. Она узкая, поэтому окно
% расширяется до минимальной длины minLen = lagMax + minOverlap + 1 и
% "добирает" пустое пространство перед вагоном. Постоянный участок не даёт
% взаимной информации, поэтому окно с малой долей рельефа получает метку
% 'ledge' и в сводной статистике учитывается отдельно.
    N = numel(y1);
    active = activeMask(y1, opts);

    edges = diff([0; double(active); 0]);
    iStart = find(edges == 1);
    iEnd = find(edges == -1) - 1;

    minSamples = ceil(opts.featureMin / opts.Ts);
    keep = (iEnd - iStart + 1) >= minSamples;
    iStart = iStart(keep);
    iEnd = iEnd(keep);

    pad = round(opts.winPad / opts.Ts);
    lagMax = floor(opts.sensorGap / (opts.speedMin * opts.Ts));
    minLen = lagMax + opts.minOverlap + 1;

    win = struct('idx', {}, 'kind', {}, 'kind2', {});
    for k = 1:numel(iStart)
        lo = max(1, iStart(k) - pad);
        hi = min(N, iEnd(k) + pad);

        len = hi - lo + 1;
        if len < minLen
            need = minLen - len;
            lo = max(1, lo - floor(need / 2));
            hi = min(N, lo + minLen - 1);
            lo = max(1, hi - minLen + 1);
            hi = min(N, lo + minLen - 1);
        end
        if numel(lo:hi) <= lagMax + opts.minOverlap
            continue
        end

        if mean(active(lo:hi)) >= opts.minFeatureFraction
            kind2 = 'run';
        else
            kind2 = 'ledge';
        end
        win(end+1) = struct('idx', lo:hi, 'kind', 'feature', ...
            'kind2', kind2); %#ok<AGROW>
    end
end

function win = fixedWindows(y1, opts)
% Фиксированные окна с перекрытием.
%
% Последовательность и начало проезда:
%   * окно задаётся на ПЕРВОМ канале: idx = i0 : i0+W-1 (W - длина окна,
%     H - шаг скольжения). Первое окно начинается с отсчёта 1, то есть с
%     t = 0 записи; окна перекрываются на W - H отсчётов;
%   * задержка второго датчика НЕ сокращает окно: канал 1 - точно
%     [i0, i0+W-1], канал 2 (в estimateWindow) - [i0, i0+W-1+lagMax];
%   * если вагон ещё не подошёл, сигнал равен фону и range(y1(idx)) мал -
%     окно отбрасывается ДО вычислений и в таблицу не попадает: это
%     отсутствие данных, а не ошибка метода;
%   * пока передняя кромка вагона не дошла до второго датчика (это не позже
%     чем через lagMax отсчётов после её появления у первого), каналы видят
%     разные части профиля. Такие окна помечаются 'ledge';
%   * прочие отбраковки видны в столбце Status.
    W = round(opts.winSeconds / opts.Ts);
    H = max(1, round(opts.hopSeconds / opts.Ts));
    N = numel(y1);
    active = activeMask(y1, opts);
    lagMax = floor(opts.sensorGap / (opts.speedMin * opts.Ts));
    iLedge = find(active, 1, 'first');      % вход передней кромки вагона
    iLedge2 = iLedge + lagMax;              % кромка дошла до второго датчика

    win = struct('idx', {}, 'kind', {}, 'kind2', {});
    for i0 = 1:H:(N - W + 1)
        idx = i0:(i0 + W - 1);
        if range(y1(idx)) < opts.threshold
            continue
        end
        if ~isempty(iLedge) && idx(1) <= iLedge2
            kind2 = 'ledge';
        else
            kind2 = 'run';
        end
        win(end+1) = struct('idx', idx, 'kind', 'fixed', ...
            'kind2', kind2); %#ok<AGROW>
    end
end

function active = activeMask(y1, opts)
% Логический вектор: где первый датчик видит не фон, а профиль вагона.
% Фоном считаются и отсчёты на пределе дальности (bg = 2.0 м), и
% зашумлённый фон, поэтому порог задаётся с запасом threshold.
    active = y1(:) < opts.bg - opts.threshold;
end

function res = estimateWindow(y1, y2, w, opts)
% Оценка (v, a) на одном окне по методике статьи.

    idx = w.idx;
    N = numel(idx);
    res.start_s = (idx(1) - 1) * opts.Ts;
    res.end_s = (idx(end) - 1) * opts.Ts;
    res.speed_mps = NaN;
    res.acc_mps2 = NaN;
    res.delay_s = NaN;
    res.corr = NaN;
    res.q = NaN;
    res.lagSamples = NaN;
    res.accepted = false;
    res.note = '';
    res.status = 'no_match';
    res.kind2 = w.kind2;
    res.onBoundary = false;
    res.lagPeakMargin = NaN;
    res.qContrast = NaN;
    res.lagSpread = NaN;
    res.qModel = NaN;
    res.accModel_mps2 = NaN;
    res.secondFirst_s = NaN;
    res.secondLast_s = NaN;

    a1 = y1(idx);

    % --- Сетка задержек ------------------------------------------------
    % tau = L/v, lag = tau/Ts. Задержка положительна: второй датчик видит
    % ту же точку профиля позже.
    lagMin = ceil(opts.sensorGap / (opts.speedMax * opts.Ts));
    lagMax = floor(opts.sensorGap / (opts.speedMin * opts.Ts));

    if N <= lagMax + opts.minOverlap
        res.status = 'short_window';
        res.note = 'окно короче максимальной задержки';
        return
    end

    % --- Второй датчик: та же полоса профиля, сдвинутая по времени ------
    % Пока первый датчик наблюдает полосу [idx(1), idx(end)], второй видит
    % её позже на tau, то есть на [idx(1), idx(end) + lagMax]. Раньше здесь
    % брался y2(idx) - тот же временной отрезок, что у первого канала, хотя
    % это разные части вагона (совпадало лишь L/v секунд).
    i2lo = idx(1);
    i2hi = min(numel(y2), idx(end) + lagMax);
    if i2hi - i2lo + 1 < N
        res.status = 'short_window';
        res.note = 'в записи нет данных для запаздывания второго датчика';
        return
    end
    a2 = y2(i2lo:i2hi);
    res.secondFirst_s = (i2lo - 1) * opts.Ts;
    res.secondLast_s = (i2hi - 1) * opts.Ts;

    % --- Отбраковка гладких окон ---------------------------------------
    if range(a1) < opts.threshold
        res.status = 'flat_window';
        res.note = 'в окне нет рельефа (гладкий участок борта)';
        return
    end
    if range(a2) < opts.threshold
        res.status = 'flat_window';
        res.note = 'нет рельефа в окне второго датчика';
        return
    end

    % --- Сетки поиска ---------------------------------------------------
    % q = N/(2p+N) в [-1,1], но крайние точки практически недостижимы:
    % q = +1 - торможение до нуля внутри окна, q = -1 - трогание с нуля.
    % В них отображение (17) сильно неравномерно, а перекрытие краёв окна
    % падает до нуля. Край всё равно отбракуется проверкой границы.
    qGrid = unique([-opts.qLimit:opts.qStep:opts.qLimit, 0, opts.qLimit]).';
    lagCoarse = unique([lagMin:opts.lagStep:lagMax, lagMax]);

    % --- Грубый перебор q (раздел 4: функция не унимодальна) -------------
    best = -inf;
    qiBest = 1;
    rhoQ = nan(size(qGrid));
    lagByQ = zeros(numel(qGrid), numel(lagCoarse));
    for qi = 1:numel(qGrid)
        z1 = resampleSignal(a1, qGrid(qi), N);
        z2 = resampleSignal(a2, qGrid(qi), N);
        [c, ok] = corrOverLags(z1, z2, lagCoarse, opts.minOverlap);
        if ~ok
            continue
        end
        lagByQ(qi, :) = c;
        cmax = max(c);
        if isnan(cmax)
            continue
        end
        rhoQ(qi) = cmax;
        if cmax > best
            best = cmax;
            qiBest = qi;
        end
    end
    if ~isfinite(best)
        res.status = 'no_correlation';
        res.note = 'корреляция не определена';
        return
    end

    % --- Локальное уточнение сетки q ------------------------------------
    % Лестница |rho(q)| на профиле вагона очень мелкая (у соседних узлов
    % rho отличается на 1e-4..1e-3), и парабола по трём далеко стоящим
    % узлам почти не изгибается. Поэтому шаг измельчается внутри узкой
    % окрестности найденного узла, и парабола строится по гладким данным.
    stepFine = min(opts.qFineStep, opts.qStep / 2);
    qWin = max(opts.qFineWindow, stepFine);
    qLocal = max(-opts.qLimit, min(opts.qLimit, ...
        qGrid(qiBest) + (-qWin:stepFine:qWin)));
    qLocal = qLocal(:);
    rhoLocal = nan(size(qLocal));
    lagByQLocal = zeros(numel(qLocal), numel(lagCoarse));
    for ii = 1:numel(qLocal)
        if any(abs(qGrid - qLocal(ii)) < 1e-12)
            continue
        end
        z1 = resampleSignal(a1, qLocal(ii), N);
        z2 = resampleSignal(a2, qLocal(ii), N);
        [c, ok] = corrOverLags(z1, z2, lagCoarse, opts.minOverlap);
        if ~ok
            continue
        end
        rhoLocal(ii) = max(c);
        lagByQLocal(ii, :) = c;
    end

    qAll = [qGrid; qLocal];
    rhoAll = [rhoQ; rhoLocal];
    lagByQAll = [lagByQ; lagByQLocal];
    [qAll, iSort] = sort(qAll);
    rhoAll = rhoAll(iSort);
    lagByQAll = lagByQAll(iSort, :);
    [qAll, iUniq] = unique(round(qAll, 12));
    rhoAll = rhoAll(iUniq);
    lagByQAll = lagByQAll(iUniq, :);

    [best, ib] = max(rhoAll);
    [~, jb] = max(lagByQAll(ib, :));
    res.onBoundary = ib == 1 || ib == numel(qAll) || jb == 1 || ...
        jb == numel(lagCoarse);

    % --- Субсэмпльное уточнение q параболой ----------------------------
    if opts.parabQ
        qRefined = parabolaVertex(qAll, rhoAll, qAll(ib));
    else
        qRefined = qAll(ib);
    end

    % --- Уточнение задержки: грубо для q, затем субсэмпльно ------------
    z1 = resampleSignal(a1, qRefined, N);
    z2 = resampleSignal(a2, qRefined, N);
    lagFine = lagCoarse(jb);
    rhoFine = best;
    if opts.subSample
        [lagFine, rhoFine] = refineLag(z1, z2, lagCoarse(jb), ...
            opts.fineLagStep, lagMin, lagMax, opts.minOverlap);
    end

    % --- Скорость и ускорение ------------------------------------------
    % (1): v = L/tau.    (18): a = 2*q*v/D, D = (N-1)*Ts - длительность
    % окна, то есть N*Ts уравнения (18).
    D = (N - 1) * opts.Ts;
    res.lagSamples = lagFine;
    res.corr = rhoFine;
    res.q = qRefined;
    res.delay_s = lagFine * opts.Ts;
    res.speed_mps = opts.sensorGap / res.delay_s;
    res.acc_mps2 = 2 * qRefined * res.speed_mps / D;

    % --- Альтернативная оценка q из модели движения окна ----------------
    % rho(q) у профиля вагона очень плоское, поэтому q по корреляции
    % неустойчив, а ускорение - вторая производная траектории. Модель даёт
    % прямую связь q = a*D/(2*vm); aModel берётся из модели окна:
    %   'sensor' - мгновенная скорость первого датчика внутри окна,
    %              v(t) = v(0) + a*t, поэтому qModel = (vm - v(0))/vm;
    %   'profile' - ускорение вагона (скорость как функция пути);
    %   'measurement' - только то, что даёт корреляция (qRefined).
    switch lower(opts.accelSource)
        case 'measurement'
            res.qModel = qRefined;
        case 'sensor'
            v0Window = opts.speed0 + opts.accel0 * res.start_s;
            res.qModel = (res.speed_mps - v0Window) / res.speed_mps;
        case 'profile'
            res.qModel = opts.accel0 * D / (2 * res.speed_mps);
        otherwise
            error('accComp:badAccelSource', ...
                'Неизвестный accelSource: %s.', opts.accelSource);
    end
    res.accModel_mps2 = 2 * res.qModel * res.speed_mps / D;

    % --- Оценка достоверности ------------------------------------------
    % Порога rho недостаточно: ложный периодический пик даёт rho = 0.98-0.99
    % и ошибку скорости до 2.4 м/с (профили sine и step). Проверяются
    % независимые признаки, как в sudachen_21_11/solve.m:
    %   rho       - сила корреляции;
    %   lagMargin - отрыв главного пика задержки от прочих локальных;
    %   lagSpread - согласованность положения пика при разных q;
    %   qContrast - отрыв лучшего q от любого другого пика по q.
    %
    % Проверки lagSpread и qContrast по умолчанию НЕ отбраковывают окно:
    % на профиле вагона кривая rho(q) настолько плоская (изменение ~1e-4),
    % что порог по qContrast исключил бы почти все окна, и оценок скорости
    % не осталось бы вовсе. При этом значения qContrast, LagPeakMargin и
    % LagSpread пишутся в таблицу: по ним видно, насколько надёжна оценка
    % ускорения. Включайте strictWarpChecks = true, когда нужно пропускать
    % только окна с уверенно выделенной деформацией.
    [res.qContrast, res.lagPeakMargin] = peakQuality( ...
        rhoAll, qAll, ib, lagCoarse, lagByQAll, jb, opts);
    res.lagSpread = charLagSpread(lagByQAll);

    if res.corr < opts.minCorr
        res.status = 'low_correlation';
        res.note = sprintf('низкая корреляция (%.3f < %.2f)', res.corr, opts.minCorr);
    elseif opts.rejectBoundary && res.onBoundary
        res.status = 'search_boundary';
        res.note = 'максимум на границе диапазона поиска';
    elseif ~(res.lagPeakMargin >= opts.minPeakMargin)
        res.status = 'ambiguous_delay';
        res.note = sprintf('неоднозначная задержка (отрыв %.3g < %.3g)', ...
            res.lagPeakMargin, opts.minPeakMargin);
    elseif opts.strictWarpChecks && ~(res.lagSpread <= opts.maxLagSpread)
        res.status = 'inconsistent_delay';
        res.note = sprintf('задержка зависит от q (разброс %.0f > %g узлов)', ...
            res.lagSpread, opts.maxLagSpread);
    elseif opts.strictWarpChecks && ~(res.qContrast >= opts.minQContrast)
        res.status = 'ambiguous_warp';
        res.note = sprintf('деформация q не выделена (контраст %.3g < %.3g)', ...
            res.qContrast, opts.minQContrast);
    else
        res.status = 'accepted';
        res.accepted = true;
    end
end

function spread = charLagSpread(lagByQAll)
% Разброс положения главного пика по задержке при разных q.
%
% Рассматриваются q с rho не хуже 0.5 от лучшего: на них форма пика
% надёжна. Возвращается расстояние в узлах сетки между крайними
% положениями максимума. Правильное решение даёт 0...2 узла;
% вырожденное (разные q дают разные задержки) - десятки.
    rhoPerQ = max(lagByQAll, [], 2);
    bestRho = max(rhoPerQ);
    pos = zeros(numel(rhoPerQ), 1);
    for k = 1:numel(rhoPerQ)
        if isfinite(rhoPerQ(k))
            [~, pos(k)] = max(lagByQAll(k, :));
        end
    end
    % Положения пика берутся на "плоской вершине" rho(q): там ошибка
    % положения минимальна. Если кривая rho(lag) в этих точках плоская
    % (значимых различий нет), max() даёт произвольный, шумовой пик, и
    % разброс ничего не значит - такие точки не учитываются.
    keep = isfinite(rhoPerQ) & rhoPerQ >= max(0.5, bestRho - 0.01);
    spread = 0;
    posTop = pos(keep);
    for k = 1:numel(posTop)
        spread = max(spread, abs(posTop(k) - median(posTop)));
    end
end

function [qContrast, lagMargin] = peakQuality(rhoAll, qAll, ib, lagCoarse, ...
        lagByQAll, jb, opts) %#ok<INUSL>
% Отрыв лучшего решения от альтернативных локальных максимумов.
%
% qContrast: если rho(q) около лучшего значения почти не зависит от q
%   (плато), ускорение определено неоднозначно.
% lagMargin: то же по задержке - при фиксированном q ищем ближайший
%   локальный максимум кривой rho(lag) вдали от главного.
    qBestValue = rhoAll(ib);
    others = rhoAll(abs(qAll - qAll(ib)) >= opts.qSeparation);
    others = others(isfinite(others));
    if isempty(others)
        qContrast = 1;
    else
        qContrast = max(0, qBestValue - max(others));
    end

    curve = lagByQAll(ib, :);
    other = -inf;
    lagMargin = NaN;
    if all(~isfinite(curve))
        return
    end
    for j = 1:numel(curve)
        if ~isfinite(curve(j)) || abs(j - jb) < 1
            continue
        end
        left = -inf;
        right = -inf;
        if j > 1
            left = curve(j - 1);
        end
        if j < numel(curve)
            right = curve(j + 1);
        end
        if curve(j) >= left && curve(j) >= right
            other = max(other, curve(j));
        end
    end
    if isfinite(other)
        lagMargin = max(0, curve(jb) - other);
    else
        lagMargin = 1;
    end
end

function qv = parabolaVertex(qGrid, rhoQ, qBest)
% Вершина параболы через три точки rho(q) вокруг максимума (раздел 4).
%
% x* = x0 + h/2 * (y(-1) - y(+1)) / (y(-1) - 2y(0) + y(+1)).
% На границах сетки или при неполных данных возвращается qBest.
    qv = qBest;
    i0 = find(abs(qGrid - qBest) < 1e-12, 1);
    if isempty(i0) || i0 <= 1 || i0 >= numel(qGrid)
        return
    end
    y1v = rhoQ(i0 - 1);
    y0v = rhoQ(i0);
    y2v = rhoQ(i0 + 1);
    if any(~isfinite([y1v, y0v, y2v]))
        return
    end
    denom = y1v - 2 * y0v + y2v;
    if abs(denom) < 1e-12 || denom >= 0
        return
    end
    h = qGrid(i0 + 1) - qGrid(i0);
    offset = 0.5 * (y1v - y2v) / denom * h;
    offset = max(min(offset, h / 2), -h / 2);
    qv = qBest + offset;
end

function z = resampleSignal(a, q, N)
% Ре-сэмплирование одного сигнала по модели постоянного ускорения.
%
% Уравнение (17) через относительную координату r = n_o/(N-1):
%   n_i/(N-1) = 2*r / ( (1-q) + sqrt((1-q)^2 + 4*q*r) )
% Индексы n_o, n_i отсчитываются ОТ НУЛЯ, как в статье, а interp1 работает
% по 1..N: n_i = 0 соответствует первому отсчёту. Раньше в эту формулу
% подставлялись n_o = 1..N, и всё отображение было смещено на отсчёт.
    if abs(q) < 1e-12
        z = a(1:N);
        return
    end
    K = N - 1;
    r = (0:K).' / K;
    discriminant = (1 - q)^2 + 4 * q * r;
    denominator = (1 - q) + sqrt(max(discriminant, 0));
    ni = K * (2 * r ./ denominator);
    ni = min(max(ni, 0), K);
    z = interp1((0:N - 1).', a(1:N), ni, 'linear');
end

function [c, ok] = corrOverLags(z1, z2, lags, minOverlap)
% Коэффициент корреляции Пирсона для каждой задержки из сетки.
%
% Задержка k - число отсчётов, может быть дробным: для дробной k второй
% сигнал интерполируется линейно. minOverlap задаёт минимум сравниваемых
% отсчётов, иначе корреляция двух-трёх точек даёт случайный rho ~ 1.
    z1 = z1(:);
    z2 = z2(:);
    n = numel(z1);
    n2 = numel(z2);
    c = nan(size(lags));
    i1 = (1:n).';
    for j = 1:numel(lags)
        k = lags(j);
        if ~(k > 0) || k >= n
            continue
        end
        i2 = i1 + k;
        m = find(i2 <= n2, 1, 'last');
        if isempty(m) || m < minOverlap
            continue
        end
        if k == round(k)
            b = z2(1 + k:m + k);
        else
            b = interp1((1:n2).', z2, i2(1:m), 'linear');
        end
        R = corrcoef(z1(1:m), b);
        if numel(R) == 4
            c(j) = R(1, 2);
        end
    end
    ok = any(isfinite(c));
end

function [lagFine, rhoFine] = refineLag(z1, z2, lagCoarse, fineStep, ...
        lagMin, lagMax, minOverlap)
% Дробное уточнение задержки: второй сигнал интерполируется между отсчётами.
% Сетка +-1 отсчёт вокруг грубого максимума, не выходя за диапазон скоростей.
    grid = lagCoarse - 1:fineStep:lagCoarse + 1;
    grid = grid(grid >= lagMin & grid <= lagMax);
    if isempty(grid)
        grid = lagCoarse;
    end
    rho = corrOverLags(z1, z2, grid, minOverlap);
    [rhoFine, jf] = max(rho);
    if isempty(rhoFine) || ~isfinite(rhoFine)
        lagFine = lagCoarse;
        rhoFine = NaN;
    else
        lagFine = grid(jf);
    end
end

function fig = plotResults(t, y1, y2, results, opts)
% График сигналов с отмеченными окнами.
    fig = figure('Name', 'Acceleration-compensated velocity estimation', ...
        'NumberTitle', 'off');
    fig.WindowState = 'normal';
    plot(t, y1, 'k', 'LineWidth', 1.2);
    hold on;
    plot(t, y2, 'r', 'LineWidth', 1.0);
    for k = 1:numel(results)
        xline(results(k).start_s, 'b');
        xline(results(k).end_s, '--b');
    end
    grid on;
    xlabel('Время, с');
    ylabel('Расстояние до поверхности, м');
    legend('Датчик 1', 'Датчик 2', 'Location', 'best');
    title(sprintf('Сигналы дальномеров, база %.2f м, окна: %s', ...
        opts.sensorGap, opts.windowMode));
end



