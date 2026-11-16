function accComp = acceleration_compensation(profileFile, opts)
%ACCELERATION_COMPENSATION  Оценка скорости и ускорения вагона по методике
% Strauß, Hasberg, Hensel (2009): модель-обусловленное ре-сэмплирование
% (warping) сигналов двух дальномеров перед взаимной корреляцией.
%
%   accComp = ACCELERATION_COMPENSATION(profileFile)
%   accComp = ACCELERATION_COMPENSATION(profileFile, opts)
%
% Суть метода (уравнения статьи в скобках):
%   1. Два датчика на базе L видят один и тот же профиль вагона с задержкой.
%   2. При ускорении задержка непостоянна вдоль окна ("signal warp"), поэтому
%      обычная корреляция даёт одну неверную задержку на всё окно.
%   3. Сигналы обоих датчиков ре-сэмплируются в постоянный пространственный
%      шаг по модели равноускоренного движения. Степень деформации описывается
%      ОДНИМ параметром q из [-1, 1] (15), поэтому поиск двумерный (q, lag),
%      но вычислительная сложность остаётся приемлемой.
%   4. После ре-сэмплирования корреляция даёт среднюю скорость vm на окне,
%      а ускорение восстанавливается как a = 2*q*vm/(N*Ts) (18).
%
% Вход:
%   profileFile - путь к MAT-файлу с профилем вагона (data/wagon_*_profile.mat)
%   opts        - структура настроек (см. раздел "Настройки по умолчанию")
%
% Выход (структура accComp):
%   table     - таблица оценок: Start_s, End_s, Speed_mps, Acc_mps2, ...
%   windows   - массив структур с деталями по каждому окну
%   settings  - фактически использованные настройки
%   figure    - дескриптор графики (если включена)
%
% Требования: базовый MATLAB (interp1, corrcoef, table). Тулбоксы не нужны.
%
% См. также MAIN.m - исходная (без warping) реализация для сравнения.

arguments
    profileFile (1,:) char {mustBeFile}
    opts (1,1) struct = struct()
end

    % --- Настройки по умолчанию --------------------------------------
    % Значения подобраны под профиль 12-6995 и физику задачи "проезд
    % датчиков": нужна оценка В ТЕЧЕНИЕ проезда, а не по всему файлу.
    def = struct( ...
        'Ts',          0.010, ...   % шаг дискретизации дальномеров, с
        'T',           5.000, ...   % длительность записи, с
        'sensorGap',   2.000, ...   % база между датчиками L, м  (было 0.50)
        'sideGap',     0.400, ...   % боковой зазор до плоскости борта, м
        'noiseStd',    0.003, ...   % стандартное отклонение шума, м
        'speed0',      5.000, ...   % начальная скорость вагона, м/с
        'accel0',      0.500, ...   % ускорение вагона, м/с^2
        'speedMin',    4.000, ...   % нижняя граница поиска скорости, м/с
        'speedMax',    8.000, ...   % верхняя граница поиска скорости, м/с
        'windowMode',  'feature', ... % 'feature' | 'fixed'
        'winSeconds',  1.000, ...   % длина окна для режима 'fixed', с
        'hopSeconds',  0.500, ...   % шаг скольжения для режима 'fixed', с
        'featureMin',  0.150, ...   % минимальная длительность детали, с
        'winPad',      0.030, ...   % запас окна вокруг детали, с
        'bg',          2.000, ...   % уровень фона дальномера, м
        'threshold',   0.100, ...   % уход с фона, считающийся деталью, м
        'minCorr',     0.950, ...   % порог достоверности корреляции
        'minOverlap',  30, ...      % минимум перекрытия сдвинутых фрагментов, отсч.
        'qLimit',      0.80, ...    % рабочий предел |q| (крайние точки вырождены)
        'qStep',       0.010, ...   % шаг сетки параметра деформации q
        'parabQ',      true, ...    % уточнять q параболой (раздел 4 статьи)
        'subSample',   true, ...    % субсэмпльное уточнение задержки
        'randomSeed',  0, ...       % 0 - без фиксации генератора шума
        'makePlot',    true);       % строить графики

    opts = applyDefaults(opts, def);

    % --- Загрузка профиля вагона -------------------------------------
    % Профиль задан в координате x_m (вдоль вагона) и глубине от плоскости
    % борта. Датчик измеряет расстояние = зазор + глубина.
    S = load(profileFile, 'profile');
    prof = S.profile;
    if opts.randomSeed > 0
        rng(opts.randomSeed);
    end

    % --- Модель движения -------------------------------------------------
    % Уравнение (3) статьи: s(t) = v0*t + a*t^2/2.
    % Здесь s(t) - путь, пройденный вагоном, то есть координата профиля,
    % оказавшаяся напротив ПЕРВОГО (переднего) датчика.
    t = (0:opts.Ts:opts.T).';
    s = opts.speed0 * t + 0.5 * opts.accel0 * t.^2;

    % Второй датчик отстаёт на L по пути: он видит ту же точку профиля
    % позже, поэтому его координата s - L. Знак "минус" соответствует
    % движению в положительном направлении оси x профиля.
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
                     'lagSamples', {}, 'accepted', {}, 'note', {});

    for k = 1:numel(win)
        res = estimateWindow(y1, y2, win(k), opts);
        results(k) = res;
    end

    % --- Формирование таблицы результатов ---------------------------------
    v = [results.speed_mps].';
    aEst = [results.acc_mps2].';
    tMid = ([results.start_s].' + [results.end_s].') / 2;

    % Истинные значения в середине окна - для контроля качества.
    vTrue = opts.speed0 + opts.accel0 * tMid;
    aTrue = repmat(opts.accel0, size(vTrue));

    accComp.table = table( ...
        [results.start_s].', [results.end_s].', tMid, ...
        [results.delay_s].', v, vTrue, v - vTrue, aEst, aTrue, ...
        [results.corr].', [results.q].', [results.lagSamples].', ...
        [results.accepted].', {results.note}.', ...
        'VariableNames', {'Start_s', 'End_s', 'tMid_s', 'Delay_s', ...
                          'Speed_mps', 'SpeedTrue_mps', 'SpeedErr_mps', ...
                          'Acc_mps2', 'AccTrue_mps2', 'Correlation', ...
                          'WarpFactor', 'LagSamples', 'Accepted', 'Note'});

    accComp.windows = results;
    accComp.settings = opts;

    % --- Графика ---------------------------------------------------------
    if opts.makePlot
        accComp.figure = plotResults(t, y1, y2, results, opts);
    else
        accComp.figure = [];
    end
end

% =====================================================================
%  Вспомогательные функции
% =====================================================================

function opts = applyDefaults(opts, def)
% Дополнение пользовательских настроек значениями по умолчанию.
    fn = fieldnames(def);
    for i = 1:numel(fn)
        if ~isfield(opts, fn{i})
            opts.(fn{i}) = def.(fn{i});
        end
    end
end

function y = measure(prof, x, opts)
% Измерение дальномера: зазор + глубина профиля + шум + ограничение дальности.
    % interp1 с 'linear' и extrap=2: вне пределов профиля возврата нет.
    y = opts.sideGap + interp1(prof.x_m, prof.depth_from_side_m, x, 'linear', 2);
    y = y + opts.noiseStd * randn(size(y));

    % Отсчёты вне дальности датчика заменяются уровнем фона.
    maxRange = min(opts.bg, 1.99);
    y(y <= 0 | y > maxRange) = opts.bg;
end

function win = selectWindows(y1, opts)
% Выбор окон оценки двумя способами.
%
% 'feature': окна строятся вокруг деталей вагона (тележки, подножки,
%   резервуар), найденных по сигналу первого датчика. Плюсы: окна
%   содержат ровно тот рельеф, который можно коррелировать. Минусы:
%   покрывают лишь ~35% времени записи.
%
% 'fixed': фиксированные окна с перекрытием. Плюсы: регулярная оценка
%   "в течение проезда" с гарантированным шагом по времени. Минусы:
%   окна попадают на гладкий участок борта и корреляция вырождается.

    switch lower(opts.windowMode)
        case 'feature'
            win = featureWindows(y1, opts);
        case 'fixed'
            win = fixedWindows(y1, opts);
        otherwise
            error('accComp:badMode', ...
                'Неизвестный windowMode: %s', opts.windowMode);
    end
end

function win = featureWindows(y1, opts)
% Окна вокруг деталей вагона.
    % Логический вектор: где датчик видит деталь (уход с фона).
    active = y1(:) < opts.bg - opts.threshold;

    % Переходы 0->1 (начало) и 1->0 (конец) дают границы деталей.
    edges = diff([0; double(active); 0]);
    iStart = find(edges == 1);
    iEnd = find(edges == -1) - 1;

    % Отбрасываем слишком короткие детали: на них корреляция неустойчива.
    minSamples = ceil(opts.featureMin / opts.Ts);
    keep = (iEnd - iStart + 1) >= minSamples;
    iStart = iStart(keep);
    iEnd = iEnd(keep);

    % Запас с каждой стороны, чтобы в окно попали фронты детали.
    pad = round(opts.winPad / opts.Ts);

    % Минимальная длина окна: максимальная задержка плюс требуемое
    % перекрытие плюс один отсчёт запаса (проверка длины строгая).
    % Деталь короче этого не позволит коррелировать сигналы, поэтому окно
    % симметрично расширяется вокруг центра детали.
    lagMax = floor(opts.sensorGap / (opts.speedMin * opts.Ts));
    minLen = lagMax + opts.minOverlap + 1;

    win = struct('idx', {}, 'kind', {});
    for k = 1:numel(iStart)
        lo = max(1, iStart(k) - pad);
        hi = min(numel(y1), iEnd(k) + pad);

        % Расширение до минимальной длины, по возможности симметрично.
        len = hi - lo + 1;
        if len < minLen
            need = minLen - len;
            lo = max(1, lo - floor(need / 2));
            hi = min(numel(y1), lo + minLen - 1);
            lo = max(1, hi - minLen + 1);
        end

        win(end+1) = struct('idx', lo:hi, 'kind', 'feature'); %#ok<AGROW>
    end
end

function win = fixedWindows(y1, opts)
% Фиксированные окна с перекрытием.
    W = round(opts.winSeconds / opts.Ts);
    H = round(opts.hopSeconds / opts.Ts);
    N = numel(y1);

    win = struct('idx', {}, 'kind', {});
    for i0 = 1:H:(N - W + 1)
        idx = i0:(i0 + W - 1);

        % Окно годится, только если в нём есть рельеф: корреляция гладкого
        % сигнала вырождена. Требуем достаточный размах показаний.
        if range(y1(idx)) < opts.threshold
            continue
        end
        win(end+1) = struct('idx', idx, 'kind', 'fixed'); %#ok<AGROW>
    end
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

    a1 = y1(idx);

    % --- Сетка задержек ------------------------------------------------
    % Из скоростного диапазона: tau = L / v, lag = tau / Ts.
    % Задержка всегда положительна: второй датчик "видит" позже.
    lagMin = ceil(opts.sensorGap / (opts.speedMax * opts.Ts));
    lagMax = floor(opts.sensorGap / (opts.speedMin * opts.Ts));
    lags = lagMin:lagMax;

    % Окно должно быть длиннее максимальной задержки с достаточным запасом,
    % иначе коррелировать нечего: перекрытие сдвинутых фрагментов мало.
    % Требуем, чтобы после сдвига оставалось не меньше minOverlap отсчётов.
    if N <= lagMax + opts.minOverlap
        res.note = 'окно короче максимальной задержки';
        return
    end

    % --- Отбраковка гладких окон ---------------------------------------
    % Главный источник ложных оценок: если профиль вагона в окне почти
    % постоянен, корреляция высока при ЛЮБОМ сдвиге, и максимум попадает
    % на край сетки lag (скорость 4 или 8 м/с). Требуем, чтобы в окне
    % был выраженный рельеф.
    if range(a1) < opts.threshold
        res.note = 'в окне нет рельефа (гладкий участок борта)';
        return
    end

    % --- Сетка параметра деформации q (уравнение 15) -------------------
    % q = N^2/(2p + N), теоретически q ∈ [-1, 1]. Но границы вырождены:
    %   q = +1 - разгон из нуля (v0 = 0), q = -1 - торможение в ноль.
    % В этих точках отображение (17) резко неравномерно, а при q -> 1
    % начало окна сжимается так, что перекрытие падает и корреляция
    % становится самокорреляцией (ложный rho -> 1).
    % Поэтому ограничиваем |q| рабочим диапазоном qLimit.
    qGrid = -opts.qLimit:opts.qStep:opts.qLimit;

    best = -inf;
    qBest = 0;
    jBest = 1;
    no = (1:N).';
    rhoQ = nan(size(qGrid));

    for qi = 1:numel(qGrid)
        q = qGrid(qi);

        % Ре-сэмплирование ОБОИХ сигналов одним и тем же отображением.
        % Статья подчёркивает: в отличие от COW/DTW деформация применяется
        % к обоим сигналам независимо, а не подгоняется один под другой.
        [z1, z2] = resamplePair(a1, y2(idx), q, N, no);

        % После ре-сэмплирования сигналы имеют постоянный пространственный
        % шаг, поэтому годится обычная корреляция (уравнение 2).
        [c, ok] = corrOverLags(z1, z2, lags, opts.minOverlap);
        if ~ok
            continue
        end

        [cmax, jb] = max(c, [], 'omitnan');
        if isnan(cmax)
            continue
        end

        % Ложный максимум на краю сетки означает, что корреляция вырождена:
        % сигнал почти гладкий и "похож сам на себя" при любом сдвиге.
        % Такая оценка недостоверна, поэтому край отвергается.
        if jb == 1 || jb == numel(lags)
            continue
        end

        % Запоминаем поверхность rho(q) целиком: она нужна для
        % субсэмпльного уточнения q параболой (раздел 4 статьи, рис. 3).
        rhoQ(qi) = cmax;

        if cmax > best
            best = cmax;
            jBest = jb;
            qBest = q;
        end
    end

    if ~isfinite(best)
        res.note = 'корреляция не определена';
        return
    end

    res.corr = best;
    res.q = qBest;

    % --- Субсэмпльное уточнение q параболой ----------------------------
    % Раздел 4 статьи и рис. 3: зависимость rho(q) вблизи максимума близка
    % к параболе, поэтому три точки сетки вокруг максимума дают вершину
    % с точностью, много лучшей шага сетки.
    %
    % Это критично для ускорения: a = 2*q*vm/(N*Ts) пропорционально q, и
    % при малом N один шаг q (0.01) сдвигает a на 0.1..0.2 м/с^2. Без
    % уточнения оценка ускорения была бы ступенчатой и шумной.
    if opts.parabQ
        qRefined = parabolaVertex(qGrid, rhoQ, qBest);
    else
        qRefined = qBest;
    end

    % --- Уточнение задержки дробным сдвигом ----------------------------
    % Сетка lag дискретна, и v = L/(lag*Ts) даёт ступеньки. Дробное
    % уточнение убирает квантование по задержке.
    [z1, z2] = resamplePair(a1, y2(idx), qRefined, N, no);

    if opts.subSample
        [lagFine, rhoFine] = refineLag(z1, z2, lags(jBest), no, N);
        res.lagSamples = lagFine;
        res.corr = rhoFine;
    else
        res.lagSamples = lags(jBest);
    end

    % --- Скорость из задержки --------------------------------------------
    % Уравнение (1): v = ds/dt = L / tau.
    res.delay_s = res.lagSamples * opts.Ts;
    res.speed_mps = opts.sensorGap / res.delay_s;

    % --- Ускорение из параметра деформации ------------------------------
    % Уравнение (18): a = 2*q*vm/(N*Ts). Оценка следует прямо из найденного
    % q, то есть ускорение измеряется, а не берётся из модели.
    % Используется уточнённый q: см. parabolaVertex выше.
    res.q = qRefined;
    res.acc_mps2 = 2 * qRefined * res.speed_mps / (N * opts.Ts);

    % --- Приёмка по качеству корреляции ----------------------------------
    if res.corr >= opts.minCorr
        res.accepted = true;
    else
        res.accepted = false;
        res.note = sprintf('низкая корреляция (%.3f < %.2f)', ...
            res.corr, opts.minCorr);
    end
end

function qv = parabolaVertex(qGrid, rhoQ, qBest)
% Вершина параболы через три точки rho(q) вокруг максимума.
%
% Раздел 4 статьи: rho(q) вблизи максимума почти параболична, что позволяет
% найти оптимум точнее шага сетки. Вершина параболы через три равноотстоящих
% точки y(-1), y(0), y(+1) есть
%       x* = x0 + h/2 * (y(-1) - y(+1)) / (y(-1) - 2y(0) + y(+1)).
%
% На границах сетки или при неполных данных возвращается qBest.

    qv = qBest;

    i0 = find(abs(qGrid - qBest) < 1e-12, 1);
    if isempty(i0) || i0 <= 1 || i0 >= numel(qGrid)
        % Оптимум на краю сетки: параболу не построить.
        return
    end

    y1v = rhoQ(i0 - 1);
    y0v = rhoQ(i0);
    y2v = rhoQ(i0 + 1);

    % Все три точки должны быть определены (не NaN).
    if any(~isfinite([y1v, y0v, y2v]))
        return
    end

    denom = y1v - 2 * y0v + y2v;

    % Знаменатель близок к нулю - вершина не определена (плоский максимум).
    if abs(denom) < 1e-12
        return
    end

    % Кривизна должна быть отрицательной (максимум, а не минимум).
    if denom >= 0
        return
    end

    h = qGrid(i0 + 1) - qGrid(i0);
    offset = 0.5 * (y1v - y2v) / denom * h;

    % Ограничение сдвига половиной шага сетки: вершина не может уходить
    % от узла больше, чем на половину шага.
    offset = max(min(offset, h / 2), -h / 2);

    qv = qBest + offset;
end

function [z1, z2] = resamplePair(a1, a2, q, N, no)
% Ре-сэмплирование пары сигналов по модели постоянного ускорения.
%
% Вывод уравнения (17): индекс в исходном сигнале ni, соответствующий
% индексу в ре-сэмплированном сигнале no.
%
%   ni = -p + p*sqrt(1 + 2*no/p + N*no/p^2),   p = N*(1-q)/(2q)
%
% Физический смысл: ре-сэмплировать сигнал с постоянным пространственным
% шагом - значит пересчитать его на путь s = vm*no*Ts. Приравнивая пути
% (уравнения 4 и 6), получаем квадратное уравнение (7) относительно ni.

    if abs(q) < 1e-12
        % q = 0: ускорения нет, деформация отсутствует - сигнал как есть.
        z1 = a1;
        z2 = a2;
        return
    end

    if abs(q - 1) < 1e-12
        % q = 1 (разгон из нуля): ветвь уравнения (17) даёт ni = sqrt(N*no).
        ni = sqrt(N * no);
    else
        % Общий случай, уравнение (17).
        p = (1 - q) / (2 * q) * N;
        ni = -p + p * sqrt(max(1 + 2 * no / p + N * no / p.^2, 0));
    end

    % Защита от выхода за границы массива при округлении.
    ni = min(max(ni, 1), N);

    z1 = interp1(no, a1, ni, 'linear');
    z2 = interp1(no, a2, ni, 'linear');
end

function [c, ok] = corrOverLags(z1, z2, lags, minOverlap)
% Коэффициент корреляции Пирсона для каждой задержки из сетки.
%
% minOverlap задаёт минимальное число сравниваемых отсчётов. Без этого
% ограничения при большой задержке перекрытие становится крошечным и
% корреляция двух-трёх точек даёт случайный rho ~ 1.
    c = nan(size(lags));
    for j = 1:numel(lags)
        k = lags(j);
        if numel(z1) - k < minOverlap
            continue
        end
        R = corrcoef(z1(1:end-k), z2(1+k:end));
        if numel(R) == 4
            c(j) = R(1, 2);
        end
    end
    ok = any(~isnan(c));
end

function [lagFine, rhoFine] = refineLag(z1, z2, lagCoarse, no, N)
% Дробное уточнение задержки: второй сигнал интерполируется между отсчётами.
%
% Корреляция считается на мелкой сетке в пределах +-1 отсчёт от грубого
% максимума. Это снимает ступенчатость v = L/(lag*Ts).
    step = 0.01;
    lo = max(lagCoarse - 1, 1);
    hi = lagCoarse + 1;
    grid = lo:step:hi;
    rho = nan(size(grid));

    for j = 1:numel(grid)
        k = grid(j);
        jj = (1:floor(N - k)).';
        if isempty(jj)
            continue
        end
        b = interp1(no, z2, jj + k, 'linear');
        ok = ~isnan(b);
        if nnz(ok) < 5
            continue
        end
        R = corrcoef(z1(jj(ok)), b(ok));
        if numel(R) == 4
            rho(j) = R(1, 2);
        end
    end

    [rhoFine, jf] = max(rho);
    if isempty(rhoFine) || isnan(rhoFine)
        lagFine = lagCoarse;
        rhoFine = NaN;
    else
        lagFine = grid(jf);
    end
end

function fig = plotResults(t, y1, y2, results, opts)
% График сигналов с отмеченными окнами и найденными скоростями.
    fig = figure('Name', 'Acceleration-compensated velocity estimation', ...
        'NumberTitle', 'off');
    fig.WindowState = 'normal';

    subplot(2, 1, 1);
    plot(t, y1, 'k', 'LineWidth', 1.2); hold on;
    plot(t, y2, 'r', 'LineWidth', 1.0);
    for k = 1:numel(results)
        xline(results(k).start_s, ':b');
        xline(results(k).end_s, ':b');
    end
    grid on;
    xlabel('Время, с');
    ylabel('Расстояние до поверхности, м');
    legend('Датчик 1', 'Датчик 2', 'Location', 'best');
    title(sprintf('Сигналы дальномеров, база %.2f м, окна: %s', ...
        opts.sensorGap, opts.windowMode));

    subplot(2, 1, 2);
    tMid = ([results.start_s] + [results.end_s]) / 2;
    vEst = [results.speed_mps];
    aEst = [results.acc_mps2];
    ok = [results.accepted];

    if any(ok)
        plot(tMid(ok), vEst(ok), 'bo', 'MarkerFaceColor', 'b'); hold on;
    end
    if any(~ok)
        plot(tMid(~ok), vEst(~ok), 'ko', 'MarkerFaceColor', 'none');
    end
    grid on;
    xlabel('Время, с');
    ylabel('Скорость, м/с');
    if any(isfinite(aEst))
        title(sprintf('v: точки; a_{сред} = %.3f м/с^2 (модель %.3f)', ...
            mean(aEst(isfinite(aEst))), opts.accel0));
    end

    % --- 3D-представление оценок --------------------------------------
    % Каждая оценка - точка в пространстве (время, скорость, ускорение).
    % Такой вид удобен тем, что связывает обе измеряемые величины: видно,
    % как растёт скорость и одновременно "плывёт" оценка ускорения, а также
    % сразу заметны выбросы по любому из двух измерений.
    if any(ok)
        fig3 = figure('Name', '3D: оценки (t, v, a)', ...
            'NumberTitle', 'off', 'WindowStyle', 'normal');
        ax3 = axes(fig3);
        hold(ax3, 'on'); grid(ax3, 'on'); box(ax3, 'on');

        % Истинная траектория из модели движения: v(t) = v0 + a*t.
        tLine = linspace(0, opts.T, 200).';
        vLine = opts.speed0 + opts.accel0 * tLine;
        aLine = repmat(opts.accel0, size(tLine));
        plot3(ax3, tLine, vLine, aLine, 'k-', 'LineWidth', 1.8, ...
            'DisplayName', 'истинная модель');

        % Принятые оценки, окрашенные по достоверности корреляции.
        vA = vEst(ok); aA = aEst(ok); tA = tMid(ok);
        cA = [results(ok).corr];
        scatter3(ax3, tA, vA, aA, 90, cA, 'filled', ...
            'MarkerEdgeColor', 'k', 'DisplayName', 'принятые оценки');
        cb = colorbar(ax3); cb.Label.String = 'correlation';

        % Отбракованные окна - полыми маркерами (если есть).
        if any(~ok)
            scatter3(ax3, tMid(~ok), vEst(~ok), aEst(~ok), 70, 'k', ...
                'Marker', 'o', 'DisplayName', 'отбракованные');
        end

        xlabel(ax3, 'Время, с');
        ylabel(ax3, 'Скорость, м/с');
        zlabel(ax3, 'Ускорение, м/с^2');
        title(ax3, 'Оценки в пространстве (время, скорость, ускорение)');
        legend(ax3, 'Location', 'northeast');
        view(ax3, -40, 25);
    end
end
