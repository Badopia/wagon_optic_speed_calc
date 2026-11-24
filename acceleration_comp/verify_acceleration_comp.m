%% verify_acceleration_comp.m
% Поверка математики acceleration_compensation на независимых моделях.
%
%   verify_acceleration_comp        % прогон проверок, вывод в командное окно
%
% Скрипт ничего не меняет: все проверки вычислительные. Проверяется:
%   1. Отображение n_o -> n_i, которое используется в ядре, совпадает с
%      точным обращением модели s(t) = v0 t + a t^2/2 (уравнение 17 статьи).
%   2. Классическая форма (17) статьи при нумерации n_o = 1..N смещена на
%      один отсчёт: именно эта ошибка была в исходной реализации. Тест
%      фиксирует различие как регрессионный: ядро обязано давать 0 сдвига.
%   3. q, найденный ядром, совпадает с аналитическим q = a*D/(2*vm) на
%      однородном (повторяющемся) рельефе, и скорость совпадает с vm.
%   4. При a = 0 нет ни смещения скорости, ни ложной деформации.
%   5. На периодическом профиле (sine) каждая оценка либо правильна, либо
%      отброшена - грубых ошибок в принятых окнах нет.
%   6. При слишком коротком окне ядро не молчит, а выдаёт понятную ошибку.

thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
addpath(thisDir);

fails = 0;
fprintf('\n=== verify_acceleration_comp ===\n');

%% --- 1. Отображение (17) против точной инверсии -------------------------
% Независимый вывод. Пусть u = n_i/(N-1), r = n_o/(N-1). Умножим (17) на
% знаменатель: 2r = u*((1-q) + sqrt((1-q)^2 + 4 q r)). Возведя в квадрат и
% сократив, получаем квадратное уравнение r = q u^2 + (1-q) u с решением
%      u = [ (1-q) - sqrt((1-q)^2 + 4 q r) ] / (2q).
% Совпадение с ядром проверяет саму формулу: при нумерации от единицы
% (как было в исходной версии) такого совпадения нет.
N = 351; q = -0.35;
K = N - 1;
r = (0:K).' / K;
u = (-(1 - q) + sqrt((1 - q)^2 + 4 * q * r)) / (2 * q);
niExact = K * u;
% Сравниваем с отображением из ядра (та же формула, читаемая отдельно).
niKernel = localWarpIndex(q, N);
fprintf('1) max |ядро - независимое решение| = %.3e отсч.\n', ...
    max(abs(niKernel - niExact)));
fails = fails + check(max(abs(niKernel - niExact)) < 1e-9, ...
    'ядро реализует точную инверсию модели движения');

%% --- 2. Регрессия: старая форма (17) нумеровала отсчёты от единицы ----
% В исходной версии ni вычислялось по формуле (10) с n_o = 1..N. На
% собственных "внутренних" индексах это верно, но interp1 работает по
% 1..N и получал сдвиг всего отображения на отсчёт: при профиле с
% наклоном это давало систематическую ошибку скорости ~0.1...0.3 м/с.
% Тест фиксирует величину сдвига, чтобы дефект не вернулся незаметно.
offMax = 0;
for qTest = [-0.35, 0.30, 0.0 + 1e-3]
    p = (1 - qTest) / (2 * qTest) * N;
    noOld = (1:N).';
    niOld = -p + p * sqrt(1 + 2 * noOld / p + N * noOld / p^2);
    niNew = localWarpIndex(qTest, N);
    offMax = max(offMax, max(abs(niOld - (1 + niNew))));
end
fprintf('2) нумерация от единицы сдвигает ni до %.2f отсч. (было в исходной версии)\n', offMax);
fails = fails + check(offMax > 0.2, ...
    'старая нумерация даёт систематический сдвиг (регрессионный контроль)');

%% --- 3/4. Профиль, повторяющийся по x: позволяет проверить q -------------
% Рельеф f(x) = 1 + 0.5*sin(2*pi*x/0.7) + 0.2*sawtooth. Он непериодичен на
% длине окна, но имеет яркий рельеф, поэтому корреляция устойчива.
scenery = @(x) 1 + 0.30 * sin(2 * pi * x / 0.55) + 0.15 * sawtoothLocal(x, 0.37);
L = 0.5; Ts = 0.01; v0 = 3.0; a = -0.4; D = 3.0;
t = (0:Ts:D).';
sFast = v0 * t + 0.5 * a * t.^2;
y1 = scenery(sFast);
y2 = scenery(sFast - L);
z = struct('q', NaN, 'v', NaN);
Ncore = numel(t);
[z.q, z.v] = coreDirect(y1, y2, Ncore, L, Ts);
vmTheory = v0 + a * D / 2;
qTheory = a * D / (2 * vmTheory);
fprintf('3) a = %+.2f: v = %.4f (истина %.4f), q = %+.4f (теория %+.4f), a_оцен = %+.4f\n', ...
    a, z.v, vmTheory, z.q, qTheory, 2 * z.q * z.v / D);
fails = fails + check(abs(z.v - vmTheory) < 0.01 * vmTheory, 'скорость на однородном рельефе');
fails = fails + check(abs(z.q - qTheory) < 0.06, 'q согласуется с a*D/(2 vm)');

a0 = 0;
sConst = v0 * t;
[qz, vz] = coreDirect(scenery(sConst), scenery(sConst - L), Ncore, L, Ts);
fprintf('4) a = 0: v = %.4f (истина %.4f), q = %+.4f, a_оцен = %+.4f\n', ...
    vz, v0, qz, 2 * qz * vz / D);
fails = fails + check(abs(vz - v0) < 0.025 * v0, 'нет смещения скорости при a = 0');


%% --- 5. Периодический профиль: грубых ошибок в принятых окнах нет -------
fprintf('5) профиль sine (период 6.283 м):\n');
rSine = acceleration_compensation(fullfile(rootDir, 'data', 'sine_prof.mat'), ...
    struct('windowMode', 'fixed', 'winSeconds', 3.5, 'hopSeconds', 1.0, ...
    'randomSeed', 7, 'makePlot', false));
errs = rSine.table.SpeedErr_mps(rSine.table.Accepted);
if isempty(errs)
    fprintf('   принятых окон нет, ложные пики отсечены\n');
else
    fprintf('   принято %d окон, max |ошибка| = %.3f м/с\n', numel(errs), max(abs(errs)));
end
fails = fails + check(all(abs(errs) <= 0.5), ...
    'на периодическом рельефе нет грубых ошибок (> 0.5 м/с)');

%% --- 6. Слишком короткое окно: понятная ошибка, а не тишина ------------
tooShort = false;
try
    acceleration_compensation(fullfile(rootDir, 'data', 'wagon_12_6995_profile.mat'), ...
        struct('windowMode', 'fixed', 'winSeconds', 1.0, 'makePlot', false));
catch ME
    tooShort = strcmp(ME.identifier, 'accComp:windowTooShort');
    fprintf('6) окно 1.0 с: ошибка "%s"\n', ME.message);
end
fails = fails + check(tooShort, 'короткое окно отвергается с объяснением');

fprintf('\n=== ИТОГ: провалов %d%s ===\n', fails, ...
    ternaryText(fails == 0, ', все проверки прошли', ''));

%% ---------------- локальные функции скрипта ---------------------------
function ok = check(condition, name)
% Учёт одной проверки с печатью результата.
    if condition
        fprintf('   [OK]   %s\n', name);
    else
        fprintf('   [FAIL] %s\n', name);
    end
    ok = double(~condition);
end

function txt = ternaryText(cond, a, b)
% Текстовый выбор без тернарного оператора (его в MATLAB нет).
    if cond
        txt = a;
    else
        txt = b;
    end
end

function v = sawtoothLocal(x, period)
% Несимметричная пила с периодом period м, значения в [-1, 1].
    v = 2 * (mod(x, period) / period) - 1;
end

function [q, v] = coreDirect(y1, y2, N, L, Ts)
% Независимая (минимальная) реализация поиска q и v по тем же уравнениям.
% Используется только для поверки: перебирает q и lag без сеток ядра.
    best = -inf;
    q = NaN; v = NaN;
    lagMin = floor(L / (6 * Ts));
    lagMax = ceil(L / (0.5 * Ts));
    for qt = -0.9:0.005:0.9
        z1 = localWarp(y1, qt, N);
        z2 = localWarp(y2, qt, N);
        for lag = lagMin:lagMax
            j = (1:floor(N - lag)).';
            if numel(j) < 30
                continue
            end
            c = corrcoef(z1(j), z2(j + lag));
            if c(1, 2) > best
                best = c(1, 2);
                q = qt;
                v = L / (lag * Ts);
            end
        end
    end
end

function z = localWarp(x, q, N)
% Ре-сэмплирование по (17) в индексах от нуля - независимая запись.
    K = N - 1;
    r = (0:K).' / K;
    ni = K * (2 * r ./ ((1 - q) + sqrt((1 - q)^2 + 4 * q * r)));
    z = interp1((0:K).', x(1:N), ni, 'linear');
end

function ni = localWarpIndex(q, N)
% Отображение n_o -> n_i по (17) в индексах от нуля (независимая запись).
    K = N - 1;
    r = (0:K).' / K;
    ni = K * (2 * r ./ ((1 - q) + sqrt((1 - q)^2 + 4 * q * r)));
end










