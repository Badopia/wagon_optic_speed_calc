function [output1,output2,accepted] = solve(arg1,arg2,arg3,arg4,arg5,arg6)
%SOLVE Окна с перекрытием, ожидание второго датчика, warping и обновление движения.
%
% Работа с записями или порциями с отдельными временными метками:
%   state=[];
%   [result,state]=solve(t1,y1,t2,y2,analysis,state);
%   result.windows   — границы всех окон и причины отказов;
%   result.estimates — принятые средние скорости и ускорения;
%   result.current   — скорость по последней принятой оценке и её ускорению.
%   Ускорение удерживается, скорость изменяется линейно со временем.
%   analysis.printWindows=true — вывод новых окон и изменений их статуса
%   в командное окно MATLAB; false — отключение этого вывода.
%   При следующем вызове передавайте только новые отсчёты и прежний state.
% Работа в существующем блоке Simulink, один отсчёт на вызов:
%   [v,a,accepted]=solve(y1,y2,valid1,valid2,time,analysis);
%   Здесь обе временные метки равны time. Состояние хранится внутри функции.
%   При новом запуске времени или изменении настроек оно сбрасывается.
%   Также можно сбросить его вручную: clear solve.
%
% analysis создаётся в parametrs.m при ручном запуске main.m.
% Первое окно: [t0,t0+windowSeconds], второе: [t0,t0+windowSeconds+L/vMin].
% Полная обработка расположена ниже в этом файле. Расчёт по паре участков
% выполняется локальной функцией estimate_motion в этом же файле.
% Ничего автоматически не загружается из файлов.
%
% Два интерфейса различаются пятым аргументом:
% структура analysis — порции данных; числовое time — один отсчёт Simulink.

if nargin<6, arg6=[]; end
if isstruct(arg5)
    [output1,output2]=processChunks(arg1,arg2,arg3,arg4,arg5,arg6);
    accepted=output1.current.valid;
    return;
end

persistent stepState stepTime stepConfig
validateattributes(arg1,{'numeric'},{'scalar','real'});
validateattributes(arg2,{'numeric'},{'scalar','real'});
validateattributes(arg5,{'numeric'},{'scalar','real','finite'});
validateattributes(arg3,{'numeric','logical'},{'scalar','real','finite'});
validateattributes(arg4,{'numeric','logical'},{'scalar','real','finite'});
cfg=arg6;
% Совместимость с прежним набором p из старого блока.
if isstruct(cfg) && ~isfield(cfg,'windowSeconds') && isfield(cfg,'windowSamples')
    p=cfg;
    cfg=struct('L',p.L,'Ts',p.Ts,'vMin',p.vMin,'vMax',p.vMax,'bg',p.bg, ...
        'windowSeconds',p.windowSamples*p.Ts, ...
        'overlapFraction',1-p.hopSamples/p.windowSamples, ...
        'motion',struct('minCorr',p.minCorr,'qLimit',p.qLimit,'qStep',p.qStep));
end
cfg=normalizeConfig(cfg);
if isempty(stepTime) || arg5<=stepTime || ~isequaln(cfg,stepConfig)
    stepState=[];
    stepConfig=cfg;
end
value1=double(arg1); value2=double(arg2);
if ~logical(arg3) || ~isfinite(value1), value1=cfg.bg; end
if ~logical(arg4) || ~isfinite(value2), value2=cfg.bg; end
[result,stepState]=processChunks(arg5,value1,arg5,value2,cfg,stepState);
stepTime=arg5;
output1=result.current.speed_mps;
output2=result.current.acc_mps2;
accepted=result.current.valid;
end

function [result,state] = processChunks(t1,y1,t2,y2,config,state)
%PROCESSCHUNKS Полная обработка двух потоков с временными метками.
% [result,state]=solve(t1,y1,t2,y2,config,state)
% Первый вызов: state=[] или шестой аргумент не задан.
% Последующие вызовы: передавайте ТОЛЬКО НОВЫЕ отсчёты и прежний state.
% Можно передать всю запись одним вызовом. Пустые порции разрешены.
% Время — числовые секунды общей шкалы; оба потока строго возрастают.
%
% config: Ts, L, vMin, vMax, windowSeconds, overlapFraction.
% Дополнительные настройки и их значения по умолчанию — normalizeConfig.
% Интервалы: [start,start+windowSeconds] и [start,start+windowSeconds+L/vMin].
% vMin предполагается нижней границей фактической скорости прохода,
% чтобы запас L/vMin соответствовал условию задачи.
% Окно первого канала регистрируется, когда оно полностью поступило.
% Расчёт выполняется только после полного поступления второго интервала.
% Исходные метки могут различаться. Первый канал интерполируется на Ts;
% второй передаётся в estimate_motion с собственными временными метками.
% Большие пропуски и нечисловые показания не интерполируются.
%
% result.windows — все окна, причины отказов, ожидающие и принятые оценки.
% result.estimates — только принятые окна; средняя скорость относится к Center_s.
% result.current: v(t)=vm+a*(t-Center_s); ускорение равно последнему принятому a.
% Принятое окно уточняет опорные vm и a; отказы эту опору не меняют.
% Перекрывающиеся окна не считаются независимыми измерениями и не усредняются.
% До первой принятой оценки выходы NaN. После неё прогноз не обрывается
% по времени; current.referenceTime_s сохраняет время середины исходного окна.
% state хранит накопленные записи и результаты; ничего не пишет в файлы.
%
% Пример после main.m:
% state=[];
% [result,state]=solve(t1,y1,t2,y2,analysis,state);
% disp(result.windows(:,{'Start_s','End1_s','End2_s','Status'}));
% disp(result.estimates(:,{'Center_s','SpeedMean_mps','Acc_mps2'}));

    cfg=normalizeConfig(config);
    if nargin<6 || isempty(state)
        state=struct('config',cfg,'t1',zeros(0,1),'y1',zeros(0,1), ...
            't2',zeros(0,1),'y2',zeros(0,1),'origin',NaN,'nextId',1, ...
            'windows',repmat(windowTemplate(),0,1), ...
            'printedStatuses',strings(0,1),'printedHeader',false);
    else
        assert(isequaln(state.config,cfg), ...
            'При изменении настроек начните новый поток с state=[].');
    end
    [state.t1,state.y1]=appendData(state.t1,state.y1,t1,y1,'1');
    [state.t2,state.y2]=appendData(state.t2,state.y2,t2,y2,'2');
    if isnan(state.origin) && ~isempty(state.t1), state.origin=state.t1(1); end
    duration=round(cfg.windowSeconds/cfg.Ts)*cfg.Ts;
    hop=max(1,round(duration*(1-cfg.overlapFraction)/cfg.Ts))*cfg.Ts;
    tail=cfg.L/cfg.vMin;
    if ~isempty(state.t1)
        % Регистрируем завершённые окна первого канала, даже если второй ещё ждём.
        while true
            start=state.origin+(state.nextId-1)*hop;
            stop=start+duration;
            if state.t1(end)<stop-timeTolerance(stop), break; end
            row=windowTemplate();
            row.Id=state.nextId; row.Start_s=start; row.End1_s=stop;
            row.End2_s=stop+tail; row.Center_s=(start+stop)/2;
            [~,x,reason]=extractInterval(state.t1,state.y1,start,stop,cfg,true);
            if reason=="ok"
                [suitable,q,reason]=quality(x,cfg);
                row.SuitableFirst=suitable;
                row.ValidFraction1=q.fraction; row.Range1_m=q.range;
                row.Std1_m=q.std; row.Changes1=q.changes;
            end
            if row.SuitableFirst
                row.Status="waiting_sensor2";
            else
                row.Status="rejected_first: "+reason;
            end
            state.windows(end+1,1)=row;
            state.nextId=state.nextId+1;
        end
    end

    for k=1:numel(state.windows)
        row=state.windows(k);
        if row.Status~="waiting_sensor2", continue; end
        if isempty(state.t2) || state.t2(end)<row.End2_s-timeTolerance(row.End2_s)
            continue;
        end
        row.Ready=true;
        i=find(state.t2>=row.End2_s-timeTolerance(row.End2_s),1);
        row.Available_s=state.t2(i);
        [secondTime,second,reason]=extractInterval( ...
            state.t2,state.y2,row.Start_s,row.End2_s,cfg,false);
        if reason~="ok"
            row.Status="rejected_second: "+reason;
        else
            [row.SuitableSecond,q,reason]=quality(second,cfg);
            row.ValidFraction2=q.fraction; row.Range2_m=q.range;
            row.Std2_m=q.std; row.Changes2=q.changes;
            if ~row.SuitableSecond
                row.Status="rejected_second: "+reason;
            else
                [~,first,reason]=extractInterval(state.t1,state.y1, ...
                    row.Start_s,row.End1_s,cfg,true);
                assert(reason=="ok");
                motion=cfg.motion;
                motion.sampleTimes2=secondTime-row.Start_s;
                [vm,a,info]=estimate_motion(first,second,cfg.L,cfg.Ts,cfg.vMin,cfg.vMax,motion);
                row.Correlation=info.correlation; row.Q=info.q;
                row.LagPeakMargin=info.lagPeakMargin; row.QContrast=info.qContrast;
                row.OnSearchBoundary=info.onSearchBoundary;
                row.Status="rejected_motion: "+string(info.status);
                if info.accepted && info.correlation>=cfg.minCorrelation
                    if cfg.rejectBoundary && info.onSearchBoundary
                        row.Status="rejected_motion: search_boundary";
                    elseif isnan(info.lagPeakMargin) || info.lagPeakMargin<cfg.minLagPeakMargin
                        row.Status="rejected_motion: ambiguous_delay";
                    elseif ~isfinite(info.qContrast) || info.qContrast<cfg.minQContrast
                        row.Status="rejected_motion: uncertain_acceleration";
                    else
                        row.Accepted=true; row.Status="accepted";
                        row.SpeedMean_mps=vm; row.Acc_mps2=a;
                        row.SpeedAtAvailable_mps=vm+a*(row.Available_s-row.Center_s);
                    end
                elseif info.accepted
                    row.Status="rejected_motion: low_correlation";
                end
            end
        end
        state.windows(k)=row;
    end
    if cfg.printWindows, state=printWindowUpdates(state,cfg); end
    allWindows=struct2table(state.windows);
    acceptedWindows=allWindows(allWindows.Accepted,:);
    now=max([state.t1;state.t2]);
    if isempty(now), now=NaN; end
    current=struct('time_s',now,'valid',false,'speed_mps',NaN,'acc_mps2',NaN, ...
        'referenceTime_s',NaN,'availableTime_s',NaN,'referenceSpeed_mps',NaN, ...
        'windowId',NaN,'status',"no_estimate");
    if ~isempty(acceptedWindows)
        ref=acceptedWindows(end,:);
        current.referenceTime_s=ref.Center_s;
        current.availableTime_s=ref.Available_s;
        current.referenceSpeed_mps=ref.SpeedMean_mps;
        current.windowId=ref.Id;
        current.valid=true;
        current.speed_mps=ref.SpeedMean_mps+ref.Acc_mps2*(now-ref.Center_s);
        current.acc_mps2=ref.Acc_mps2;
        if now>ref.Available_s+timeTolerance(now)
            current.status="predicted_from_last_accepted";
        else
            current.status="accepted_update";
        end
    end
    result=struct('windows',allWindows,'estimates',acceptedWindows,'current',current, ...
        'windowSeconds',duration,'hopSeconds',hop,'tailSeconds',tail);
end

function cfg=normalizeConfig(cfg)
    assert(isstruct(cfg)&&isscalar(cfg),'config должен быть структурой.');
    for name={'Ts','L','vMin','vMax','windowSeconds'}
        assert(isfield(cfg,name{1}),'Не задан config.%s.',name{1});
        validateattributes(cfg.(name{1}),{'double'},{'scalar','real','finite','positive'});
    end
    assert(cfg.vMax>cfg.vMin);
    assert(isfield(cfg,'overlapFraction'));
    validateattributes(cfg.overlapFraction,{'double'},{'scalar','real','finite','>=',0,'<',1});
    assert(round(cfg.windowSeconds/cfg.Ts)>=4,'Окно должно содержать минимум пять отсчётов.');
    defaults=struct('bg',2,'backgroundTolerance',1e-6,'minValidFraction',0.3, ...
        'minValidSamples',10,'minRange',0.03,'minStd',0.01,'minChanges',2, ...
        'changeThreshold',0.02,'smoothSamples',3,'maxGapSeconds',2.5*cfg.Ts, ...
        'minCorrelation',0.9,'minLagPeakMargin',0.001,'minQContrast',0.0001, ...
        'rejectBoundary',true, ...
        'motion',struct(),'printWindows',true);
    for name=fieldnames(defaults).'
        if ~isfield(cfg,name{1}), cfg.(name{1})=defaults.(name{1}); end
    end
    validateattributes(cfg.bg,{'double'},{'scalar','real','finite','positive'});
    validateattributes(cfg.minValidFraction,{'double'},{'scalar','real','finite','>=',0,'<=',1});
    validateattributes(cfg.minCorrelation,{'double'},{'scalar','real','finite','>=',-1,'<=',1});
    for name={'backgroundTolerance','minLagPeakMargin','minQContrast'}
        validateattributes(cfg.(name{1}),{'double'},{'scalar','real','finite','nonnegative'});
    end
    for name={'minRange','minStd','changeThreshold','maxGapSeconds'}
        validateattributes(cfg.(name{1}),{'double'},{'scalar','real','finite','positive'});
    end
    for name={'minValidSamples','minChanges','smoothSamples'}
        validateattributes(cfg.(name{1}),{'double'},{'scalar','real','finite','integer','positive'});
    end
    validateattributes(cfg.rejectBoundary,{'logical'},{'scalar'});
    validateattributes(cfg.printWindows,{'logical'},{'scalar'});
    assert(isstruct(cfg.motion)&&isscalar(cfg.motion));
end

function [times,values]=appendData(times,values,t,y,label)
    validateattributes(t,{'numeric'},{'real','finite'});
    validateattributes(y,{'numeric'},{'real'});
    assert((isvector(t)||isempty(t))&&(isvector(y)||isempty(y))&&numel(t)==numel(y), ...
        'Время и значения датчика %s должны быть векторами одинаковой длины.',label);
    t=double(t(:)); y=double(y(:));
    assert(all(diff(t)>0),'Времена датчика %s должны строго возрастать.',label);
    if ~isempty(t)&&~isempty(times)
        assert(t(1)>times(end),'Передавайте только новые отсчёты датчика %s.',label);
    end
    times=[times;t]; values=[values;y];
end

function [time,x,reason]=extractInterval(t,y,start,stop,cfg,uniform)
    time=[]; x=[]; reason="missing_start";
    tol=timeTolerance(stop);
    if isempty(t)||t(1)>start+tol, return; end
    reason="missing_end";
    if t(end)<stop-tol, return; end
    left=find(t<=start+tol,1,'last'); right=find(t>=stop-tol,1);
    spanT=t(left:right); spanY=y(left:right);
    reason="sample_gap";
    if any(diff(spanT)>cfg.maxGapSeconds), return; end
    reason="nonfinite_data";
    if any(~isfinite(spanY)), return; end
    if uniform
        time=start+(0:round((stop-start)/cfg.Ts)).'*cfg.Ts;
    else
        time=[start;spanT(spanT>start+tol & spanT<stop-tol);stop];
    end
    query=max(spanT(1),min(spanT(end),time));
    x=interp1(spanT,spanY,query,'linear');
    if any(~isfinite(x)), return; end
    reason="ok";
end

function [ok,q,reason]=quality(x,cfg)
    mask=x>0 & abs(x-cfg.bg)>cfg.backgroundTolerance;
    sorted=sort(x); n=numel(x);
    p5=interp1((1:n).',sorted,1+0.05*(n-1));
    p95=interp1((1:n).',sorted,1+0.95*(n-1));
    smooth=movmean(x,min(n,cfg.smoothSamples));
    anchor=smooth(1); changes=0;
    for i=2:n
        if abs(smooth(i)-anchor)>=cfg.changeThreshold
            changes=changes+1; anchor=smooth(i);
        end
    end
    q=struct('fraction',mean(mask),'range',p95-p5,'std',std(x,1),'changes',changes);
    ok=false;
    if nnz(mask)<cfg.minValidSamples || q.fraction<cfg.minValidFraction
        reason="too_few_returns";
    elseif q.range<cfg.minRange || q.std<cfg.minStd
        reason="flat_or_noise";
    elseif changes<cfg.minChanges
        reason="too_few_features";
    else
        ok=true; reason="ok";
    end
end

function row=windowTemplate()
    row=struct('Id',0,'Start_s',NaN,'End1_s',NaN,'End2_s',NaN,'Center_s',NaN, ...
        'SuitableFirst',false,'Ready',false,'SuitableSecond',false,'Accepted',false, ...
        'Status',"waiting_sensor2",'Available_s',NaN,'ValidFraction1',NaN, ...
        'Range1_m',NaN,'Std1_m',NaN,'Changes1',NaN,'ValidFraction2',NaN, ...
        'Range2_m',NaN,'Std2_m',NaN,'Changes2',NaN,'Correlation',NaN, ...
        'LagPeakMargin',NaN,'QContrast',NaN,'Q',NaN,'OnSearchBoundary',false, ...
        'SpeedMean_mps',NaN,'Acc_mps2',NaN,'SpeedAtAvailable_mps',NaN);
end

function tolerance=timeTolerance(time)
    tolerance=max(1e-9,32*eps(max(1,abs(time))));
end

function state=printWindowUpdates(state,cfg)
% Печатаем только новые окна и изменения статуса, не на каждом такте.
    if ~state.printedHeader && ~isempty(state.windows)
        fprintf('\nОКНА АНАЛИЗА: окно %.3f с, шаг %.3f с, добавка датчика 2 %.3f с\n', ...
            cfg.windowSeconds,cfg.windowSeconds*(1-cfg.overlapFraction),cfg.L/cfg.vMin);
        state.printedHeader=true;
    end
    for k=1:numel(state.windows)
        row=state.windows(k);
        if numel(state.printedStatuses)>=k && state.printedStatuses(k)==row.Status, continue; end
        fprintf('#%02d  датчик 1 [%.3f; %.3f], датчик 2 [%.3f; %.3f]: ', ...
            row.Id,row.Start_s,row.End1_s,row.Start_s,row.End2_s);
        if row.Status=="accepted"
            fprintf('ПРИНЯТО | v_ср=%.5f м/с, a=%.5f м/с^2, R=%.6f, q=%.6f, t_ср=%.3f\n', ...
                row.SpeedMean_mps,row.Acc_mps2,row.Correlation,row.Q,row.Center_s);
        elseif row.Status=="waiting_sensor2"
            fprintf('ОЖИДАНИЕ данных второго датчика до %.3f с\n',row.End2_s);
        else
            if startsWith(row.Status,"rejected_first")
                source='датчик 1'; fraction=row.ValidFraction1;
                range=row.Range1_m; deviation=row.Std1_m; changes=row.Changes1;
            else
                source='датчик 2'; fraction=row.ValidFraction2;
                range=row.Range2_m; deviation=row.Std2_m; changes=row.Changes2;
            end
            parts=split(row.Status,': '); code=parts(end);
            switch code
                case "too_few_returns"
                    fprintf('ОТКАЗ (%s): мало возвратов, %.1f%%; нужно >=%.1f%% и >=%d отсчётов\n', ...
                        source,100*fraction,100*cfg.minValidFraction,cfg.minValidSamples);
                case "flat_or_noise"
                    fprintf('ОТКАЗ (%s): мало рельефа, P95-P5=%.5f м (порог %.5f), std=%.5f м (порог %.5f)\n', ...
                        source,range,cfg.minRange,deviation,cfg.minStd);
                case "too_few_features"
                    fprintf('ОТКАЗ (%s): изменений %g, требуется >=%d\n',source,changes,cfg.minChanges);
                case "sample_gap"
                    fprintf('ОТКАЗ (%s): разрыв записи больше %.5f с\n',source,cfg.maxGapSeconds);
                case "nonfinite_data"
                    fprintf('ОТКАЗ (%s): в данных есть NaN или Inf\n',source);
                case {"missing_start","missing_end"}
                    fprintf('ОТКАЗ (%s): не записан полный интервал (%s)\n',source,char(code));
                case "low_correlation"
                    fprintf('ОТКАЗ: низкая корреляция R=%.6f, требуется >=%.6f\n',row.Correlation,cfg.minCorrelation);
                case "search_boundary"
                    fprintf('ОТКАЗ: максимум на границе поиска, q=%.6f, R=%.6f\n',row.Q,row.Correlation);
                case "ambiguous_delay"
                    fprintf('ОТКАЗ: неоднозначная задержка, разница пиков %.6g, требуется >=%.6g\n',row.LagPeakMargin,cfg.minLagPeakMargin);
                case "uncertain_acceleration"
                    fprintf('ОТКАЗ: ускорение не определено однозначно, контраст q %.6g, требуется >=%.6g\n',row.QContrast,cfg.minQContrast);
                otherwise
                    fprintf('ОТКАЗ: %s\n',char(code));
            end
        end
        state.printedStatuses(k,1)=row.Status;
    end
end

function [vMean,acceleration,info]=estimate_motion(window1,data2,L,Ts,vMin,vMax,options)
% Расчёт пары участков методом модельного пересэмплирования Strauß et al.
% D=(N-1)*Ts; q=a*D/(2*vMean); vMean=L/(lag*Ts).
% Второй сигнал может иметь собственные относительные временные метки.
    if nargin<7, options=struct(); end
    validateattributes(window1,{'numeric'},{'vector','real','finite','nonempty'});
    validateattributes(data2,{'numeric'},{'real','finite'});
    assert(isvector(data2)||isempty(data2),'data2 должен быть вектором.');
    validateattributes(L,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(Ts,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(vMin,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(vMax,{'numeric'},{'scalar','real','finite','>',vMin});
    assert(isstruct(options)&&isscalar(options));
    defaults=struct('qLimit',0.8,'qStep',0.01,'lagStep',0.25, ...
        'fineLagStep',0.01,'minCorr',0.8,'sampleTimes2',[], ...
        'qSeparation',0.05,'peakSeparation',1);
    names=fieldnames(defaults);
    for i=1:numel(names)
        if ~isfield(options,names{i}), options.(names{i})=defaults.(names{i}); end
    end
    validateattributes(options.qLimit,{'numeric'},{'scalar','real','finite','>',0,'<',1});
    validateattributes(options.qStep,{'numeric'},{'scalar','real','finite','positive','<=',options.qLimit});
    validateattributes(options.lagStep,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(options.fineLagStep,{'numeric'},{'scalar','real','finite','positive','<=',options.lagStep});
    validateattributes(options.minCorr,{'numeric'},{'scalar','real','finite','>=',-1,'<=',1});
    validateattributes(options.qSeparation,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(options.peakSeparation,{'numeric'},{'scalar','real','finite','positive'});
    first=double(window1(:)); second=double(data2(:));
    L=double(L); Ts=double(Ts); vMin=double(vMin); vMax=double(vMax);
    N=numel(first); assert(N>=5,'В окне нужно минимум пять отсчётов.');
    K=N-1; D=K*Ts;
    vMean=NaN; acceleration=NaN;
    info=struct('accepted',false,'status','no_match','correlation',NaN, ...
        'q',NaN,'lagSamples',NaN,'delay_s',NaN,'duration_s',D, ...
        'timeMidRelative_s',D/2,'onSearchBoundary',false, ...
        'lagPeakMargin',NaN,'qContrast',NaN);
    if isempty(options.sampleTimes2)
        secondIndex=(0:numel(second)-1).';
    else
        validateattributes(options.sampleTimes2,{'numeric'},{'vector','real','finite'});
        secondIndex=double(options.sampleTimes2(:))/Ts;
        assert(numel(secondIndex)==numel(second)&&all(diff(secondIndex)>0)&& ...
            ~isempty(secondIndex)&&abs(secondIndex(1))<1e-7, ...
            'sampleTimes2 должны соответствовать data2, возрастать и начинаться с 0.');
    end
    if numel(second)<2 || secondIndex(end)<K-1e-9
        info.status='need_more_sensor2_data'; return;
    end
    if norm(first-mean(first))<=eps(max(1,norm(first)))
        info.status='flat_first_window'; return;
    end
    lagLo=L/(vMax*Ts); lagHi=L/(vMin*Ts);
    qGrid=unique([-options.qLimit:options.qStep:options.qLimit,0,options.qLimit]);
    rhoGrid=-inf(size(qGrid)); lagGrid=nan(size(qGrid));
    for i=1:numel(qGrid)
        [rhoGrid(i),lagGrid(i)]=searchLag(first,second,secondIndex,K,qGrid(i), ...
            lagLo,lagHi,options.lagStep);
    end
    [best,bestIndex]=max(rhoGrid);
    if ~isfinite(best), return; end
    q=qGrid(bestIndex); lag=lagGrid(bestIndex);
    qFine=parabolicPeak(qGrid,rhoGrid,bestIndex);
    if qFine~=q
        [rhoFine,lagFine]=searchLag(first,second,secondIndex,K,qFine,lagLo,lagHi,options.lagStep);
        if rhoFine>best, q=qFine; lag=lagFine; best=rhoFine; end
    end
    [rhoFine,lagFine]=searchLag(first,second,secondIndex,K,q, ...
        max(lagLo,lag-options.lagStep),min(lagHi,lag+options.lagStep),options.fineLagStep);
    if rhoFine>best, best=rhoFine; lag=lagFine; end
    info.correlation=best; info.q=q; info.lagSamples=lag; info.delay_s=lag*Ts;
    info.onSearchBoundary=abs(abs(q)-options.qLimit)<1e-10 || ...
        abs(lag-lagLo)<1e-10 || abs(lag-lagHi)<1e-10;
    alternatives=rhoGrid(abs(qGrid-q)>=options.qSeparation & isfinite(rhoGrid));
    if ~isempty(alternatives), info.qContrast=max(0,best-max(alternatives)); end
    info.lagPeakMargin=lagPeakMargin(first,second,secondIndex,K,q,lag,best, ...
        lagLo,lagHi,options.lagStep,options.peakSeparation);
    if best<options.minCorr
        info.status='low_correlation'; return;
    end
    vMean=L/(lag*Ts);
    acceleration=2*q*vMean/D;
    info.accepted=true; info.status='ok';
end

function [best,bestLag]=searchLag(first,second,secondIndex,K,q,lo,hi,step)
    grid=(0:K).';
    [z1,ok]=warpedSamples(first,grid,q,K);
    best=-inf; bestLag=NaN;
    if ~ok, return; end
    z1=z1-mean(z1); norm1=norm(z1);
    if norm1<=eps(max(1,norm(first))), return; end
    lags=unique([lo:step:hi,hi]); rho=-inf(size(lags));
    for j=1:numel(lags)
        rho(j)=scoreAtLag(z1,norm1,second,secondIndex,grid,q,K,lags(j));
    end
    [best,index]=max(rho);
    if ~isfinite(best), return; end
    bestLag=lags(index);
    refined=parabolicPeak(lags,rho,index);
    if refined~=bestLag
        value=scoreAtLag(z1,norm1,second,secondIndex,grid,q,K,refined);
        if value>best, best=value; bestLag=refined; end
    end
end

function value=scoreAtLag(z1,norm1,second,secondIndex,grid,q,K,lag)
    value=-inf;
    [z2,ok]=warpedSamples(second,grid+lag,q,K,secondIndex);
    if ~ok, return; end
    z2=z2-mean(z2); norm2=norm(z2);
    if norm2<=eps(max(1,norm(z2))), return; end
    value=max(-1,min(1,dot(z1,z2)/(norm1*norm2)));
end

function [samples,ok]=warpedSamples(signal,newIndex,q,K,sampleIndex)
    r=newIndex/K;
    discriminant=(1-q)^2+4*q*r;
    samples=nan(size(newIndex)); ok=false;
    if any(discriminant<0), return; end
    denominator=(1-q)+sqrt(discriminant);
    if any(denominator<=0), return; end
    sourceIndex=K*(2*r./denominator);
    if nargin<5, sampleIndex=(0:numel(signal)-1).'; end
    first=sampleIndex(1); last=sampleIndex(end);
    if any(sourceIndex<first-1e-7 | sourceIndex>last+1e-7), return; end
    sourceIndex=max(first,min(last,sourceIndex));
    samples=interp1(sampleIndex,signal,sourceIndex,'linear');
    ok=all(isfinite(samples));
end

function margin=lagPeakMargin(first,second,secondIndex,K,q,lag,best,lo,hi,step,separation)
    grid=(0:K).';
    [z1,ok]=warpedSamples(first,grid,q,K);
    margin=NaN; if ~ok, return; end
    z1=z1-mean(z1); norm1=norm(z1);
    lags=unique([lo:step:hi,hi]); values=-inf(size(lags));
    for j=1:numel(lags)
        values(j)=scoreAtLag(z1,norm1,second,secondIndex,grid,q,K,lags(j));
    end
    other=-inf;
    for j=1:numel(lags)
        if ~isfinite(values(j)) || abs(lags(j)-lag)<separation, continue; end
        left=-inf; right=-inf;
        if j>1, left=values(j-1); end
        if j<numel(lags), right=values(j+1); end
        if values(j)>=left && values(j)>=right, other=max(other,values(j)); end
    end
    if isfinite(other), margin=max(0,best-other); else, margin=inf; end
end

function x=parabolicPeak(grid,values,index)
    x=grid(index);
    if index<=1 || index>=numel(grid), return; end
    h=grid(index+1)-grid(index);
    if abs((grid(index)-grid(index-1))-h)>1e-9*max(1,h), return; end
    left=values(index-1); center=values(index); right=values(index+1);
    if any(~isfinite([left,center,right])), return; end
    d=left-2*center+right;
    if d>=-1e-12, return; end
    offset=0.5*(left-right)/d*h;
    x=x+max(-h/2,min(h/2,offset));
end
