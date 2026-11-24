clear

mainRoot = fileparts(mfilename('fullpath'));
addpath(mainRoot); % solve.m должен быть доступен при последующем запуске Simulink.
run(fullfile(mainRoot, 'parametrs.m'));

%% Загрузка одного профиля вагона
clear profile profileX profileDepth profileMin profileMax
assert(isfile(profileFile), 'Не найден файл профиля: %s', profileFile);
profileData = load(profileFile, 'profile');
assert(isfield(profileData, 'profile') && isstruct(profileData.profile) && ...
    isscalar(profileData.profile) && isfield(profileData.profile, 'x_m') && ...
    isfield(profileData.profile, 'depth_from_side_m'), ...
    'MAT-файл должен содержать profile.x_m и profile.depth_from_side_m.');
profile = profileData.profile;
validateattributes(profile.x_m, {'numeric'}, {'vector','real','finite','nonempty'});
validateattributes(profile.depth_from_side_m, {'numeric'}, {'vector','real','finite','nonempty'});
profileX = double(profile.x_m(:));
profileDepth = double(profile.depth_from_side_m(:));
assert(numel(profileX) >= 2 && numel(profileX) == numel(profileDepth) && ...
    all(diff(profileX) > 0), ...
    'Нужны одинаковые длины векторов и строго возрастающие координаты профиля.');
profileMin = profileX(1);
profileMax = profileX(end);
clear profileData
fprintf('Загружен профиль: %s\n', profileFile);
fprintf('Точек: %d. Координата: %.3f .. %.3f м.\n', ...
    numel(profileX), profileMin, profileMax);
fprintf('Теперь запустите calculation.slx кнопкой Run.\n');
