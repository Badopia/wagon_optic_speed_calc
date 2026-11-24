% Просмотр продольного профиля полувагона 12-6995 на высоте датчиков.
scriptFolder = fileparts(mfilename('fullpath'));
profileFile = fullfile(scriptFolder, 'data', 'wagon_12_6995_profile.mat');
if ~isfile(profileFile)
    error('Файл профиля не найден: %s', profileFile);
end

data = load(profileFile, 'profile');
profile = data.profile;

profileFigure = figure('Name', ['Профиль полувагона ' profile.wagon_type], ...
    'NumberTitle', 'off');
ax = axes(profileFigure);
plot(ax, profile.x_m, profile.depth_from_side_m, ...
    'LineWidth', 1.6, 'DisplayName', 'Профиль');
yline(ax, profile.background_depth_from_side_m, '--', ...
    'DisplayName', 'Условный фон');
xlim(ax, [profile.x_m(1), profile.x_m(end)]);
grid(ax, 'on');

title(ax, ['Профиль полувагона ' profile.wagon_type ' на высоте датчиков']);
xlabel(ax, 'Координата от левой концевой балки, м');
ylabel(ax, 'Глубина относительно нуля профиля, м');
legend(ax, 'Location', 'best');
