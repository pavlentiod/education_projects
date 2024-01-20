% Загрузка данных
fromf = open('1_5.mat'); % Считываем файл

data = sort(fromf.z);
xlswrite('filename.xlsx',data);
% Гистограмма выборки
figure;
histogram(data);
title('Гистограмма выборки');

% Эмпирическая функция распределения
figure;
ecdf(data);
title('Эмпирическая функция распределения');

% Минимальный и максимальный элементы, размах
min_value = min(data);
max_value = max(data);
range = max_value - min_value;
disp(['Минимальный элемент: ', num2str(min_value)]);
disp(['Максимальный элемент: ', num2str(max_value)]);
disp(['Размах: ', num2str(range)]);

% Квантили
quantiles = quantile(data, [0.25, 0.5, 0.75]);
disp(['25% квантиль: ', num2str(quantiles(1))]);
disp(['50% квантиль (медиана): ', num2str(quantiles(2))]);
disp(['75% квантиль: ', num2str(quantiles(3))]);

% Оценка математического ожидания и дисперсии
mean_value = mean(data);
variance_value = var(data);
disp(['Математическое ожидание: ', num2str(mean_value)]);
disp(['Дисперсия: ', num2str(variance_value)]);

% Оценка коэффициентов асимметрии и эксцесса
skewness_value = skewness(data);
kurtosis_value = kurtosis(data);
disp(['Коэффициент асимметрии: ', num2str(skewness_value)]);
disp(['Эксцесс: ', num2str(kurtosis_value)]);

% Уровни значимости
levels = [0.01, 0.05, 0.1];

for alpha = levels
    edges_normal = 0:0.1:1; 
    edges_uniform_01 = 0:0.05:1; % Интервалы для U(0;1)
    edges_uniform_symmetric = -1:0.05:1; % Интервалы для U(–1;1)
    
    % Тест хи-квадрат
    observed_counts_normal = histcounts(data, edges_normal);
    
    % Оценка параметров нормального распределения по выборке
    mu_hat_normal = mean(data);
    sigma_hat_normal = std(data);
    
    % Ожидаемые частоты для стандартного нормального распределения
    expected_counts_normal = numel(data) * diff(normcdf(edges_normal, mu_hat_normal, sigma_hat_normal));
    
    % Статистика и тест хи-квадрат для стандартного нормального распределения
    [h_normal, p_normal, stats_normal] = chi2gof(data, 'Expected', expected_counts_normal, 'Edges', edges_normal, 'Alpha', alpha);
    
    % Тест хи-квадрат для стандартного равномерного распределения U(0;1)
    [h_uniform_01, p_uniform_01, stats_uniform_01] = chi2gof(data, 'Expected', diff(edges_uniform_01) * numel(data), 'Edges', edges_uniform_01, 'Alpha', alpha);
    
    % Тест хи-квадрат для симметричного равномерного распределения U(–1;1)
    [h_uniform_symmetric, p_uniform_symmetric, stats_uniform_symmetric] = chi2gof(data, 'Expected', diff(edges_uniform_symmetric) * numel(data), 'Edges', edges_uniform_symmetric, 'Alpha', alpha);
    
    % Code 2: Тесты Колмогорова-Смирнова
    [h_ks_gaussian, p_ks_gaussian, ksstat_gaussian] = kstest(data, 'CDF', makedist('Normal', 'mu', 0, 'sigma', 1), 'Alpha', 0.01);
    % Тест на соответствие стандартному равномерному распределению U(0;1)
    [h_ks_uniform, p_ks_uniform, ksstat_uniform] = kstest(data, 'CDF', makedist('Uniform', 'lower', 0, 'upper', 1), 'Alpha', 0.01);
    % Тест на соответствие симметричному равномерному распределению U(–1;1)
    [h_ks_symmetric, p_ks_symmetric, ksstat_symmetric] = kstest(data, 'CDF', makedist('Uniform', 'lower', -1, 'upper', 1), 'Alpha', 0.01);
    
    % Code 3: Дополнительные тесты
    a = 1;
    A = 1;

    % Генерация данных
    num_samples = 1000;
    z_values = linspace(-a, a, num_samples);
    w_values = A * cos(pi * z_values / (2 * a));
    CDF_matrix = [z_values', cumsum(w_values') / sum(w_values)];

    % Тест Колмогорова-Смирнова
    [h_ks, p_ks, ksstat] = kstest(data, 'CDF', CDF_matrix, 'Alpha', alpha);

    % Генерация выборки данных
    sampleSize = 1000;
    % Оценка плотности вероятности для заданного распределения
    x = linspace(-a, a, 100);
    pdf_estimate = cos(pi * x / (2 * a));
    
    % Тест Пирсона
    [h_custom, p_custom, stats_custom] = chi2gof(data, 'Expected', pdf_estimate);

    distributions = {'Равномерное U(0,1)','Гауссовское N(0,1)', 'Симметричное равномерное U(-1,1)','w(z) = A*cos(pi*z/2a)'};
    h_values = [h_ks_uniform, h_ks_gaussian, h_ks_symmetric, h_ks];
    p_values = [p_ks_uniform, p_ks_gaussian, p_ks_symmetric    ,p_ks];
    ksstat_values = [ksstat_uniform, ksstat_gaussian, ksstat_symmetric, ksstat];
    pirs_h = [h_uniform_01, h_normal, h_uniform_symmetric, h_custom]; 
    pirs_p = [p_uniform_01, p_normal, p_uniform_symmetric, p_custom];
    pirs_s = [stats_normal.chi2stat, stats_uniform_01.chi2stat, stats_normal.chi2stat, stats_uniform_symmetric.chi2stat, stats_custom.chi2stat];
    
    disp(['УРОВЕНЬ ЗНАЧИМОСТИ:', num2str(alpha)])
    
    for i = 1:length(distributions)
        disp(' ')
        disp([char(distributions(i)), ' распределение:']);
        disp(['Тест Колмогорова-Смирнова: h = ', num2str(h_values(i)), ', p = ', num2str(p_values(i)), ', статистика = ', num2str(ksstat_values(i))]);
        disp(['Тест Пирсона: h = ', num2str(pirs_h(i)), ', p = ', num2str(pirs_p(i)), ', статистика = ', num2str(pirs_s(i))]);
        
        % Вывод результата теста Пирсона
        display_chi_square_test_result(pirs_h(i));
        
        if h_values(i) == 0
            disp('Выборка соответствует распределению по Колмагорову-Смирнову.');
        else
            disp('Выборка не соответствует распределению по Колмагорову-Смирнову.');
        end
    end
end

% График модуля разницы между эмпирической и модельной ФР
figure;
subplot(2, 1, 1);
ecdf(data);
hold on;
plot(edges_normal, normcdf(edges_normal, mu_hat_normal, sigma_hat_normal), 'r', 'LineWidth', 2);
legend('Эмпирическая ФР', 'Модельная ФР (нормальное распределение)');
title('Эмпирическая и Модельная ФР (нормальное распределение)');
hold off;

subplot(2, 1, 2);
modulus_difference = abs(ecdf(data) - normcdf(edges_normal, mu_hat_normal, sigma_hat_normal));
plot(edges_normal, modulus_difference, 'b', 'LineWidth', 2);
title('Модуль Разницы между Эмпирической и Модельной ФР');
xlabel('Значения данных');
ylabel('Модуль Разницы');
legend('Модуль Разницы');


% Функция для отображения результата теста
function display_chi_square_test_result(hypothesis_rejected)
    if hypothesis_rejected
        disp('Выборка не соответствует распределению по критерию Пирсона');
    else
        disp('Выборка соответствует распределению по критерию Пирсона');
    end
end

