% Loading data
fromf = open('1_5.mat'); % Assuming it's a file containing some vector data

data = sort(fromf.z);
xlswrite('filename.xlsx', data);

% Histogram of the dataset
figure;
histogram(data);
title('Histogram of the Dataset');

% Empirical Cumulative Distribution Function (ECDF)
figure;
ecdf(data);
title('Empirical Cumulative Distribution Function');

% Minimum and maximum elements, range
min_value = min(data);
max_value = max(data);
range = max_value - min_value;
disp(['Minimum element: ', num2str(min_value)]);
disp(['Maximum element: ', num2str(max_value)]);
disp(['Range: ', num2str(range)]);

% Quantiles
quantiles = quantile(data, [0.25, 0.5, 0.75]);
disp(['25th quantile: ', num2str(quantiles(1))]);
disp(['50th quantile (median): ', num2str(quantiles(2))]);
disp(['75th quantile: ', num2str(quantiles(3))]);

% Estimation of mean and variance
mean_value = mean(data);
variance_value = var(data);
disp(['Mean: ', num2str(mean_value)]);
disp(['Variance: ', num2str(variance_value)]);

% Skewness and kurtosis estimation
skewness_value = skewness(data);
kurtosis_value = kurtosis(data);
disp(['Skewness: ', num2str(skewness_value)]);
disp(['Kurtosis: ', num2str(kurtosis_value)]);

% Significance levels
levels = [0.01, 0.05, 0.1];

for alpha = levels
    edges_normal = 0:0.1:1;
    edges_uniform_01 = 0:0.05:1; % Intervals for U(0;1)
    edges_uniform_symmetric = -1:0.05:1; % Intervals for U(–1;1)
    
    % Chi-square test
    observed_counts_normal = histcounts(data, edges_normal);
    
    % Parameter estimation for normal distribution
    mu_hat_normal = mean(data);
    sigma_hat_normal = std(data);
    
    % Expected frequencies for standard normal distribution
    expected_counts_normal = numel(data) * diff(normcdf(edges_normal, mu_hat_normal, sigma_hat_normal));
    
    % Chi-square test for standard normal distribution
    [h_normal, p_normal, stats_normal] = chi2gof(data, 'Expected', expected_counts_normal, 'Edges', edges_normal, 'Alpha', alpha);
    
    % Chi-square test for standard uniform distribution U(0;1)
    [h_uniform_01, p_uniform_01, stats_uniform_01] = chi2gof(data, 'Expected', diff(edges_uniform_01) * numel(data), 'Edges', edges_uniform_01, 'Alpha', alpha);
    
    % Chi-square test for symmetric uniform distribution U(–1;1)
    [h_uniform_symmetric, p_uniform_symmetric, stats_uniform_symmetric] = chi2gof(data, 'Expected', diff(edges_uniform_symmetric) * numel(data), 'Edges', edges_uniform_symmetric, 'Alpha', alpha);
    
    % Code 2: Kolmogorov-Smirnov tests
    [h_ks_gaussian, p_ks_gaussian, ksstat_gaussian] = kstest(data, 'CDF', makedist('Normal', 'mu', 0, 'sigma', 1), 'Alpha', 0.01);
    % Kolmogorov-Smirnov test for standard uniform distribution U(0;1)
    [h_ks_uniform, p_ks_uniform, ksstat_uniform] = kstest(data, 'CDF', makedist('Uniform', 'lower', 0, 'upper', 1), 'Alpha', 0.01);
    % Kolmogorov-Smirnov test for symmetric uniform distribution U(–1;1)
    [h_ks_symmetric, p_ks_symmetric, ksstat_symmetric] = kstest(data, 'CDF', makedist('Uniform', 'lower', -1, 'upper', 1), 'Alpha', 0.01);
    
    % Code 3: Additional tests
    a = 1;
    A = 1;

    % Data generation
    num_samples = 1000;
    z_values = linspace(-a, a, num_samples);
    w_values = A * cos(pi * z_values / (2 * a));
    CDF_matrix = [z_values', cumsum(w_values') / sum(w_values)];

    % Kolmogorov-Smirnov test
    [h_ks, p_ks, ksstat] = kstest(data, 'CDF', CDF_matrix, 'Alpha', alpha);

    % Generating data sample
    sampleSize = 1000;
    % Probability density estimate for the given distribution
    x = linspace(-a, a, 100);
    pdf_estimate = cos(pi * x / (2 * a));
    
    % Pearson chi-square test
    [h_custom, p_custom, stats_custom] = chi2gof(data, 'Expected', pdf_estimate);

    distributions = {'Uniform U(0,1)','Gaussian N(0,1)', 'Symmetric Uniform U(-1,1)','w(z) = A*cos(pi*z/2a)'};
    h_values = [h_ks_uniform, h_ks_gaussian, h_ks_symmetric, h_ks];
    p_values = [p_ks_uniform, p_ks_gaussian, p_ks_symmetric    ,p_ks];
    ksstat_values = [ksstat_uniform, ksstat_gaussian, ksstat_symmetric, ksstat];
    pirs_h = [h_uniform_01, h_normal, h_uniform_symmetric, h_custom]; 
    pirs_p = [p_uniform_01, p_normal, p_uniform_symmetric, p_custom];
    pirs_s = [stats_normal.chi2stat, stats_uniform_01.chi2stat, stats_normal.chi2stat, stats_uniform_symmetric.chi2stat, stats_custom.chi2stat];
    
    disp(['SIGNIFICANCE LEVEL:', num2str(alpha)])
    
    for i = 1:length(distributions)
        disp(' ')
        disp([char(distributions(i)), ' distribution:']);
        disp(['Kolmogorov-Smirnov test: h = ', num2str(h_values(i)), ', p = ', num2str(p_values(i)), ', statistic = ', num2str(ksstat_values(i))]);
        disp(['Pearson chi-square test: h = ', num2str(pirs_h(i)), ', p = ', num2str(pirs_p(i)), ', statistic = ', num2str(pirs_s(i))]);
        
        % Displaying the result of the Pearson chi-square test
        display_chi_square_test_result(pirs_h(i));
        
        if h_values(i) == 0
            disp('The sample corresponds to the distribution according to Kolmogorov-Smirnov.');
        else
            disp('The sample does not correspond to the distribution according to Kolmogorov-Smirnov.');
        end
    end
end

% Plotting the absolute difference between empirical and model CDF
figure;
subplot(2, 1, 1);
ecdf(data);
hold on;
plot(edges_normal, normcdf(edges_normal, mu_hat_normal, sigma_hat_normal), 'r', 'LineWidth', 2);
legend('Empirical CDF', 'Model CDF (normal distribution)');
title('Empirical and Model CDF (normal distribution)');
hold off;

subplot(2, 1, 2);
modulus_difference = abs(ecdf(data) - normcdf(edges_normal, mu_hat_normal, sigma_hat_normal));
plot(edges_normal, modulus_difference, 'b', 'LineWidth', 2);
title('Modulus of the Difference between Empirical and Model CDF');
xlabel('Data Values');
ylabel('Modulus of the Difference');
legend('Modulus of the Difference');

% Function to display the result of the chi-square test
function display_chi_square_test_result(hypothesis_rejected)
    if hypothesis_rejected
        disp('The sample does not correspond to the distribution according to the Pearson chi-square test');
    else
        disp('The sample corresponds to the distribution according to the Pearson chi-square test');
    end
end
