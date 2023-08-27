% Define constants
GAMMA_THRESHOLD_LOWER = 3.6;
GAMMA_THRESHOLD_UPPER = 5;
GAMMA_COEFFICIENT = 5.75;
GAMMA_EXPONENT = 1.15;
A_GAMMA_FACTOR = 0.287;
SIGMA_LOWER = 0.07;
SIGMA_UPPER = 0.09;
SMALL_THRESHOLD_FACTOR = 0.0005;

% Input wave characteristics
H_s = input('Please insert the significant wave height in meters: '); % Significant wave height (in meters)
T_p = input('Please insert the spectral peak period in seconds: '); % Spectral peak period (in seconds)
storm_duration = input('Please insert the spectral peak period in seconds: '); % Storm duration (in seconds)
deltaT = round(10*T_p)/200; % Sampling deltaT, approximately 20 samples per period

% Calculate parameters
[gamma, A_gamma, omega_p, omega] = calculate_parameters(H_s, T_p, storm_duration,...
    GAMMA_THRESHOLD_LOWER, GAMMA_THRESHOLD_UPPER, GAMMA_COEFFICIENT, GAMMA_EXPONENT, A_GAMMA_FACTOR);
[S_pm, S_J] = calculate_spectra(H_s, omega_p, omega, gamma, A_gamma, SIGMA_LOWER, SIGMA_UPPER);
[omega_Max, N1, EPS, t, eta, DeltaOmega, a, b, omega_W, omega_E] ...
    = calculate_eta(S_J, omega, storm_duration, deltaT,SMALL_THRESHOLD_FACTOR);
plot_results(gamma, omega, S_pm, S_J, omega_E, omega_Max, H_s, T_p, t, eta, storm_duration);
export_data(t, eta);
create_waveDict(N1, EPS, a, b, S_J, omega, DeltaOmega);

%% functions

function [gamma, A_gamma, omega_p, omega] = calculate_parameters(H_s, T_p, storm_duration,...
    GAMMA_THRESHOLD_LOWER, GAMMA_THRESHOLD_UPPER, GAMMA_COEFFICIENT, GAMMA_EXPONENT, A_GAMMA_FACTOR)
    X1 = T_p / (H_s^0.5);
    if X1 <= GAMMA_THRESHOLD_LOWER
        gamma = 5;
    elseif X1 >= GAMMA_THRESHOLD_UPPER
        gamma = 1;
    else
        gamma = exp(GAMMA_COEFFICIENT - GAMMA_EXPONENT * X1);
    end
    A_gamma = 1 - A_GAMMA_FACTOR * log(gamma); % Normalizing factor
    omega_p = 2*pi / T_p; % Peak angular frequency
    deltaOmega_minimum = 2*pi / storm_duration;
    omega = 2*pi*1e-4 : deltaOmega_minimum : 20*pi; % Angular spectral frequency distribution [long period waves - capillary waves]
end

function [S_pm, S_J] = calculate_spectra(H_s, omega_p, omega, gamma, A_gamma, SIGMA_LOWER, SIGMA_UPPER)
    S_pm = (5/16) * (H_s^2) * (omega_p^4) .* ((omega).^-5) .* exp((-1.25) .* (((1/omega_p) .* omega).^(-4))); % Pierson-Moskowitz spectrum
    Sigma = zeros(size(omega));
    Sigma(omega <= omega_p) = SIGMA_LOWER;
    Sigma(omega > omega_p) = SIGMA_UPPER;
    S_J = S_pm .* A_gamma .* gamma.^(exp(-0.5.*((omega-omega_p)./(Sigma.*omega_p)).^2)); % JONSWAP spectrum
end

function [omega_Max, N1, EPS, t, eta, DeltaOmega, a, b, omega_W, omega_E] ...
    = calculate_eta(S_J, omega, storm_duration, deltaT,SMALL_THRESHOLD_FACTOR)
    SMALL_THRESHOLD = SMALL_THRESHOLD_FACTOR * max(S_J); % Ignore small values
    S_J(S_J < SMALL_THRESHOLD) = 0; % Set values below the threshold to zero
    % Find indices where spectrum transitions across the threshold
    indices = find(diff(S_J <= SMALL_THRESHOLD));
    % Extract the relevant indices
    a = indices(1);
    b = indices(2);
    omega_Max = omega(S_J == max(S_J));
    omega_W = omega(a);
    omega_E = omega(b);
    N1 = b - a + 1; % Renumbering the spectrum
    DeltaOmega = (omega_E - omega_W) / (2 * N1);
    EPS = 2 * pi * rand(1, N1); % Random phases uniformly distributed between 0 and 2*pi
    t = 0 : deltaT : storm_duration; % Time series
    eta = NaN(1, length(t));    % the surface elevation 
    for i = 1 : length(t)
        sum = 0;
        for j = 1 : N1
            sum = sum + sqrt(2) * (2 * S_J(j) * DeltaOmega)^0.5 * cos(omega(j) * t(i) + EPS(j));
        end
        eta(i) = sum;
    end
end

function plot_results(gamma, omega, S_pm, S_J, omega_E, omega_Max, H_s, T_p, t, eta, storm_duration)
    f1 = figure;
    plot(omega, S_pm, '--', 'LineWidth', 2, 'Color', [1 0 0]);
    hold on
    plot(omega, S_J, 'LineWidth', 2, 'Color', [0 0 1]);
    title('Pierson-Moskowitz and JONSWAP spectra');
    ylabel('S(\omega)', 'FontSize', 16);
    xlabel('\omega', 'FontSize', 16);
    text(0, max(S_pm), 'S_{PM} \rightarrow', 'Color', [1 0 0], 'FontSize', 18);
    text(find(S_J == max(S_J)) * 0.008, max(S_J), ...
        ['\leftarrow S_{JONSWAP} (\gamma = ', num2str(round(gamma, 2)), ')'], 'Color', [0 0 1], 'FontSize', 18);
    text(0.7 * omega_E, 0.4 * max(S_J), ['\omega_{max} = ', ...
        num2str(round(omega_Max, 3)), ' rad/s', newline, 'H_s = ', ...
        num2str(H_s), ' m', newline, 'T_p = ', ...
        num2str(T_p), ' s'], 'Color', 'black', 'FontSize', 10);
    axis([0 omega_E 0 1.1 * max(S_J)]);
    print(f1, '-r300', '-djpeg', 'JONSWAP_spectrum');
    close(f1);

    f2 = figure;
    plot(t, eta);
    title('A possible realization of the storm');
    ylabel('\eta (m)', 'FontSize', 16);
    xlabel('t (s)', 'FontSize', 16);
    axis([0 storm_duration 1.1 * min(eta) 1.1 * max(eta)]);
    print(f2, '-r300', '-djpeg', 'possible_realization_storm');
    close(f2);

    f3 = figure;
    hist(eta, 100);
    [counts, centers] = hist(eta, 100);
    title('Surface elevation histogram');
    ylabel('Frequency', 'FontSize', 16);
    xlabel('\eta (m)', 'FontSize', 16);
    [meanValue, medianValue, modeValue, skewnessValue, kurtosisValue, varianceValue] =...
    realization_statistics(eta, counts, centers);
    text(min(centers), max(counts), ['mean = ' num2str(meanValue)], 'FontSize', 10);
    text(min(centers), 0.95 * max(counts), ['median = ' num2str(medianValue)], 'FontSize', 10);
    text(min(centers), 0.9 * max(counts), ['mode = ' num2str(modeValue)], 'FontSize', 10);
    text(min(centers), 0.85 * max(counts), ['variance = ' num2str(varianceValue)], 'FontSize', 10);
    text(min(centers), 0.8 * max(counts), ['skewness = ' num2str(skewnessValue)], 'FontSize', 10);
    text(min(centers), 0.75 * max(counts), ['kurtosis = ' num2str(kurtosisValue)], 'FontSize', 10);
    axis([1.05 * min(centers) 1.05 * max(centers) 0 1.1 * max(counts)]);
    print(f3, '-r300', '-djpeg', 'Histogram');
    close(f3);
end

function [meanValue, medianValue, modeValue, skewnessValue, kurtosisValue, varianceValue] =...
    realization_statistics(eta, counts, centers)
    % Statistics parameters
    meanValue = mean(eta);
    medianValue = median(eta);
    modeValue = centers(counts == max(counts));
    skewnessValue = skewness(eta);
    kurtosisValue = kurtosis(eta);
    varianceValue = var(eta);
end

function export_data(t, eta)
    csvwrite('t_eta.csv', [t.', eta.']);
end

function create_waveDict(N1, EPS, a, b, S_J, omega, DeltaOmega)
    % OpenFOAM waveDict
    omega_final = omega(a:b);
    S_final = S_J(a:b);
    period_final = 2 * pi ./ omega_final;
    height_final = 2 * sqrt(2 * S_final .* DeltaOmega);
    % Define the file content as a string
    fileContent_1 = ['/*---------------------------------------------------------------------------*\' newline ...
        '| =========                 |                                                 |' newline ...
        '| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |' newline ...
        '|  \\    /   O peration     | Version:  1.3                                   |' newline ...
        '|   \\  /    A nd           | Web:      http://www.openfoam.org               |' newline ...
        '|    \\/     M anipulation  |                                                 |' newline ...
        '\*---------------------------------------------------------------------------*/' newline ...
        'FoamFile' newline ...
        '{' newline ...
        '    version     2.0;' newline ...
        '    format      ascii;' newline ...
        '    class       dictionary;' newline ...
        '    location    "constant";' newline ...
        '    object      waveDict;' newline ...
        '}' newline ...
        '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //' newline newline ...
        'waveType        irregular;' newline newline ...
        'genAbs          1;' newline newline ...
        'absDir          0.0;' newline newline ...
        'nPaddles        1;' newline newline ...
        'secondOrder    0; // 0/1, true/false, on/off' newline newline ...
        'waveHeights' newline];
        fileContent_2 = [newline '(' newline];
        fileContent_3 = [');' newline newline ...
        'wavePeriods' newline];
        fileContent_4 = [newline '(' newline];
        fileContent_5 = [');' newline newline ...
        'wavePhases' newline];
        fileContent_6 = [newline '(' newline];
        fileContent_7 = [');' newline newline ...
        'waveDirs' newline];
        fileContent_8 = [newline '{' newline ...
        '0' newline ...
        '};' newline newline ...
        '// ************************************************************************* //'];
    fileContent = [fileContent_1 num2str(N1) fileContent_2 num2str(height_final) fileContent_3 num2str(N1) fileContent_4...
        num2str(period_final) fileContent_5 num2str(N1) fileContent_6 num2str(EPS) fileContent_7 num2str(N1) fileContent_8 ];
    % Specify the file path and name
    base_directory = pwd;
    filePath = base_directory;
    filename = 'waveDict';
    % Create the file
    fid = fopen(fullfile(filePath, filename), 'w');
    fprintf(fid, '%s', fileContent);
    fclose(fid);
    disp('File created successfully.');
end