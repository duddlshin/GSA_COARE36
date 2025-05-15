%% 12.870 final project: global sensitivity analysis for coefficients in COARE 3.6
% Author: Ethan YoungIn Shin
% Last updated: 04/30/2025

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% README %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code runs a global sensitivity analysis using a variance-based
% method described Saltelli 2010 or Jansen 1999. The only values that need
% to be configured are the sample size N found in line 15 and the model
% output for which the sensivity to coefficients is being assessed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

N = 10000;   % Sample size = [1000, 10000]
qoi_option = 1;   % Optoins for model output: option 1: fluxes, option 2: transfer coefficients


% coefficients: [a1, a2, gamma]
a1=0.0017;
a2=-0.0050;
gamma=0.11;
coefficients = [a1, a2, gamma];
[output, qoi] = run_forward_model(coefficients, qoi_option);

y=load('VP_test_data_1.txt');
jdy=y(:,1)-y(1,1)+datenum(2000,1,8)-datenum(1999,12,31);%julian day
u=y(:,5);%COLUMN E: Wind speed (m/sec, 2m)
usr=output(:,1);%%   usr = friction velocity (m/s)
U10n=output(:,31);%   U10n = 10-m neutral wind (m/s)

yy=load('VP_test_results_1.txt');
%COLUMN A: DATE & TIME
Rsday=yy(:,2);%COLUMN B: DAILY INSOLATION (8 same values)
Rl=yy(:,3);%COLUMN C: LONG-WAVE RADIATION
Hs=-yy(:,4);%COLUMN D: SENSIBLE HEAT
Hl=-yy(:,5);%COLUMN E: LATENT HEAT
Tau=yy(:,6);%COLUMN F: WIND STRESS

% Plot original data from exercise;
figure; hold on;
plot(jdy,output(:,2))
plot(jdy,Tau,'x');
xlabel('Julian Day');
ylabel('Stress ($N/m^2$)');
legend('Prediction','Data');
plot_settings;

figure; hold on;
plot(jdy,output(:,3));
plot(jdy,Hs,'x');
xlabel('Julian Day');
ylabel('Sensible Heat ($W/m^2$)');
legend('Prediction','Data');
plot_settings;

figure; hold on;
plot(jdy,output(:,4));
plot(jdy,Hl,'x');
xlabel('Julian Day');
ylabel('Latent heat ($W/m^2$)');
legend('Prediction','Data');
plot_settings;

figure; hold on;
plot(jdy,u);
plot(jdy,U10n,'x');
xlabel('Julian Day');
ylabel('Wind Speed ($m/s$)');
legend('u input height','U10n');
plot_settings;


%% Make matrices 

% bounds for the coefficient distributions
a1_bounds = [a1*0.5, a1*1.5];
a2_bounds = [a2*0.5, a2*1.5];
gamma_bounds = [gamma*0.5, gamma*1.5];

% matrix A
rng(42); % random seed
A_a1 = a1_bounds(1) + (a1_bounds(2) - a1_bounds(1)) * rand(N, 1); 
A_a2 = a2_bounds(1) + (a2_bounds(2) - a2_bounds(1)) * rand(N, 1); 
A_gamma = gamma_bounds(1) + (gamma_bounds(2) - gamma_bounds(1)) * rand(N, 1); 
A = [A_a1, A_a2, A_gamma];

% matrix B
rng(43); % random seed
B_a1 = a1_bounds(1) + (a1_bounds(2) - a1_bounds(1)) * rand(N, 1); 
B_a2 = a2_bounds(1) + (a2_bounds(2) - a2_bounds(1)) * rand(N, 1); 
B_gamma = gamma_bounds(1) + (gamma_bounds(2) - gamma_bounds(1)) * rand(N, 1); 
B = [B_a1, B_a2, B_gamma];

% matrix A_B_{coefficient}  (matrix B with only {coefficient} values from matrix A)
A_B_a1 = [A_a1, B_a2, B_gamma];
A_B_a2 = [B_a1, A_a2, B_gamma];
A_B_gamma = [B_a1, B_a2, A_gamma];



%% Generate f(A), f(B) and f(A_B_{coefficient})

fA = zeros(N,3);
fB = zeros(N,3);
fA_B_a1 = zeros(N,3);
fA_B_a2 = zeros(N,3);
fA_B_gamma = zeros(N,3);

% loop
for i = 1:N
    [~,fA(i,:)] = run_forward_model(A(i,:), qoi_option);
    [~,fB(i,:)] = run_forward_model(B(i,:), qoi_option);
    [~,fA_B_a1(i,:)] = run_forward_model(A_B_a1(i,:), qoi_option);
    [~,fA_B_a2(i,:)] = run_forward_model(A_B_a2(i,:), qoi_option);
    [~,fA_B_gamma(i,:)] = run_forward_model(A_B_gamma(i,:), qoi_option);
end


%% Calculate sensitivity indices (Jansen 1999 version)

% total variance
var_AB = (1/(2*N)) * sum((fA - fB).^2);

% first-order sensitivity indices
S_a1 = abs(var_AB - (1/(2*N)) * sum((fA - fA_B_a1).^2)) ./ var_AB;
S_a2 = abs(var_AB - (1/(2*N)) * sum((fA - fA_B_a2).^2)) ./ var_AB;
S_gamma = abs(var_AB - (1/(2*N)) * sum((fA - fA_B_gamma).^2)) ./ var_AB;

% total sensitivity indices
S_T_a1 = (1/(2*N)) * sum((fB - fA_B_a1).^2) ./ var_AB;
S_T_a2 = (1/(2*N)) * sum((fB - fA_B_a2).^2) ./ var_AB;
S_T_gamma = (1/(2*N)) * sum((fB - fA_B_gamma).^2) ./ var_AB;



%% Plot sensitivity indices for fluxes or transfer coefficients

% grouping sensitivities indices
S_values = [S_a1; S_a2; S_gamma];
S_T_values = [S_T_a1; S_T_a2; S_T_gamma];

% coefficients labels
% Note: a1 in the code is m, a2 in the code is b, gamma is gamma
labels = {'$m$', '$b$', '$\gamma$'};


if qoi_option == 1
    figure("Position",[100,100,1100,400]); hold on;
    
    % first-order sensitivity indices
    subplot(1, 2, 1);
    bar(S_values);
    set(gca, 'XTickLabel', labels, 'XTick', 1:numel(labels));
    title('First Order Sensitivity Indices');
    xlabel('Coefficients');
    ylabel('$S_i$', 'Interpreter', 'latex');
    yticks([0, 0.5, 1.0]);
    yticklabels({'0.0', '0.5', '1.0'});
    ylim([0,1.1]);
    legend("Momentum", "Sensible Heat", "Latent Heat", "Location","northeast");
    plot_settings;
    
    % total sensitivity indices
    subplot(1, 2, 2);
    bar(S_T_values);
    set(gca, 'XTickLabel', labels, 'XTick', 1:numel(labels));
    title('Total Sensitivity Indices');
    xlabel('Coefficients');
    ylabel('$S_{T_i}$', 'Interpreter', 'latex');
    yticks([0, 0.5, 1.0]);
    yticklabels({'0.0', '0.5', '1.0'});
    ylim([0,1.1]);
    plot_settings;

elseif qoi_option == 2
    figure("Position",[100,100,1100,400]); hold on;
    
    % first-order sensitivity indices
    subplot(1, 2, 1);
    bar(S_values);
    set(gca, 'XTickLabel', labels, 'XTick', 1:numel(labels));
    title('First Order Sensitivity Indices');
    xlabel('Coefficients');
    ylabel('$S_i$', 'Interpreter', 'latex');
    yticks([0, 0.5, 1.0]);
    yticklabels({'0.0', '0.5', '1.0'});
    ylim([0,1.1]);
    legend("$C_D$", "$C_H$", "$C_E$", "Location","northeast");
    plot_settings;
    
    % total sensitivity indices
    subplot(1, 2, 2);
    bar(S_T_values);
    set(gca, 'XTickLabel', labels, 'XTick', 1:numel(labels));
    title('Total Sensitivity Indices');
    xlabel('Coefficients');
    ylabel('$S_{T_i}$', 'Interpreter', 'latex');
    yticks([0, 0.5, 1.0]);
    yticklabels({'0.0', '0.5', '1.0'});
    ylim([0,1.1]);
    plot_settings;
end