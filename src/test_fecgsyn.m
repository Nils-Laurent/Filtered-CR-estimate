addpath(genpath('./fecgsyn-master'));
% close all;
clear workspace;

Fs = 1000;
C_Lx = 15;
Lx = C_Lx*Fs;

param.fs = Fs; % sampling frequency [Hz]
param.n = Lx;

param.mvcg = 5; % choose mother VCG (if empty then the simulator randomly
                %    choose one within the set of available VCGs)
param.posdev = 0; % slight deviation from default hearts and electrodes
                  %   positions (0: hard coded values, 1: random deviation
                  %   and phase initialisation)

param.mres = 0.25; % mother respiration frequency
param.mhr = 80;
param.macc = 20;
param.mtypeacc = 'tanh';
param.flo = 0;
param.fhi = 0;
param.hrstd = 0;
param.lfhfr = 0;


% strhrv.hr = param.mhr;
% strhrv.lfhfr = 0.6;
% strhrv.hrstd = 2;%use to be 2
% strhrv.flo = param.mres;
% strhrv.fhi = 0.25; %use to be 0.25
% strhrv.acc = param.macc;
% strhrv.typeacc = param.mtypeacc;
% strhrv.accmean = param.maccmean;
% strhrv.accstd = param.maccstd;
% [theta_m,w_m] = generate_hrv(strhrv,param.n,param.fs,theta0_m);



% debug = 0;
% FM_syn_ecg = run_ecg_generator(param,debug);
% % out=clean_compress(out); % compression du signal
% 
% sr_ecg_A = FM_syn_ecg.mixture(1, :);
% sr_ecg_M = FM_syn_ecg.mecg(1, :);
% sr_ecg_F = FM_syn_ecg.fecg{1, 1}(1, :);

load("sim_fecgsyn_1.mat");

% return;
%% Time frequency

% Fs=1000;
Nx = 1:Lx;
t = (Nx - 1)/Fs;
s_ecg_m = sr_ecg_M(Nx);
m_qrs = FM_syn_ecg.mqrs;

prec_bpm = 0.25; % frequency bin per bpm
nBpm = Fs/2*60;
Nfft = ceil(2*nBpm*prec_bpm);

%% filter design

% FRACTION = .1 for 1% of maximum energy%
% FRACTION = .3162 for 10% of maximum energy
FRACTION = .3162;
frec_corte = sqrt(-lambertw(-FRACTION*exp(-1)))/(sigma_w*sqrt(pi));

%% comparison

[STFT_GT,omega,omega2,phi2_hat,phi2_bar] =...
    q_bar_ecg(s_ecg_m,sigma_w,Nfft,max_f_norm);

[c,e] = exridge(STFT(1:120, :),0,0,2);


index = Fs/10;
Li = (index:Lx-index);
y_hat_GT = zeros(1, length(t));
X_hat = zeros(1, length(t));
for k=1:Lx
    y_hat_GT(k) = phi2_hat(c(1,k),k);
    X_hat(k) = omega2(c(1,k),k);
end

%% butter filter
% we do not symmetrize because we exclude borders
[B4, A4] = butter(4, 2*frec_corte/Fs,'low');
y_hfilt_GT = filtfilt(B4, A4, y_hat_GT);

%% realisation over 100 noises

Nrand = 100;
snr_in = 10;

qfilt_data = zeros(Nrand, Lx);
qhat_data = zeros(Nrand, Lx);

for n_rand=1:Nrand
    noise_r = randn(1, Lx);
    x = sigmerge(s_ecg_m', noise_r', snr_in);

    max_f = 10;
    max_f_norm = max_f*C_Lx;
    % max_f_norm

    [g, Lh, sigma_w] = create_gaussian_window_len(Lx, Nfft, 2.5*Fs + 1, 10^(-3));
    length(g)

    [STFT,omega,omega2,phi2_hat,phi2_bar] =...
        q_bar_ecg(x,sigma_w,Nfft,max_f_norm);

    [c,e] = exridge(STFT(1:120, :),0,0,2);


    index = Fs/10;
    Li = (index:Lx-index);
    y_hat = zeros(1, length(t));
    X_hat = zeros(1, length(t));
    for k=1:Lx
        y_hat(k) = phi2_hat(c(1,k),k);
        X_hat(k) = omega2(c(1,k),k);
    end

    %% filter estimation of q

    % FRACTION = .1 for 1% of maximum energy%
    % FRACTION = .3162 for 10% of maximum energy
    FRACTION = .3162;
    frec_corte = sqrt(-lambertw(-FRACTION*exp(-1)))/(sigma_w*sqrt(pi));

    %% butter filter
    % we do not symmetrize because we exclude borders
    [B4, A4] = butter(4, 2*frec_corte/Fs,'low');
    y_hfilt = filtfilt(B4, A4, y_hat);
    
    qhat_data(n_rand, :) = y_hat;
    qfilt_data(n_rand, :) = y_hfilt;
end

save("data_fecgsyn.mat");

%% figures

figure;
hold on;
plot(y_hat_GT, 'k-');
plot(quantile(qhat_data, .95), 'b--');
plot(quantile(qhat_data, .05), 'b--');
hold off;

figure;
hold on;
plot(y_hfilt_GT, 'k-');
plot(quantile(qfilt_data, .95), 'b--');
plot(quantile(qfilt_data, .05), 'b--');
hold off;

% TFRsc(abs(STFT));
% hold on;
% plot(c, 'r');
% hold off;
% 
% figure;
% plot(t, X_hat);
% 
% figure;
% hold on
% % plot(t, X_hat);
% plot(t, y_hat_GT, 'g:');
% plot(t, y_hfilt, '--');




