close all;
% clear all;

x = importdata('62-a_lhl.wav');
fs = x.fs;
x = x.data;
x = resample(x,1,10);
% x = [x; zeros(2*8192-length(x),1)];

Fs = fs/10;
Lx = length(x);

C_Lx = Lx/Fs;

%% filter design

sigma_w = 0.005; % 0.05 para N = 8192, 0.025 para N = 4096 - revisar bien los valores de sigma

% FRACTION = .1 for 1% of maximum energy%
% FRACTION = .3162 for 10% of maximum energy
FRACTION = .3162;
frec_corte = sqrt(-lambertw(-FRACTION*exp(-1)))/(sigma_w*sqrt(pi));

% we do not symmetrize because we exclude borders
[B4, A4] = butter(4, 2*frec_corte/Lx,'low');

%% Time frequency

Nfft = 4096;
Nx = 1:Lx;
t = (Nx - 1)/Fs;

max_f_norm = Nfft/2;

% snr_in = 10;
noise_r = randn(1, Lx);
xn = sigmerge(x, noise_r', snr_in);

K0 = 350;
%% comparison
[STFT_GT,omega,omega2,phi2_hat,phi2_bar] =...
    q_bar_ecg(x,sigma_w,Nfft,max_f_norm);

[c_GT,e] = exridge(STFT_GT(1:K0, :),0,0,2);

index = round(Fs/10);
Li = (index:Lx-index);
y_hat_GT = zeros(1, length(t));
X_hat = zeros(1, length(t));
for k=1:Lx
    y_hat_GT(k) = phi2_hat(c_GT(1,k),k);
    X_hat(k) = omega2(c_GT(1,k),k);
end

y_hfilt_GT = filtfilt(B4, A4, y_hat_GT);

%% realisation over 100 noises

Nrand = 100;
snr_in = 10;

qhat_data = zeros(Nrand, Lx);
qfilt_data = zeros(Nrand, Lx);

for n_rand=1:Nrand
    fprintf("%u/100\n", n_rand);
    noise_r = randn(1, Lx);
    xn = sigmerge(x, noise_r', snr_in);

    [STFT,omega,omega2,phi2_hat,phi2_bar] =...
        q_bar_ecg(xn,sigma_w,Nfft,max_f_norm);

    [c,e] = exridge(STFT(1:K0, :),0,0,2);


    index = Fs/10;
    Li = (index:Lx-index);
    y_hat = zeros(1, length(t));
    for k=1:Lx
        y_hat(k) = phi2_hat(c(1,k),k);
    end

    y_hfilt = filtfilt(B4, A4, y_hat);
    
    qhat_data(n_rand, :) = y_hat;
    qfilt_data(n_rand, :) = y_hfilt;
end

save(['data_voice_', suffix, '.mat'], 'qhat_data', 'qfilt_data');

% This file weights 1.6 Gb
% load("/home/nils/Documents/data_voice.mat");

%% figures
F_vec = ((1:Nfft/2) - 1)*Fs/Nfft;
TFRsc_Ismall(t, F_vec, abs(STFT_GT(1:Nfft/2, :)));
hold on;
plot(t, (c_GT - 1)*Fs/Nfft, 'r', "DisplayName", "Ridge");
hold off;
legend_Ismall("northeast");
saveas(gcf,"fig_voice_TFR",'epsc');
close all;

Y0 = 15000;
plot_Ismall("time", "chirp rate");
hold on;
plot(t(Li), quantile(qhat_data(:, Li), .95), 'b--', "DisplayName", "$Q_{0.95}(\widehat{q}_{f + n})$");
plot(t(Li), y_hat_GT(Li), 'k-', "DisplayName", "$\widehat{q}_{f}$");
plot(t(Li), quantile(qhat_data(:, Li), .05), 'b--', "DisplayName", "$Q_{0.05}(\widehat{q}_{f + n})$");
hold off;
ylim([-Y0, Y0]);
legend_Ismall("northeast");
fname = ['fig_voice_nofilt_', suffix];
savefig(fname);
saveas(gcf,fname,'epsc');
close all;

plot_Ismall("time", "chirp rate");
hold on;
plot(t(Li), quantile(qfilt_data(:, Li), .95), 'b--', "DisplayName", "$Q_{0.95}(F(\widehat{q}_{f + n}))$");
plot(t(Li), y_hfilt_GT(Li), 'k-', "DisplayName", "$F(\widehat{q}_{f})$");
plot(t(Li), quantile(qfilt_data(:, Li), .05), 'b--', "DisplayName", "$Q_{0.05}(F(\widehat{q}_{f + n}))$");
hold off;
ylim([-Y0, Y0]);
legend_Ismall("northeast");
fname = ['fig_voice_filt_', suffix];
savefig(fname);
saveas(gcf,fname,'epsc');
close all;

% figure;
% hold on;
% plot(y_hat_GT, 'k-');
% plot(quantile(qhat_data, .95), 'b--');
% plot(quantile(qhat_data, .05), 'b--');
% plot(y_hfilt_GT, 'g-');
% plot(quantile(qfilt_data, .95), 'g--');
% plot(quantile(qfilt_data, .05), 'g--');
% hold off;

