%  close all;
%  cas = 1;

 N = 1024;
 t = (0:N-1)/N;
 a = 2;

 %% signal definition
 % different type of monocomponent signals 
 if cas == 1,
  s   = a.*exp(2*pi*1i*(100*t+200*t.^2));
  phi1_ref = 100 + 400*t;
  phi2_ref = 400*ones(size(t));
  sigma_opt = 0.03;
 elseif cas == 2,
  A = 50;
  B = 50;
  C = 50;
  D = 50;
  s =  exp(2*pi*1i*(A*t + B*t.^2 + C*t.^3 + D*t.^4));
  phi1_ref = A + 2*B*t + 3*C*t.^2 + 4*D*t.^3;
  phi2_ref = 2*B + 6*C*t + 12*D*t.^2;
  phi3_ref = 6*C + 24*D*t;
  phi4_ref = 24*D*ones(size(t));
  sigma_opt = 0.03; %estaba en 0.03
 elseif cas == 3,
  a1 = 1+ 5*t.^3 + 7*(1-t).^6;
  f_sin = 10;
  phi1 = 340*t-2.*exp(-2*(t-0.2)).*sin(f_sin*pi.*(t-0.2));
  phi1_ref = 340+4*exp(-2*(t-0.2)).*sin(f_sin*pi.*(t-0.2))-2*f_sin*pi.*exp(-2*(t-0.2)).*cos(f_sin*pi.*(t-0.2)); 
  phi2_ref = -8*exp(-2*(t-0.2)).*sin(f_sin*pi.*(t-0.2))+ 4*f_sin*pi*exp(-2*(t-0.2)).*cos(f_sin*pi.*(t-0.2))...
      +4*f_sin*pi.*exp(-2*(t-0.2)).*cos(f_sin*pi.*(t-0.2))+2*(f_sin*pi)^2*exp(-2*(t-0.2)).*sin(f_sin*pi.*(t-0.2));
  s = a1.*exp(2*pi*1i*(phi1));
  sigma_opt = 0.01; %estaba en 0.01
 else
  return;
 end

%% TF parameters
gamma =10^(-2);
Nfft = N;
[STFT_GT,omega,omega2,phi2sec] = q_bar(s,sigma_opt,Nfft,gamma);
 
 Nrand = 100;
%  Nrand = 3;
SNRs = 0:5:50;
% SNRs = inf;

Nsnr = length(SNRs);
 Rq_hat_data = zeros(Nrand, Nsnr);
%  Rq_bar_data = zeros(Nrand, Nsnr);
 Rq_hfilt_data = zeros(Nrand, Nsnr);
%  Rq_bfilt_data = zeros(Nrand, Nsnr);
%  Rq_h3filt_data = zeros(Nrand, Nsnr);
%  Rq_b3filt_data = zeros(Nrand, Nsnr);
 R_LCR_data = zeros(6, Nrand, Nsnr);
 index = 100;
 Li = (index:N-index);
%  for n_snr = 1:Nsnr
%      snr_in = SNRs(n_snr);
%      snr_in
%      for n_rand=1:Nrand
%          n  = randn(N,1)+1i*randn(N,1);
%          [sn]  = sigmerge(s(:),n,snr_in);
% 
%          [STFT,omega1,omega2,phi2_hat,~] =...
%              q_bar(sn,sigma_opt,Nfft,gamma);
% 
%          [c,e] = exridge(STFT,0,0,10);
%          y_hat = zeros(1, length(t));
%          X_hat = zeros(1, length(t));
%          for k=1:N
%           y_hat(k) = phi2_hat(c(1,k),k);
%           X_hat(k) = omega2(c(1,k),k);
%          end 
% 
%          %% filter estimation of q
% 
%          % FRACTION = .1 for 1% of maximum energy%
%          % FRACTION = .3162 for 10% of maximum energy
%          FRACTION = .3162;
%          frec_corte = sqrt(-lambertw(-FRACTION*exp(-1)))/(sigma_opt*sqrt(pi));
% 
%          %% test symmetrization of y_bar
%          % Results are not satisfactory:
%          % it improves global SNR (with phi2_ref) at borders
%          % but SNR excluding border gets worse
%          
% %          [B, A] = butter(1,2*frec_corte/N,'low');
% %          [B2, A2] = butter(2,2*frec_corte/N,'low');
% %          [B4, A4] = butter(4, frec_corte/N,'low');
% %          y_filt = filter(B4, A4, [y_bar fliplr(y_bar)]);
% %          y_filt = y_filt(1:N);
%          
% %          [B4, A4] = butter(4, 2*frec_corte/(3*N),'low');
% %          y_h3filt = filtfilt(B4, A4, [fliplr(y_hat) y_hat fliplr(y_hat)]);
% %          y_h3filt = y_h3filt(N+1:2*N);
% %          y_b3filt = filtfilt(B4, A4, [fliplr(y_bar) y_bar fliplr(y_bar)]);
% %          y_b3filt = y_b3filt(N+1:2*N);
% 
% %          [B6, A6] = butter(6,2*frec_corte/N,'low');
% 
% %          d = designfilt('lowpassiir','FilterOrder',2, ...
% %                  'PassbandFrequency',frec_corte/N,'PassbandRipple',0.2,'SampleRate',1); 
% 
% %          fvtool(B,A, B2,A2, B4,A4, B6,A6, d);
%          
%         %% butter filter
%         % we do not symmetrize because we exclude borders
%          [B4, A4] = butter(4, 2*frec_corte/N,'low');
%          y_hfilt = filtfilt(B4, A4, y_hat);
%          
%          snr_hat = snr(phi2_ref(Li), phi2_ref(Li) - y_hat(Li));
%          snr_hfilt = snr(phi2_ref(Li), phi2_ref(Li) - y_hfilt(Li));
%          Rq_hat_data(n_rand, n_snr) = snr_hat;
%          Rq_hfilt_data(n_rand, n_snr) = snr_hfilt;
% 
% %          fprintf("snr in %u, snr hat %f, snr hfilt %f, snr bfilt %f\n",...
% %              snr_in, snr_hat, snr_hfilt, snr_bfilt);
% %          
% %          figure;
% %          hold on;
% %          plot(abs(phi2_ref), 'k');
% %          plot(abs(y_hat));
% %          plot(abs(y_filt), 'r--');
% %          hold off;
%          %% apply LCR
% 
% %          t = -0.5:1/N:0.5-1/N;
% %          t=t';
% %          g =  1/sigma_opt*exp(-pi/sigma_opt^2*(t.^2));
% 
%          g = create_gaussian_window(N, Nfft, sigma_opt, 10^(-3));
%          [modes_nf, T_LCR_nf, Lh] = ...
%              LCR(STFT, X_hat, y_hat, g, sigma_opt, N, Nfft, 1);
% 
%          [modes_nf_ref, T_LCR_nf_ref, Lh] = ...
%              LCR(STFT, phi1_ref, y_hat, g, sigma_opt, N, Nfft, 1);
% 
%          [modes_filt, T_LCR_filt, Lh] = ...
%              LCR(STFT, X_hat, y_hfilt, g, sigma_opt, N, Nfft, 1);
% 
%          [modes_filt_ref, T_LCR_filt_ref, Lh] = ...
%              LCR(STFT, phi1_ref, y_hfilt, g, sigma_opt, N, Nfft, 1);
% 
%          [modes_ref, T_LCR_filt_ref, Lh] = ...
%              LCR(STFT, phi1_ref, phi2_ref, g, sigma_opt, N, Nfft, 1);
% 
%          [modes_GT, T_LCR_filt_GT, Lh] = ...
%              LCR(STFT_GT, phi1_ref, phi2_ref, g, sigma_opt, N, Nfft, 1);
% 
%          R_LCR_data(1, n_rand, n_snr) =...
%              snr(s(Li), s(Li) - modes_nf(Li));
% 
%          R_LCR_data(2, n_rand, n_snr) =...
%              snr(s(Li), s(Li) - modes_nf_ref(Li));
%  
%          R_LCR_data(3, n_rand, n_snr) =...
%              snr(s(Li), s(Li) - modes_filt(Li));
% 
%          R_LCR_data(4, n_rand, n_snr) =...
%              snr(s(Li), s(Li) - modes_filt_ref(Li));
% 
%          R_LCR_data(5, n_rand, n_snr) =...
%              snr(s(Li), s(Li) - modes_ref(Li));
% 
%          R_LCR_data(6, n_rand, n_snr) =...
%              snr(s(Li), s(Li) - modes_GT(Li));
%      end
%  end
 
%  save(['data_', prefix, '.mat'], 'Rq_hat_data', 'Rq_hfilt_data', 'R_LCR_data');
load(['data_', prefix, '.mat']);
 
 %% fig cmp q
 plot_Ismall("input SNR", "output SNR");
 hold on;
 plot(SNRs, mean(Rq_hat_data), 'k--', "DisplayName", "$\widehat{q}_{f + n}$");
 plot(SNRs, mean(Rq_hfilt_data), 'g-', "DisplayName", "$F(\widehat{q}_{f + n})$");
 hold off;
 legend_Ismall("southeast");
 savefig(['fig_', prefix, '_cmp_q']);
saveas(gcf,['fig_', prefix, '_cmp_q'],'epsc');
 close all;
 
 %% fig LCR
 disp_ = ["k", "k:", "g-", "g--o", "r--s", "b--*"];
 name_ = ["LCR, $\widehat{q}_{f + n}, \widehat{\omega}^{[2]}_{f + n}$",...
     "LCR, $\widehat{q}_{f + n}, \phi'$",...
     "LCR, $F(\widehat{q}_{f + n}), \widehat{\omega}^{[2]}_{f + n}$",...
     "LCR, $F(\widehat{q}_{f + n}), \phi'$",...
     "LCR, $\phi'', \phi'$",...
     "LCR, ground truth"];
 plot_Ismall("input SNR", "\% of improvement");
 hold on;
 old_ = mean(squeeze(R_LCR_data(1, :, :)));
 for ind=[3, 4, 5]
     mean_ = mean(squeeze(R_LCR_data(ind, :, :)));
     plot(SNRs, 100*(mean_ - old_)./old_, disp_(ind), "DisplayName", name_(ind));
 end
 hold off;
 legend_Ismall("southeast");
 savefig(['fig_', prefix, '_LCR']);
saveas(gcf,['fig_', prefix, '_LCR'],'epsc');
 close all;

 %% fig boxplot
 plot_Ismall("input SNR", "output SNR");
 hold on;
 boxplot(Rq_hat_data, SNRs, 'Color', [0, 0, 0]);
 boxplot(Rq_hfilt_data, SNRs, 'Color', [0, 1, 0]);
 plot(median(Rq_hat_data), 'k--', "DisplayName", "$\widehat{q}_{f + n}$");
 plot(median(Rq_hfilt_data), 'g-', "DisplayName", "$F(\widehat{q})_{f + n}$");
 hold off;
 MIN_v = min(min(min(Rq_hat_data)), min(min(Rq_hfilt_data))) - 2;
 MAX_v = max(max(max(Rq_hat_data)), max(max(Rq_hfilt_data))) + 2;
 ylim([MIN_v, MAX_v]);
 legend_Ismall("southeast");
 savefig(['fig_boxplot_', prefix, '_cmp_q']);
saveas(gcf,['fig_boxplot_', prefix, '_cmp_q'],'epsc');
 close all;