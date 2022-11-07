 close all;
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
[STFT_GT,omega,omega2,phi2_GT] = q_bar(s,sigma_opt,Nfft,gamma);
 
 snr_in = 10;
 index = 100;
 Li = (index:N-index);

 n  = randn(N,1)+1i*randn(N,1);
 [sn]  = sigmerge(s(:),n,snr_in);

 [STFT,omega1,omega2,phi2_hat,~] =...
     q_bar(sn,sigma_opt,Nfft,gamma);

 [c,e] = exridge(STFT,0,0,10);
 y_hat = zeros(1, length(t));
 X_hat = zeros(1, length(t));
 for k=1:N
  y_hat(k) = phi2_hat(c(1,k),k);
  X_hat(k) = omega2(c(1,k),k);
 end 

 %% filter estimation of q

 % FRACTION = .1 for 1% of maximum energy
 % FRACTION = .3162 for 10% of maximum energy
 FRACTION = .3162;
 frec_corte = sqrt(-lambertw(-FRACTION*exp(-1)))/(sigma_opt*sqrt(pi));

%% butter filter
% we do not symmetrize because we exclude borders
 [B4, A4] = butter(4, 2*frec_corte/N,'low');
 y_hfilt = filtfilt(B4, A4, y_hat);
 
 fc2 = frec_corte;
 fft_y = fft(y_hat);
 fft_LP = fft(y_hfilt);
 fft_GT = fft(phi2_ref);
 f_vec = 0:(length(phi2_ref) - 1);
 
 plot_Ismall("frequency", "magnitude");
 hold on;
 plot(f_vec, abs(fft_GT), "DisplayName", "ground truth");
 Z = ylim;
 plot([fc2, fc2], [Z(1), Z(2)], '--',...
     "DisplayName", "cut-off frequency");
 hold off;
 
 xlim([0, 2*fc2]);
 legend_Ismall("northeast");
 savefig(['fig_', prefix, '_cutoff']);
saveas(gcf,['fig_', prefix, '_cutoff'],'epsc');
close all;
 
 return;
 
 plot_Ismall("time", "estimation of $\phi''$");
 hold on;
 plot(t(Li), phi2_ref(Li), 'r-.', "DisplayName", "$\phi''$");
 plot(t(Li), y_hat(Li), 'k--', "DisplayName", "$\widehat{q}_{f+n}$");
 plot(t(Li), y_hfilt(Li), 'g-', "DisplayName", "$F(\widehat{q}_{f+n})$");
 hold off;
 legend_Ismall("southeast");
 savefig(['fig_', prefix, '_q']);
saveas(gcf,['fig_', prefix, '_q'],'epsc');