 close all;

 cas = 1;
 N = 1024;
 t = (0:N-1)/N;
 a = 2;

 %% signal definition
 %different type of monocomponent signals 
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
  phi1 = 340*t-2.*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2));
  phi1_ref = 340+4*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2))-28*pi.*exp(-2*(t-0.2)).*cos(14*pi.*(t-0.2)); 
  phi2_ref = -8*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2))+ 56*pi*exp(-2*(t-0.2)).*cos(14*pi.*(t-0.2))...
      +56*pi.*exp(-2*(t-0.2)).*cos(14*pi.*(t-0.2))+28*14*pi^2*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2));
  s = a1.*exp(2*pi*1i*(phi1));
  sigma_opt = 0.01; %estaba en 0.01
 else
  return;
 end
 
 snr_in = 10;
 index = 100;
 Li = (index:N-index);
 %% TF parameters

for iter=1:5
 gamma =10^(-2);
 Nfft = N;
 n  = randn(N,1)+1i*randn(N,1);
 [sn]  = sigmerge(s(:),n,snr_in);%ac� estaba en 10dB

 [STFT0,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0,...
     phi2sec,phi2sec_simple,vtg,phi3sec,phi4sec] =...
     sstn_det(sn,sigma_opt,Nfft,gamma,0);

 global STFT_GT;
 [STFT_GT, ~] = sstn_det(s,sigma_opt,Nfft,gamma,0);

 [c,e] = exridge(SST2d0,0,0,10);
 X = zeros(size(t));
%  kk =index:N-index;
 for k=1:N
  X(k) = phi2sec_simple(c(1,k),k);
  X_phi1(k) = omega2d0(c(1,k),k);
 end 

 %% filter estimator q
 FRACTION = .1; %for 1% of maximum energy or .3162 for 10% of maximum energy
 frec_corte = sqrt(-lambertw(-FRACTION*exp(-1)))/(sigma_opt*sqrt(pi));

 d = designfilt('lowpassiir','FilterOrder',2, ...
         'PassbandFrequency',frec_corte/N,'PassbandRipple',0.2,'SampleRate',1); 
 y = filtfilt(d,[fliplr(X) X fliplr(X)]);
 y = y(N+1:end-N);

 %% apply LCR

 t = -0.5:1/N:0.5-1/N;
 t=t';
 sigma = sigma_opt;
 g =  1/sigma*exp(-pi/sigma^2*(t.^2));

 X = X(:);
 X = X';
 
    [modes, T_amp_old, Lh] = ...
        LCR_old(STFT0, phi1_ref, phi2_ref, g, sigma_opt, N, Nfft, 1);
    [modes, T_amp_new, Lh] = ...
        LCR(STFT0, phi1_ref, phi2_ref, g, sigma_opt, N, Nfft, 1);
    
    q_old = norm(abs(T_amp_old(:, Li) - STFT_GT(:, Li)), 'fro');
    q_new = norm(abs(T_amp_new(:, Li) - STFT_GT(:, Li)), 'fro');
    
    fprintf("old = %f, new = %f, diff = %f\n", q_old, q_new, q_new - q_old);
end
 return;

%  return;
 [modes, T_LCR_nf, Lh] = ...
     LCR(STFT0, X_phi1, X, g, sigma_opt, N, Nfft, 1);

 [modes, T_LCR_nf2, Lh] = ...
     LCR(STFT0, phi1_ref, X, g, sigma_opt, N, Nfft, 1);

 [modes, T_LCR_filt, Lh] = ...
     LCR(STFT0, X_phi1, y, g, sigma_opt, N, Nfft, 1);

 [modes, T_LCR_filt2, Lh] = ...
     LCR(STFT0, phi1_ref, y, g, sigma_opt, N, Nfft, 1);
 
norm(abs(T_LCR_filt(:, Li) - STFT_GT(:, Li)), 'fro')

norm(abs(T_LCR_filt2(:, Li) - STFT_GT(:, Li)), 'fro')

norm(abs(T_LCR_nf(:, Li) - STFT_GT(:, Li)), 'fro')

norm(abs(T_LCR_nf2(:, Li) - STFT_GT(:, Li)), 'fro')
 
 %% figures
 figure 
 imagesc(t(1:N),0:Nfft/2-1,abs(STFT0));
 xlabel('time');
 ylabel('frequency');
 set(gca,'ydir','normal');
 hold on
 plot(t(1:N),c(1,:),'LineWidth',2);
 hold off;
 
 figure 
 imagesc(t(1:N),0:Nfft/2-1,abs(T_LCR_filt2));
 xlabel('time');
 ylabel('frequency');
 set(gca,'ydir','normal');
 hold on
 plot(t(1:N),c(1,:),'LineWidth',2);
 hold off;
 return;
 
 figure;
 hold on;
 plot(t(index:N-index),X(index:N-index))
 plot(t(index:N-index),phi2_ref(index:N-index))
 plot(t(index:N-index),y(index:N-index))
 hold off;
 legend({'$\bar{q}$','ground truth','filtered $\bar{q}$'},'Interpreter','latex','FontSize',30,'location','best')
 ylabel('estimation of $\phi^{(2)} $ ','Interpreter','latex','FontSize',30);
 xlabel('time','FontSize',30);
  
 figure 
 imagesc(t(Li),0:Nfft/2-1,...
     abs(T_LCR_nf(:, Li) - STFT_GT(:, Li)));
 xlabel('time');
 ylabel('frequency');
 set(gca,'ydir','normal');
%  hold on
%  plot(t(Li),c(1,:),'LineWidth',2);
%  hold off;
 title("D nofilt GT");
  
 figure 
 imagesc(t(index:N-index),0:Nfft/2-1,...
     abs(T_LCR_filt(:, Li) - STFT_GT(:, Li)));
 xlabel('time');
 ylabel('frequency');
 set(gca,'ydir','normal');
%  hold on
%  plot(t(Li),c(1,:),'LineWidth',2);
%  hold off;
 title("D filt GT");
 
 norm(abs(T_LCR_filt(:, Li) - STFT_GT(:, Li)), 'fro')
 
 norm(abs(T_LCR_nf(:, Li) - STFT_GT(:, Li)), 'fro')
 
%  figure 
%  imagesc(t(1:N),0:Nfft/2-1,abs(T_LCR_2 - T_LCR_1));
%  xlabel('time');
%  ylabel('frequency');
%  set(gca,'ydir','normal');
%  hold on
%  plot(t(1:N),c(1,:),'LineWidth',2);
%  hold off;
%  title("D filt nofilt");
 
%  figure 
%  imagesc(t(1:N),0:Nfft/2-1,abs(STFT_GT));
%  xlabel('time');
%  ylabel('frequency');
%  set(gca,'ydir','normal');
%  hold on
%  plot(t(1:N),c(1,:),'LineWidth',2);
%  hold off;
%  title("Ground truth");