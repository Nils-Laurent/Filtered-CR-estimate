 close all
 N = 1024;
 t = (0:N-1)/N;
 a  = 2;
 s  = a.*exp(2*pi*1i*(100*t+200*t.^2));
 phi1_ref = 100 + 400*t;
 phi2_ref = 400*ones(size(t));
 
 gamma =10^(-2);
 Nfft = N;
 SNR = 0:20;
 nb_real = 1;
 val_SNR = zeros(1,length(SNR));
 for k=1:length(SNR)
  k
  val =[];
  for p = 1:nb_real,
  n  = randn(N,1)+1i*randn(N,1);
  sigma_opt = 0.03; 
  [sn]  = sigmerge(s(:),n,SNR(k));
  [STFT0,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0,phi2sec,phi2sec_simple,vt2gn,phi3sec,phi4sec] = sstn_det(sn-s(:),sigma_opt,Nfft,gamma,0);
  val_ridge = a^2*sigma_opt^4/(4*pi^2*(1+(400*sigma_opt^2)^2)^(3/2));
  val_noise = mean(mean(abs(vt2gn(:,:).^2)));
  val = [val val_ridge/val_noise];
  end
  val_SNR(k) = mean(val); 
 end
  