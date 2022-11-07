 close all;
 
% cas: for the type of signal
% cas1: for noise removal
% cas2: proportion detection over the whole signal (0 else 1)

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
 
 Nrand = 3;
 SNRs = [4, 6, 8, 10, 15, 20];
 Nsnr = length(SNRs);
 R_data = zeros(Nrand, Nsnr);
 Li = (index:N-index);
 for n_snr = 1:Nsnr
     snr_in = SNRs(n_snr);
     snr_in
     for n_rand=1:Nrand
         %% TF parameters
         gamma =10^(-2);
         Nfft = N;
         n  = randn(N,1)+1i*randn(N,1);
         [sn]  = sigmerge(s(:),n,snr_in);

         [STFT0,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0,...
             phi2sec,phi2sec_simple,vtg,phi3sec,phi4sec] =...
             sstn_det(sn,sigma_opt,Nfft,gamma,0);

         [STFT_GT, ~] = sstn_det(s,sigma_opt,Nfft,gamma,0);

         [c,e] = exridge(SST2d0,0,0,10);
         X = zeros(size(t));
         index = 100;
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
         
         R_data(n_rand, n_snr) = snr(phi2_ref(Li), y(Li));
     end
 end
 
 %% test filt vs nofilt
 
 figure;
 boxplot(R_data, SNRs);