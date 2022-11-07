function [SNR,val_SNR,val_SNR1] = test_phi2sec(cas)
 
%cas: for the type of signal
%cas1: for noise removal
%cas2: proportion detection over the whole signal (0 else 1)
% close all;
 N = 1024;
 t = (0:N-1)/N;
 a  = 2;

 %different type of monocomponent signals 
 % test signal -------------------------
 if cas == 0,
  A = 50;
  B = 50;
  C = 50;
  D = 50;
  s =  exp(2*pi*1i*(A*t + B*t.^2 + C*t.^3 + D*t.^4));
  sigma_opt = 0.03;
%-------------------------------------
 elseif cas == 1,
  s   = a.*exp(2*pi*1i*(50*t+50*t.^2+50*t.^3+50*t.^4));
  phi1_ref = 50+100*t+150*t.^2+200*t.^3;
  phi2_ref = 100+ 300*t+600*t.^2;
  phi3_ref = 300+ 1200 * t;
  phi4_ref = 1200*ones(size(t));
  sigma_opt = 0.03;
 elseif cas == 2,  
  s   = a.*exp(2*pi*1i*(100*t+200*t.^2));
  phi1_ref = 100 + 400*t;
  phi2_ref = 400*ones(size(t));
  sigma_opt = 0.03;
 else  
  a1 = 1+ 5*t.^3 + 7*(1-t).^6;
  phi1 = 340*t-2.*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2));
  phi1_ref = 340+4*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2))-28*pi.*exp(-2*(t-0.2)).*cos(14*pi.*(t-0.2)); 
  phi2_ref = -8*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2))+ 56*pi*exp(-2*(t-0.2)).*cos(14*pi.*(t-0.2))...
      +56*pi.*exp(-2*(t-0.2)).*cos(14*pi.*(t-0.2))+28*14*pi^2*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2));
  s = a1.*exp(2*pi*1i*(phi1));
  sigma_opt = 0.01;
 end
 
 Nfft = N;
 SNR = 0:30;
 nbreal = 20;
 val = zeros(nbreal,length(SNR));
 val1 = zeros(nbreal,length(SNR));
 
 for k = 1:length(SNR)
  for nb = 1:nbreal   
   n  = randn(N,1)+1i*randn(N,1);
   [sn]  = sigmerge(s(:),n,SNR(k));
   [STFT0,phi2sec,phi2sec_simple,~] = compute_phi2sec(sn,zeros(N,1),sigma_opt,Nfft);
   [~,~,phi2sec_simple_nonoise,extra_term] = compute_phi2sec(s,sn-s,sigma_opt,Nfft);
 
   [c,~] = exridge(STFT0,0,0,10);
   index = 100;
   kk =index:N-index;
   X = [];
   Y = [];
   Z = [];
   for p=kk
    X = [X phi2sec(c(1,p),p)];
    Y = [Y phi2sec_simple(c(1,p),p)];
    Z = [Z phi2sec_simple_nonoise(c(1,p),p)+extra_term(c(1,p),p)];
   end
   val(nb,k)  = norm(X-Y)/norm(X);
   val1(nb,k) = norm(X-Z)/norm(X);
  end
  val_SNR  = mean(val); 
  val_SNR1 = mean(val1);
 end