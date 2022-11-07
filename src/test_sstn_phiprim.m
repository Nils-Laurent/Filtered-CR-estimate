function test_sstn_phiprim(cas,cas1,cas2)
 
%cas: for the type of signal
%cas1: for noise removal
%cas2: proportion detection over the whole signal (0 else 1)

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
  % x = cos(2*pi*250*t + 100*cos(2*pi*t));
  % x = cos(2*pi*100*t + 2*pi*150*t.^2);
  % x = cos(2*pi*300*t);

  phi1_ref = A + 2*B*t + 3*C*t.^2 + 4*D*t.^3;
  phi2_ref = 2*B + 6*C*t + 12*D*t.^2;
  phi3_ref = 6*C + 24*D*t;
  phi4_ref = 24*D*ones(size(t));
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
 
 gamma =10^(-2);
 Nfft = N;
  
 s = s(:);
 SNR = 10;
 n  = randn(N,1)+1i*randn(N,1);
 [sn]  = sigmerge(s,n,SNR);

 Nfft = N;
 
 gamma =10^(-2);
 %computation of SST,SST2, SST3 and SST4, and of the different reassignment operators
 [STFT0,STFT_thresh,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0] = sstn_det_simple_prec(sn,sigma_opt,Nfft,gamma,cas1);
 
 figure 
 imagesc(abs(STFT0));
 set(gca,'ydir','normal');
 figure;
 
 [c,e] = exridge(SSTd0,0,0,10);
 k = 500;    
 plot(-5:5,omegad0(c(1,k)-5:c(1,k)+5,k));
 hold on;
 [c,e] = exridge(SST2d0,0,0,10);   
 plot(-5:5,omega2d0(c(1,k)-5:c(1,k)+5,k),'--');
 [c,e] = exridge(SST3d0,0,0,10);
 plot(-5:5,omega3d0(c(1,k)-5:c(1,k)+5,k),'-.');
 [c,e] = exridge(SST4d0,0,0,10);
 plot(-5:5,omega4d0(c(1,k)-5:c(1,k)+5,k),':');
 
 plot(-5:5,phi1_ref(k)*ones(1,11),'-s');  
 hold off;
end