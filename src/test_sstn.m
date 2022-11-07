function test_sstn(cas,cas1,cas2)
 
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
  y =  exp(2*pi*i1*(A*t + B*t.^2 + C*t.^3 + D*t.^4));
  % x = cos(2*pi*250*t + 100*cos(2*pi*t));
  % x = cos(2*pi*100*t + 2*pi*150*t.^2);
  % x = cos(2*pi*300*t);

  phi1_ref = A + 2*B*t + 3*C*t.^2 + 4*D*t.^3;
  phi2_ref = 2*B + 6*C*t + 12*D*t.^2;
  phi3_ref = 6*C + 24*D*t;
  phi4_ref = 24*D*ones(size(t));
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
 [STFT0,STFT_thresh,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0] = sstn_det_simple_prec(s,sigma_opt,Nfft,gamma,cas1);
 
 
 s = s(:);
 nb_real = 30; %number of realizations
 SNR = 20:-5:-10;
 l_SNR =length(SNR);
 Count0 = zeros(nb_real,l_SNR);
 Count10 = zeros(nb_real,l_SNR);
 Count20 = zeros(nb_real,l_SNR);
 Count30 = zeros(nb_real,l_SNR);
 Count40 = zeros(nb_real,l_SNR);

 for k0 = 1:length(SNR),  
  k0
  for nb = 1:nb_real,    
   n  = randn(N,1)+1i*randn(N,1);
   [sn]  = sigmerge(s,n,SNR(k0));

   Nfft = N;
 
   gamma =10^(-2);
   %computation of SST,SST2, SST3 and SST4, and of the different reassignment operators
   [STFT0,STFT_thresh,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0] = sstn_det_simple_prec(sn,sigma_opt,Nfft,gamma,cas1);

   %we compute the proportion of good detection with the different
   %representation
 
   %count for good detection
   count0=0;
   count10=0;
   count20=0;
   count30=0;
   count40=0;
   
   %count for detection
   count=0;
   count1=0;
   count2=0;
   count3=0;
   count4=0;
   
   %We only consider interior points
   index = 50:N-50;
   
   for k = index,
    [a,ind0] = max(abs(STFT_thresh(:,k))); 
    if a > 0
     count = count+1;
    end
    
    [b,ind10] = max(abs(SSTd0(:,k)));
    if b > 0
     count1 = count1+1;
    end
    
    [c,ind20] = max(abs(SST2d0(:,k)));
    if c > 0
     count2 = count2+1;
    end 

    [d,ind30] = max(abs(SST3d0(:,k)));
    if d > 0
     count3 = count3+1;
    end 

    [e,ind40] = max(abs(SST4d0(:,k)));
    if e > 0
     count4 = count4+1;
    end 

     %first technique
     xx = round(if1(k)*Nfft/N)+1;
     if (abs(ind0-xx)<=20)
      count0 =count0+1;
     end
     %second technique
     if (abs(ind10-xx)<=20)
      count10 = count10+1;
     end
     if (abs(ind20-xx)<=20)
      count20 =count20+1;
     end
     if (abs(ind30-xx)<=20)
      count30 =count30+1;
     end
     if (abs(ind40-xx)<=20)
      count40 =count40+1;
     end
   end
   if cas2 == 1
    Count0(nb,k0) = count0/count;
    Count10(nb,k0) = count10/count1;
    Count20(nb,k0) = count20/count2;
    Count30(nb,k0) = count30/count3;
    Count40(nb,k0) = count40/count4;
   else
    LL = length(index);
    Count0(nb,k0) = count0/LL;
    Count10(nb,k0) = count10/LL;
    Count20(nb,k0) = count20/LL;
    Count30(nb,k0) = count30/LL;
    Count40(nb,k0) = count40/LL;  
   end
  end 
  mean_Count0(k0)  = mean(Count0(:,k0));
  mean_Count10(k0) = mean(Count10(:,k0));
  mean_Count20(k0) = mean(Count20(:,k0));
  mean_Count30(k0) = mean(Count30(:,k0));
  mean_Count40(k0) = mean(Count40(:,k0));
 end
 hold on;
 plot(SNR,mean_Count0);
 plot(SNR,mean_Count10,'--');
 plot(SNR,mean_Count20,'-.');
 plot(SNR,mean_Count30,':');
 plot(SNR,mean_Count40,'--s');
 hold off;
 
end