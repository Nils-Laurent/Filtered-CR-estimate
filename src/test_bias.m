 close all
 N = 1024;
 t = (0:N-1)/N;
 a  = 2;
 s  = a.*exp(2*pi*1i*(100*t+200*t.^2));
 phi1_ref = 100 + 400*t;
 phi2_ref = 400*ones(size(t));
 
 gamma =10^(-2);
 Nfft = N;
 n  = randn(N,1)+1i*randn(N,1);
 sigman_2 =[];

 for k =5:30 
  [sn]  = sigmerge(s(:),n,k);
  %we compute the variance of the noise 
  sigman_2 = [sigman_2 var(sn-s(:))];
 end
 SNR = 5:30;
 
 %sigma_1 minimizing the bias on the ridge
 sigma_1 = 1/(5^(1/4)*sqrt(400));
 val = (a^2*sigma_1*4*sqrt(2))/(3*(1+(400*sigma_1^2)^2)^(3/2)); 
 gamma =10^(-2);
 Nfft = N;
 
 sigma_opt = 0.03; 
 [sn]  = sigmerge(s(:),n,10);
 [STFT0,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0,phi2sec,phi2sec_simple,vt2g,phi3sec,phi4sec] = sstn_det(s,sigma_opt,Nfft,gamma,0);
 [STFT0,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0,phi2sec,phi2sec_simple,vt2gn,phi3sec,phi4sec] = sstn_det(sn-s(:),sigma_opt,Nfft,gamma,0);
 [STFT0,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0,phi2sec,phi2sec_simple,vt2gn,phi3sec,phi4sec] = sstn_det(n(:),sigma_opt,Nfft,gamma,0);
 %plot(1:Nfft,abs(vt2g(:,300)).^2,1:Nfft,abs(vt2gn(:,300).^2))
 sigman = var(n)
 mean(mean(abs(vt2gn(:,:).^2)))
 3*sigma_opt^3/(4*pi^2*4*sqrt(2))*sigman
 pause
 [c,e] = exridge(STFT0,0,0,10);
 Y = abs(STFT0);
 Z = abs(vt2g);
 X  = [];
 X1 = [];
 for k = 1:N
   X  = [X Y(c(1,k),k)];
   X1 = [X1 Z(c(1,k),k)];
 end
 [sn]  = sigmerge(s(:),n,0);
 [STFT0n,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0,phi2sec,phi2sec_simple,vt2gn,phi3sec,phi4sec] = sstn_det(sn,sigma_opt,Nfft,gamma,0);
 [c1,e] = exridge(STFT0n,0,0,10);
 Y = abs(STFT0n);
 Z = abs(vt2gn);
 X2  = [];
 X3 = [];
 for k = 1:N
   X2  = [X2 Y(c(1,k),k)];
   X3 = [X3 Z(c(1,k),k)];
 end
 
 Y = abs(STFT0n);
 Z = abs(vt2gn);
 X4  = [];
 X5 = [];
 for k = 1:N
   X4  = [X4 Y(c1(1,k),k)];
   X5 = [X5 Z(c1(1,k),k)];
 end
 figure
 plot(1:N,X,1:N,X2,'--',1:N,X4,'-.')
  
 figure
 plot(1:N,X1,1:N,X3,'--',1:N,X5,'-.')
 pause
 %plot(1:N,X,1:N,2/(1+(400*sigma_opt^2)^2)^(1/4)*ones(1,N),'--')
 %figure;
 plot(SNR,4*sigma_opt*4*sqrt(2)/(3*(1+(400*sigma_opt^2)^2)^(3/2))*ones(1,length(SNR)),'--')
 pause
 
 plot(5:30,exp(-ones(1,26)*2*a^2./(3*sqrt(400)*sigman_2)),5:30,exp(-ones(1,26)*val./sigman_2),'--','LineWidth',2);  
 legend({'$e^{-F(\sigma_r)/\sigma_n^2}$','$e^{-F(\hat \sigma)/\sigma_n^2}$'},'Interpreter','latex','FontSize',30)
 ax.fontsize = 30;
 %ylabel('estimation of $\phi^{(2)} $ ','Interpreter','latex','FontSize',30);
 xlabel('SNR','FontSize',30);
 