close all
 N = 1024;
 t = (0:N-1)/N;
 a  = 2;
 q = 50;
 s  = a.*exp(2*pi*1i*(100*t+q*t.^2));
 %s  = a.*exp(2*pi*1i*(100*t));
 phi1_ref = 100 + 2*q*t;
 phi2_ref = 2*q*ones(size(t));
 
 gamma =10^(-2);
 Nfft = N;
 n  = randn(N,1)+1i*randn(N,1);
 [sn]  = sigmerge(s(:),n,30);
 
 %sigma maximizing F
 %sigma_1 = 1/(5^(1/4)*sqrt(400));
 sigma_1 = 0.01;
 [STFT0,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0,phi2sec,phi2sec_simple,vtg,phi3sec,phi4sec] = sstn_det(sn,sigma_1,Nfft,gamma,0);
 
 [omega,omega_simple,omega2,omega2_simple] = omega_comput(sn,s(:),sn-s(:),sigma_1,Nfft);


 figure;
 [c,e] = exridge(SST2d0,0,0,10);
 X1 =[];
 index = 100;
 kk =index:N-index;
 for k=kk
  X1 = [X1 omega(c(1,k),k)];
 end 
 plot(t(index:N-index),X1,'LineWidth',2);
 hold on;
 kk =index:N-index;
 X =[];
 for k=kk
  X = [X omega_simple(c(1,k),k)];
 end
 plot(t(index:N-index),X,'-.','LineWidth',2);
 %plot(t(index:N-index),c(1,index:N-index)-1,'-o','LineWidth',2);
 plot(t(index:N-index),phi1_ref(index:N-index),'-s','LineWidth',2);  
 legend({'$\hat{w}$','$\bar{w}$','ground truth'},'Interpreter','latex','FontSize',30)
 xlabel('time','FontSize',30);
 xlim([0.3 0.5]);
 %axis off;
 hold off;

 figure;
 [c,e] = exridge(SST2d0,0,0,10);
 X1 =[];
 index = 100;
 kk =index:N-index;
 for k=kk
  X1 = [X1 omega2(c(1,k),k)];
 end 
 plot(t(index:N-index),X1,'LineWidth',2);
 hold on;
 kk =index:N-index;
 X =[];
 for k=kk
  X = [X omega2_simple(c(1,k),k)];
 end
 plot(t(index:N-index),X,'-.','LineWidth',2);
 %plot(t(index:N-index),c(1,index:N-index)-1,'-o','LineWidth',2);
 plot(t(index:N-index),phi1_ref(index:N-index),'-s','LineWidth',2);  
 legend({'$\hat{w}$','$\bar{w}$','ground truth'},'Interpreter','latex','FontSize',30)
 xlabel('time','FontSize',30);
 xlim([0.3 0.5]);
 %axis off;
 hold off;

