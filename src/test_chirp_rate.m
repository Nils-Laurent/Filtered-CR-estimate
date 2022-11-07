close all
 N = 1024;
 t = (0:N-1)/N;
 a  = 2;
 q = 0;
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
 
 [phi2sec,phi2sec_simple,phi2sec_simple2,val_car] = chirp_rate_comput(sn,s(:),sn-s(:),sigma_1,Nfft);


 figure;
 [c,e] = exridge(SST2d0,0,0,10);
 X =[];
 index = 100;
 kk =index:N-index;
 for k=kk
  X = [X phi2sec(c(1,k),k)];
 end 
 plot(t(index:N-index),X,'LineWidth',2);
 hold on;
 kk =index:N-index;
 X =[];
 for k=kk
  X = [X phi2sec_simple(c(1,k),k)];
 end
 plot(t(index:N-index),X,'-.','LineWidth',2);
 X =[];
 val_car_real = [];
 val_car_imag = [];
 val_car_modulus = [];
 for k=kk
  X = [X phi2sec_simple2(c(1,k),k)];
  val_car_real = [val_car_real real(val_car(c(1,k),k))];
  val_car_imag = [val_car_imag imag(val_car(c(1,k),k))];
  val_car_modulus = [val_car_modulus abs(val_car(c(1,k),k))];
 end 
 plot(t(index:N-index),X,'--','LineWidth',2);
 
 plot(t(index:N-index),phi2_ref(index:N-index),'-s','LineWidth',2);  
 legend({'$\hat{q}$','$\bar{q}$','approx $\bar{q}$','ground truth'},'Interpreter','latex','FontSize',30)
 ylabel('estimation of $\phi^{(2)} $ ','Interpreter','latex','FontSize',30);
 xlabel('time','FontSize',30);
 xlim([0.3 0.5]);
 %axis off;
 hold off;
 
 figure;
 %plot(t(index:N-index),val_car_real,'-s','LineWidth',2);
 %hold on; 
 %plot(t(index:N-index),val_car_imag,'-o','LineWidth',2);
  plot(t(index:N-index),val_car_imag,'-o','LineWidth',2);
 
 figure
 plot(abs(fft(val_car_imag)))
 hold on;
  plot(abs(fft(val_car_real)),'--')
  hold off
  %hold off;

 