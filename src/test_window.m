 close all
 N = 1024;
 t = (0:N-1)/N;
 a  = 2;
 s  = a.*exp(2*pi*1i*(100*t+200*t.^2));
 %s  = a.*exp(2*pi*1i*(100*t));
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
 
 %sigma maximizing F
 sigma_1 = 1/(5^(1/4)*sqrt(400));
 
 %sigma corresponding to the Rényi entropy
 sigma_ren = 1/sqrt(400);
 
 %we test the simple estimation for hat q
 [STFT0,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0,phi2sec,phi2sec_simple,vtg,phi3sec,phi4sec] = sstn_det(sn,sigma_ren,Nfft,gamma,0);
  
 %we compute q on the ridge
 [Cs,e] = exridge(STFT0,0,0,10); 
 Y =[];
 for k = 1:N
  Y = [Y phi2sec_simple(Cs(1,k),k)];
 end
 
 %we compute the extrema of q
 [indmin, indmax, indzer] = extr(Y);
 X = sort([indmin,indmax]);
 Y1 = [];
 for k = 1:length(X),
  Y1 = [Y1 phi2sec(Cs(1,X(k)),X(k))];
 end
 pp = csaps(t(X),Y1);
 
 index = 100;
 abs_spline = t(index):0.0001:t(N-index);
 val_pp = fnval(abs_spline,pp);
 [indmin_s, indmax_s, indzer] = extr(val_pp);
 figure 
 plot(abs_spline(indmax_s),val_pp(indmax_s),'*','MarkerSize',10);
 hold on;
 plot(abs_spline(indmin_s),val_pp(indmin_s),'o','MarkerSize',10);
 plot(t(index):0.0001:t(N-index),val_pp,'LineWidth',2);
 plot(t(index:N-index),Y(index:N-index),'--','LineWidth',2)
 hold off
 
 figure
 %computation of instantaneous frequency
 len = min(length(indmax_s),length(indmin_s));
 freq = [];
 if indmax_s(1) > indmin_s(1)
  for k = 1:len-1  
   freq = [freq ones(1,indmax_s(k)-indmin_s(k))/(2*(abs_spline(indmax_s(k))-abs_spline(indmin_s(k))))];
   freq = [freq ones(1,indmin_s(k+1)-indmax_s(k))/(2*(abs_spline(indmin_s(k+1))-abs_spline(indmax_s(k))))];
  end
  freq = [freq 1/(2*(abs_spline(indmin_s(len))-abs_spline(indmax_s(len-1))))];
  plot(abs_spline(indmin_s(1):indmin_s(len)),freq)
 else
  for k = 1:len-1 
   freq = [freq ones(1,indmin_s(k)-indmax_s(k))/(2*(abs_spline(indmin_s(k))-abs_spline(indmax_s(k))))];
   freq = [freq ones(1,indmax_s(k+1)-indmin_s(k))/(2*(abs_spline(indmax_s(k+1))-abs_spline(indmin_s(k))))];
  end
  freq = [freq 1/(2*(abs_spline(indmax_s(len))-abs_spline(indmax_s(len-1))))];
  plot(abs_spline(indmax_s(1):indmax_s(len)),freq) 
 end
 
 %we suppress useless points
  
 


 %
 [STFT0,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0,phi2sec,phi2sec_simple,vtg,phi3sec,phi4sec] = sstn_det(sn,0.1,Nfft,gamma,0);
 
 [Cs,e] = exridge(STFT0,0,0,10); 
 X = [];

 for k=1:N 
  X = [X abs(STFT0(Cs(1,k),k))];
 end
 plot(t(100:N-100),X(100:N-100))
 pause
 %[Loc]  = peakfinder(X);

 figure 
 imagesc(t(1:N),0:Nfft/2-1,abs(STFT0));
 xlim([0.3 0.5]);
 xlabel('time');
 ylabel('frequency');
 set(gca,'ydir','normal');
 hold on
 plot(t(1:N),Cs(1,:),'LineWidth',2);
 hold off;

 figure;
 X =[];
 Y =[];
 %index = 100;
 %kk =index:N-index;
 
 for k = 1:N
  X = [X phi2sec(Cs(1,k),k)];
 end
 [indmin, indmax, indzer] = extr(X);
 %indmax = peakfinder(X);
 %indmin = peakfinder(-X);
 
 for k = 1:length(indmax)
  Y = [Y phi2sec(Cs(1,indmax(k)),indmax(k))];
 end 
 plot(t,X,'LineWidth',2);
 hold on;
 plot(t(indmax),Y,'*','MarkerSize',10);

 Z= [];
 for k = 1:length(indmin)
  Z = [Z phi2sec(Cs(1,indmin(k)),indmin(k))];
 end 
 plot(t,X,'LineWidth',2);
 hold on;
 plot(t(indmin),Z,'o','MarkerSize',10);
 
 %  kk =index:N-index;
%  X =[];
%  for k=kk
%   X = [X phi2sec_simple(Cs(1,k),k)];
%  end 
 
 %plot(t(index:N-index),X,'-.','LineWidth',2);
 
 plot(t,phi2_ref,'LineWidth',2);  
 %legend({'$\hat{q}$','$\bar{q}$','ground truth'},'Interpreter','latex','FontSize',30)
 ylabel('estimation of $\phi^{(2)} $ ','Interpreter','latex','FontSize',30);
 xlabel('time','FontSize',30);
 xlim([0.3 0.5]);
 hold off;
 
 figure
 len = min(length(indmax),length(indmin));
 freq = [];
 if t(indmax(1)) > t(indmin(1))
  for k = 1:len  
   freq = [freq ones(1,indmax(k)-indmin(k)-1)*1/(2*(t(indmax(k))-t(indmin(k))))];
  end
  plot(freq)
 else
  for k = 1:len 
   freq = [freq ones(1,indmin(k)-indmax(k)-1)*1/(2*(t(indmin(k))-t(indmax(k))))];
  end
  length(t(indmin(1):indmax(end)))
  length(freq)
  plot(freq) 
 end  

 figure;
 X =[];
 Y =[];
 index = 100;
 
 for k = 1:length(Loc)
  Y = [Y abs(SST2d0(Cs(1,Loc(k)),Loc(k)))];
 end 
 for k = 1:N
  X = [X abs(SST2d0(Cs(1,k),k))];
 end
 plot(t(index:N-index),X(index:N-index),'LineWidth',2);
 xlim([0.3 0.5]);
 hold on;
 plot(t(Loc),Y,'*','MarkerSize',10);
 hold off;

 figure;
 X =[];
 Y =[];
 index = 100;
 kk =index:N-index;
 
 for k = kk
  X = [X omega2d0(Cs(1,k),k)];
 end
 
 for k = 1:length(Loc)
  Y = [Y omega2d0(Cs(1,Loc(k)),Loc(k))];
 end 
 plot(t(index:N-index),X,'LineWidth',2);
 plot(t(Loc),Y,'*','MarkerSize',10);
 

 hold on;
%  kk =index:N-index;
%  X =[];
%  for k=kk
%   X = [X phi2sec_simple(Cs(1,k),k)];
%  end 
 
 plot(t(index:N-index),X,'-.','LineWidth',2);
 
 plot(t(index:N-index),phi1_ref(index:N-index),'LineWidth',2);  
 %legend({'$\hat{q}$','$\bar{q}$','ground truth'},'Interpreter','latex','FontSize',30)
 ylabel('estimation of $\phi^{(1)}$','Interpreter','latex','FontSize',30);
 xlabel('time','FontSize',30);
 xlim([0.3 0.5]);
 hold off;
 
 pause
 
 [STFT0,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0,phi2sec,phi2sec_simple,vtg,phi3sec,phi4sec] = sstn_det(sn,0.05,Nfft,gamma,0);
 
 [Cs,e] = exridge(SST2d0,0,0,10); 
X = [];
for k=1:N 
 X = [X abs(SST2d0(Cs(1,k),k))];
end

[Loc]  = peakfinder(X);

figure 
imagesc(t(1:N),0:Nfft/2-1,abs(STFT0));
xlabel('time');
ylabel('frequency');
set(gca,'ydir','normal');
hold on
plot(t(1:N),Cs(1,:),'LineWidth',2);
hold off;

 figure;
 X =[];
 Y =[];
 index = 100;
 kk =index:N-index;
 
 for k = kk
  X = [X phi2sec(Cs(1,k),k)];
 end
 
 for k = 1:length(Loc)
  Y = [Y phi2sec(Cs(1,Loc(k)),Loc(k))];
 end 
 plot(t(index:N-index),X,'LineWidth',2);
 plot(t(Loc),Y,'*','MarkerSize',10);
 

 hold on;
 kk =index:N-index;
 X =[];
 for k=kk
  X = [X phi2sec_simple(Cs(1,k),k)];
 end 
 
 plot(t(index:N-index),X,'-.','LineWidth',2);
 
 plot(t(index:N-index),phi2_ref(index:N-index),'-s','LineWidth',2);  
 %legend({'$\hat{q}$','$\bar{q}$','ground truth'},'Interpreter','latex','FontSize',30)
 ylabel('estimation of $\phi^{(2)} $ ','Interpreter','latex','FontSize',30);
 xlabel('time','FontSize',30);
 xlim([0.3 0.5]);
 hold off;

 [STFT0,SSTd0,SST2d0,SST3d0,SST4d0,omegad0,omega2d0,omega3d0,omega4d0,phi2sec,phi2sec_simple,vtg,phi3sec,phi4sec] = sstn_det(sn,0.1,Nfft,gamma,0);
 
[Cs,e] = exridge(SST2d0,0); 
X = [];
for k=1:N 
 X = [X abs(SST2d0(Cs(1,k),k))];
end

[Loc]  = peakfinder(X);

figure 
imagesc(t(1:N),0:Nfft/2-1,abs(STFT0));
xlabel('time');
ylabel('frequency');
set(gca,'ydir','normal');
hold on
plot(t(1:N),Cs(1,:),'LineWidth',2);
hold off;

 figure;
 X =[];
 Y =[];
 index = 100;
 kk =index:N-index;
 
 for k = kk
  X = [X phi2sec(Cs(1,k),k)];
 end
 
 for k = 1:length(Loc)
  Y = [Y phi2sec(Cs(1,Loc(k)),Loc(k))];
 end 
 plot(t(index:N-index),X,'LineWidth',2);
 plot(t(Loc),Y,'*','MarkerSize',10);
 

 hold on;
 kk =index:N-index;
 X =[];
 for k=kk
  X = [X phi2sec_simple(Cs(1,k),k)];
 end 
 
 plot(t(index:N-index),X,'-.','LineWidth',2);
 
 plot(t(index:N-index),phi2_ref(index:N-index),'-s','LineWidth',2);  
 %legend({'$\hat{q}$','$\bar{q}$','ground truth'},'Interpreter','latex','FontSize',30)
 ylabel('estimation of $\phi^{(2)} $ ','Interpreter','latex','FontSize',30);
 xlabel('time','FontSize',30);
 xlim([0.3 0.5]);
 hold off;
