[SNR,val_SNR,val_SNR_approx]   = test_phi2sec(2);
[SNR,val_SNR1,val_SNR1_approx] = test_phi2sec(1);
[SNR,val_SNR2,val_SNR2_approx] = test_phi2sec(3);

plot(SNR,val_SNR,'Linewidth',2);
hold on; 
plot(SNR,val_SNR1,'--','Linewidth',2);
plot(SNR,val_SNR2,'-.','Linewidth',2);
plot(SNR,val_SNR_approx,'-s','Linewidth',2,'MarkerSize',10);
plot(SNR,val_SNR1_approx,'-o','Linewidth',2,'MarkerSize',10);
plot(SNR,val_SNR2_approx,'-d','Linewidth',2,'MarkerSize',10);

xlabel('input SNR','FontSize',20);
legend('$\frac{\| \hat q_f - \bar q_f \|_2}{\| \hat q_f\|_2}$, Signal Fig 1 (a)',...
       '$\frac{\| \hat q_f - \bar q_f \|_2}{\| \hat q_f\|_2}$, Signal Fig 1 (b)',...
       '$\frac{\| \hat q_f - \bar q_f \|_2}{\| \hat q_f\|_2}$, Signal Fig 1 (c)',...
       '$\frac{\| \hat q_f - \tilde q_f \|_2}{\| \hat q_f\|_2}$, Signal Fig 1 (a)',...
       '$\frac{\| \hat q_f - \tilde q_f \|_2}{\| \hat q_f\|_2}$, Signal Fig 1 (b)',...
       '$\frac{\| \hat q_f - \tilde q_f \|_2}{\| \hat q_f\|_2}$, Signal Fig 1 (c)');
ax = gca;
ax.FontSize = 20;
xlim([5 30]);
ylim([0 0.2]);
hold off;



