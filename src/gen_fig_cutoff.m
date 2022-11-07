
close all;

N = 1024;
t = (0:N-1)/N;

phi2_a = 400*ones(size(t));

B = 50;
C = 50;
D = 50;
phi2_b = 2*B + 6*C*t + 12*D*t.^2;
sigma_opt1 = 0.03;

f_sin = 10;
phi2_c = -8*exp(-2*(t-0.2)).*sin(f_sin*pi.*(t-0.2))+ 4*f_sin*pi*exp(-2*(t-0.2)).*cos(f_sin*pi.*(t-0.2))...
  +4*f_sin*pi.*exp(-2*(t-0.2)).*cos(f_sin*pi.*(t-0.2))+2*(f_sin*pi)^2*exp(-2*(t-0.2)).*sin(f_sin*pi.*(t-0.2));

sigma_opt2 = 0.01; %estaba en 0.01

FRACTION = .3162;
frec_corte_1 = sqrt(-lambertw(-FRACTION*exp(-1)))/(sigma_opt1*sqrt(pi));
frec_corte_2 = sqrt(-lambertw(-FRACTION*exp(-1)))/(sigma_opt2*sqrt(pi));

fc1 = frec_corte_1;
fc2 = frec_corte_2;

fft_GTa = fft(phi2_a);
fft_GTb = fft(phi2_b);
fft_GTc = fft(phi2_c);
f_vec = 0:(length(phi2_a) - 1);

plot_Ismall("frequency", "magnitude");
hold on;
plot(f_vec, abs(fft_GTa), "DisplayName", "DFT of $\phi''$, signal Fig. 1 (a)");
plot(f_vec, abs(fft_GTb), '-.', "DisplayName", "DFT of $\phi''$, signal Fig. 1 (b)");
Z = ylim;
plot([fc1, fc1], [Z(1), Z(2)], '--',...
 "DisplayName", "cut-off frequency");
hold off;

xlim([0, 2.5*fc1]);
legend_Ismall("northeast");
savefig(['fig_', 'LC_poly', '_cutoff']);
saveas(gcf,['fig_', 'LC_poly', '_cutoff'],'epsc');
% close all;

plot_Ismall("frequency", "magnitude");
hold on;
plot(f_vec, abs(fft_GTc), "DisplayName", "DFT of $\phi''$, signal Fig. 1 (c)");
Z = ylim;
plot([fc2, fc2], [Z(1), Z(2)], '--',...
 "DisplayName", "cut-off frequency",...
 "Color", [0.9290 0.6940 0.1250]);
hold off;

xlim([0, 2.5*fc2]);
legend_Ismall("northeast");
savefig(['fig_', 'sin', '_cutoff']);
saveas(gcf,['fig_', 'sin', '_cutoff'],'epsc');
% close all;