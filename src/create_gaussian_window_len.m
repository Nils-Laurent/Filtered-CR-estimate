function [g, Lh, sigma_w] = create_gaussian_window_len(L, Nfft, Lw_i, prec)

Lh = floor(Lw_i/2);
sigma_w = (Lh - 1)/(L*sqrt(-log(prec)/pi));

[g, Lh] = create_gaussian_window(L, Nfft, sigma_w, prec);

% t=(1:Lw)'/Fs - (Lh + 1)/Fs;
% g = exp(-(t/sigma_w).^2 * pi);
% g = g/sum(g);
% 
% if 2*Lh + 1 > Nfft
%     fprintf("[Warning] 2*Lh+1 > Nfft, simga_w = %f\n", sigma_w);
% end
% 
% 
% [~, Lh2] = create_gaussian_window(Fs, Nfft, sigma_w, prec);
% if Lh ~= Lh2
%     error("Lh ~= Lh2");
% end

end