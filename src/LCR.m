function [modes, TFR_denoised, Lh] = LCR(STFT, IFs, CRs, g, sigma_s, L, Nfft, cas)
% [N_Y, L] = size(STFT);
[Nr, ~] = size(IFs);

[g, Lh] = create_gaussian_window(L, Nfft, sigma_s, 10^(-3));
% Lh = length(g)/2;


%% MR
TFR_denoised = zeros(size(STFT));
modes = zeros(Nr, L);

for p = 1:Nr
    %% use estimate and inverse STFT
    [TFR_denoised_r] = LCR_estim_STFT(sigma_s, STFT, IFs(p, :), CRs(p, :), Nfft);
    TFR_denoised = TFR_denoised + TFR_denoised_r;
    
    modes(p, :) = L*itfrstft(TFR_denoised_r, cas, g);
end

end
