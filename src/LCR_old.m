function [modes, TFR_denoised, Lh] = LCR_old(STFT, IFs, CRs, g, sigma_s, L, Nfft, cas)
% [N_Y, L] = size(STFT);
[Nr, ~] = size(IFs);

% [g, Lg] = create_gaussian_window(L, Nfft, sigma_s, 10^(-3));
Lh = length(g)/2;
% g = g/sigma_s;


%% MR
TFR_denoised = zeros(size(STFT));
modes = zeros(Nr, L);

for p = 1:Nr
    %% use estimate and inverse STFT
    [TFR_denoised_r] = LCR_estim_STFT_old(sigma_s, STFT, IFs(p, :), CRs(p, :), Nfft);
    TFR_denoised = TFR_denoised + TFR_denoised_r;
    % modes(p, :) = L*itfrstft(TFR_denoised_r, cas, g);
%     modes(p, :) = FM_inverse(TFR_denoised_r, Fs, Nfft, g, cas);
end

end
