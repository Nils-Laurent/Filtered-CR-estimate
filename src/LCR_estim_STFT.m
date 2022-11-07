function [TFR_denoised] = LCR_estim_STFT(sigma_s, STFT, phipE, phippE, Nfft)

[N_STFT, L] = size(STFT);
TFR_denoised = zeros(size(STFT));
for n = 1:L
    Zn = pi*sigma_s^2*(1 + 1i*phippE(n)*sigma_s^2)...
        /(1 + phippE(n)^2*sigma_s^4);
    
    K0 = round(phipE(n)*Nfft/L)+1;
    K0 = min(N_STFT, K0);
    K0 = max(1, K0);
    
    F0 = ((K0-1)*L/Nfft);
    Gn = STFT(K0, n)*exp(Zn*(F0 - phipE(n))^2);
    
    %% assign coefficients
    Model = zeros(1, N_STFT);
%     Y = X;
%     Y = Gn*normpdf((0:N_STFT-1)*L/Nfft, phipE(n), Zn);
    
    for k = 1:N_STFT
        Model(k) = exp(-Zn*((k-1)*L/Nfft - phipE(n))^2);
    end
    
    TFR_denoised(:, n) = Gn*Model;
    
end

end

