function [TFR_denoised] = LCR_estim_STFT_old(sigma_s, STFT, phipE, phippE, Nfft)

global STFT_GT;

[N_STFT, L] = size(STFT);
TFR_denoised = zeros(size(STFT));
TFR2 = zeros(size(STFT));
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
    
%     LSEn = dot(Model, STFT(:, n))/dot(Model, Model);
%     LSEn = STFT_GT(K0, n)*exp(Zn*(F0 - phipE(n))^2);
    TFR_denoised(:, n) = Gn*Model;
    
%     if n == L/2
%         N = STFT(:, n);
%         T = STFT_GT(:, n);
%         X = TFR_denoised(:, n);
%         Y = LSEn*Model;
%         Y = Y.';
%         
%         figure;
%         hold on;
%         plot(real(T), 'k-', 'DisplayName', 'Ground truth');
%         plot(real(X), 'r--', 'DisplayName', 'old m.');
%         plot(real(Y), 'g--', 'DisplayName', 'new m.');
%         hold off;
%         legend;
%         
%         figure;
%         hold on;
%         plot(real(T - X), 'r--', 'DisplayName', 'old m.');
%         plot(real(T - Y), 'g--', 'DisplayName', 'new m.');
%         hold off;
%         figure;
%         hold on;
%         plot(imag(T - X), 'r--', 'DisplayName', 'old m.');
%         plot(imag(T - Y), 'g--', 'DisplayName', 'new m.');
%         hold off;
%         
%         norm(STFT_GT(:, n) - N)
%         norm(STFT_GT(:, n) - X)
%         norm(STFT_GT(:, n) - Y)
%     end
end

end

