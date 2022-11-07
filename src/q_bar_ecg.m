function  [STFT,omega,omega2,phi2_hat,phi2_bar] = q_bar_ecg(s,sigma_w,Nfft,max_f)

%% sstn_new : computes the STFT of a signal and different versions of synchrosqueezing
% adding threshold conditions for denominators of modulation operators q
% INPUTS:   
%   s: real or complex signal
%   sigma: the variance of the Gaussian window
%   Nfft: number of frequency bins
%   gamma: threshold on the STFT for reassignment
% OUTPUTS:   
%   STFT: the short-time Fourier transform
%   SST: standard synchrosqueezing
%   SST2: vertical second-order synchrosqueezing
%   omega: instantaneous frequency (vertical reassignment operator)
%   omega2: second-order instantaneous frequency
% REFERENCES

s = s(:);
N = length(s);

% Padding
% sz=zeros(N, 1);
% sleft = flipud(conj(sz(2:N/2+1)));
% sright = flipud(sz(end-N/2:end-1));
% x = [sleft; s ; sright];
x = s;
clear xleft xright;

nb = N; 
% neta = Nfft;
neta = max_f;

% Window definition
[g, Lh] = create_gaussian_window(N, Nfft, sigma_w, 10^(-3));
a   = pi/sigma_w^2;
gp  = -2*a*(-Lh:Lh)'/N.*g;


% Initialization
STFT = zeros(neta,nb);
omega = zeros(neta,nb);
omega2 = zeros(neta,nb);
phi2_hat = zeros(neta,nb);
phi2_bar = zeros(neta,nb);

%storing the different STFTs
bt = 1:N;       
ft = 1:neta; 

Vg = zeros(neta,3);
Vgp = zeros(neta,1);

 %computation of the different STFTs
 for b=1:N
    % STFT, window x^i*g
    ti = (-min([Lh,b-1]):min([Lh,N-b]))';
    for i = 0:2
%      y = fft(x(bt(b):bt(b)+N-1).*(t.^i).*g)/N;
     y = fft(x(bt(b)+ti).*((ti/N).^i).*g(Lh + ti + 1), Nfft)/N;
%      Vg((b-1)*neta+1:b*neta,i+1) = y(ft); 
     Vg(:,i+1) = y(ft); 
    end
    
%      y=fft(x(bt(b):bt(b)+N-1).*(t.^i).*gp)/N;
     y=fft(x(bt(b)+ti).*gp(Lh + ti + 1), Nfft)/N;
%      Vgp((b-1)*neta+1:b*neta,1) = y(ft);
     Vgp(:,1) = y(ft);
    
   %computation of the different determinant to compute the different
   %omegas  
   
   %determinant used in omega2
   
   detD2 =  Vg(:,1).*Vg(:,3)-Vg(:,2).*Vg(:,2);
   detU2 = -Vg(:,1).*Vg(:,2);

   detU2sec =  Vg(:,1).*Vg(:,1);
  
   omega(:,b) = (ft(:)-1)-1/(2*pi)*imag(Vgp(:,1)./Vg(:,1));
   omega2(:,b) = (ft(:)-1)-1/(2*pi)*imag(detU2./detD2);
   
   phi2_hat(:,b)        = -1/(2*pi)*imag(detU2sec./detD2);
   phi2_bar(:,b) = -1/(2*pi)*imag(Vg(:,1)./Vg(:,3));
   STFT(:,b) = Vg(:,1).*exp(2i*pi*(ft-1)'*min(Lh,b-1)/Nfft);   
 end

end