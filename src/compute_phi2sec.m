function  [STFT,phi2sec,phi2sec_simple,extra_term] = compute_phi2sec(s,n,sigma,Nfft)

%% sstn_new : computes the STFT of a signal and different versions of synchrosqueezing
% adding threshold conditions for denominators of modulation operators q
% INPUTS:   
%   s: real or complex signal
%   n: the added noise
%   sigma: the variance of the Gaussian window
%   Nfft: number of frequency bins
%   gamma: threshold on the STFT for reassignment
% OUTPUTS:   
%   STFT: the short-time Fourier transform
%   phi2sec : chirp rate estimator used in the second order
%             synchrosqueezing transform
%   phi2sec_simple : simplified chirp rate estimate 
%   vt2g     : STFT with filter t^2g
% REFERENCES

s = s(:);
n = n(:);
N = length(s);

% Padding
sz=zeros(N, 1);
sleft = flipud(conj(sz(2:N/2+1)));
sright = flipud(sz(end-N/2:end-1));
x = [sleft; s ; sright];
noise = [sleft; n ; sright];
clear xleft xright;

nb = N; 
neta = Nfft;

% Window definition
t = -0.5:1/N:0.5-1/N;
t = t';
g =  1/sigma*exp(-pi/sigma^2*(t.^2));
a = pi/sigma^2;

% Initialization
STFT           = zeros(neta,nb);
phi2sec        = zeros(neta,nb);
phi2sec_simple = zeros(neta,nb);
extra_term     = zeros(neta,nb);

%storing the different STFTs
vg   = zeros(neta,3);
vg_n = zeros(neta,3);
bt = 1:N;       
ft = 1:neta;    

 %computation of the different STFTs
 for b=1:N
    % STFT, window x^i*g  
    for i = 0:2  
     y = fft(x(bt(b):bt(b)+N-1).*(t.^i).*g)/N;
     z = fft(noise(bt(b):bt(b)+N-1).*(t.^i).*g)/N;
     vg(:,i+1)   = y(ft);            
     vg_n(:,i+1) = z(ft);  
    end

   %determinant used to compute phi2sec
   
   detD2 =  vg(:,1).*vg(:,3)-vg(:,2).*vg(:,2);
   detU2sec =  vg(:,1).*vg(:,1);
    
   phi2sec(:,b)        = -1/(2*pi)*imag(detU2sec./detD2);
   phi2sec_simple(:,b) = -1/(2*pi)*imag(vg(:,1)./vg(:,3));
   extra_term(:,b) = 1/(2*pi)*imag((vg(:,1).*vg_n(:,3))./(vg(:,3).^2)-vg_n(:,1)./vg(:,3));
   STFT(:,b) = vg(:,1).*exp(1i*pi*(ft-1)');   
 end
 
 
end