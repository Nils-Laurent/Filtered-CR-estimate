function  [STFT,phi2sec,phi2sec_simple,vt2g] = sstn_det_simple(s,sigma,Nfft,gamma,cas)

%% sstn_new : computes the STFT of a signal and different versions of synchrosqueezing
% adding threshold conditions for denominators of modulation operators q
% INPUTS:   
%   s: real or complex signal
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
N = length(s);

% Padding
sz=zeros(N, 1);
sleft = flipud(conj(sz(2:N/2+1)));
sright = flipud(sz(end-N/2:end-1));
x = [sleft; s ; sright];
clear xleft xright;

nb = N; 
neta = Nfft;

% Window definition
 t = -0.5:1/N:0.5-1/N;
 t=t';
 g =  1/sigma*exp(-pi/sigma^2*(t.^2));
 a   = pi/sigma^2;
 gp  = -2*a*t.*g; 

% Initialization
STFT = zeros(neta,nb);
SST = zeros(neta,nb);
SST2 = zeros(neta,nb);
SST3 = zeros(neta,nb);
SST4 = zeros(neta,nb);
omega = zeros(neta,nb);
omega2 = zeros(neta,nb);
omega3 = zeros(neta,nb);
omega4 = zeros(neta,nb);
phi2sec = zeros(neta,nb);
phi2sec_simple = zeros(neta,nb);
vt2g = zeros(neta,nb);
phi3sec = zeros(neta,nb);
phi4sec = zeros(neta,nb);

%storing the different STFTs
vg  = zeros(neta,7);
vgp = zeros(neta,4);
bt = 1:N;       
ft = 1:neta; 
    
 %% Computation of the different thresholds
 Vg = zeros(neta*N,7);
 Vgp = zeros(neta*N,4);

 for b=1:N   
   
    % STFT, window x^i*g  
    for i = 0:6       
     y = fft(x(bt(b):bt(b)+N-1).*(t.^i).*g)/N;
     Vg((b-1)*neta+1:b*neta,i+1) = y(ft); 
    end 
    
    % STFT, window x^i*gp
    for i = 0:3
     y=fft(x(bt(b):bt(b)+N-1).*(t.^i).*gp)/N;   
     Vgp((b-1)*neta+1:b*neta,i+1) = y(ft);
    end
 end

 Vg_thresh  = zeros(1,7);
 Vgp_thresh = zeros(1,4);

 for i=1:7,
  Vg_thresh(i) = median(abs(real(Vg(:,i))))/0.6745;
 end

 for i=1:4,
  Vgp_thresh(i) = median(abs(real(Vgp(:,i))))/0.6745;
 end
%end

 %computation of the different STFTs
 for b=1:N
    % STFT, window x^i*g  
    for i = 0:6  
     y = fft(x(bt(b):bt(b)+N-1).*(t.^i).*g)/N;
     vg(:,i+1) = y(ft);  
    end

   %determinant used to compute phi2sec
   
   detD2 =  vg(:,1).*vg(:,3)-vg(:,2).*vg(:,2);
   detU2sec =  vg(:,1).*vg(:,1);
    
   phi2sec(:,b)        = -1/(2*pi)*imag(detU2sec./detD2);
   phi2sec_simple(:,b) = -1/(2*pi)*imag(vg(:,1)./vg(:,3));
   vt2g(:,b) = vg(:,3);
   STFT(:,b) = vg(:,1).*exp(1i*pi*(ft-1)');   
 end
 
 
end