function  [STFT,SST,SST2,SST3,SST4,omega,omega2,omega3,omega4,phi2sec,phi2sec_simple,vt2g,phi3sec,phi4sec] = sstn_det(s,sigma,Nfft,gamma,cas)

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
%   SST3: vertical third-order synchrosqueezing
%   SST4: vertical fourth-order synchrosqueezing
%   omega: instantaneous frequency (vertical reassignment operator)
%   omega2: second-order instantaneous frequency
%   omega3: third-order instantaneous frequency
%   omega4: third-order instantaneous frequency
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
     if (cas == 1)
      vg(:,i+1) = vg(:,i+1).*(abs(vg(:,i+1)) > 3*Vg_thresh(i+1));
     end    
    end

    % STFT, window x^i*gp
    for i = 0:3
     y = fft(x(bt(b):bt(b)+N-1).*(t.^i).*gp)/N;
     vgp(:,i+1) = y(ft);
     if (cas == 1)
      vgp(:,i+1) = vgp(:,i+1).*(abs(vgp(:,i+1)) > 3*Vgp_thresh(i+1));
     end 
    end

   %computation of the different determinant to compute the different
   %omegas  
   
   %determinant used in omega2
   
   detD2 =  vg(:,1).*vg(:,3)-vg(:,2).*vg(:,2);
   detU2 = -vg(:,1).*vg(:,2);

   detU2sec =  vg(:,1).*vg(:,1);
   %determinant used in omega3
     
   A = -vg(:,1).*vg(:,2);
   B = -2*vg(:,2).*vg(:,2);
   C = vg(:,1).*vg(:,4)-2*vg(:,2).*vg(:,3);
   
   D = vg(:,1).*vg(:,3)-vg(:,2).*vg(:,2);
   E = vg(:,1).*vg(:,4)-vg(:,2).*vg(:,3);
   F = vg(:,2).*vg(:,4)-vg(:,3).*vg(:,3);
   
   detU3 = vg(:,5).*A-vg(:,4).*B+vg(:,3).*C;
   detD3 = vg(:,5).*D-vg(:,4).*E+vg(:,3).*F;
   
   detU3sec = vg(:,1).*(vg(:,1).*vg(:,5)-vg(:,3).^2)...
             -2*vg(:,2).*(vg(:,1).*vg(:,4)-vg(:,2).*vg(:,3));
   
   %determinant used in omega4
   A = vg(:,3).*(vg(:,5).*vg(:,7)-vg(:,6).*vg(:,6))...
       -vg(:,4).*(vg(:,4).*vg(:,7)-vg(:,6).*vg(:,5))...
       +vg(:,5).*(vg(:,4).*vg(:,6)-vg(:,5).*vg(:,5));
   
   B = vg(:,2).*(vg(:,5).*vg(:,7)-vg(:,6).*vg(:,6))...
       -vg(:,4).*(vg(:,3).*vg(:,7)-vg(:,6).*vg(:,4))...
       +vg(:,5).*(vg(:,3).*vg(:,6)-vg(:,5).*vg(:,4));
   
   C = vg(:,2).*(vg(:,4).*vg(:,7)-vg(:,6).*vg(:,5))...
       -vg(:,3).*(vg(:,3).*vg(:,7)-vg(:,6).*vg(:,4))...
       +vg(:,5).*(vg(:,3).*vg(:,5)-vg(:,4).*vg(:,4));
 
   D = vg(:,2).*(vg(:,4).*vg(:,6)-vg(:,5).*vg(:,5))...
       -vg(:,3).*(vg(:,3).*vg(:,6)-vg(:,5).*vg(:,4))...
       +vg(:,4).*(vg(:,3).*vg(:,5)-vg(:,4).*vg(:,4));
   
   detD4 = vg(:,1).*A-vg(:,2).*B+vg(:,3).*C-vg(:,4).*D;
   detU4 = -vg(:,1).*B+2*vg(:,2).*C-3*vg(:,3).*D;
   
   AA = vg(:,1).*(vg(:,5).*vg(:,7)-vg(:,6).*vg(:,6))...
       -vg(:,3).*(vg(:,3).*vg(:,7)-vg(:,6).*vg(:,4))...
       +vg(:,4).*(vg(:,3).*vg(:,6)-vg(:,5).*vg(:,4));
   
   BB = vg(:,1).*(vg(:,4).*vg(:,7)-vg(:,6).*vg(:,5))...
       -vg(:,2).*(vg(:,3).*vg(:,7)-vg(:,6).*vg(:,4))...
       +vg(:,4).*(vg(:,3).*vg(:,5)-vg(:,4).*vg(:,4));

   CC = vg(:,1).*(vg(:,4).*vg(:,6)-vg(:,5).*vg(:,5))...
       -vg(:,2).*(vg(:,3).*vg(:,6)-vg(:,5).*vg(:,4))...
       +vg(:,3).*(vg(:,3).*vg(:,5)-vg(:,4).*vg(:,4));
   
   detU4sec = vg(:,1).*AA-2*vg(:,2).*BB+3*vg(:,3).*CC;
  
   omega(:,b) = (ft(:)-1)-1/(2*pi)*imag(vgp(:,1)./vg(:,1));
   omega2(:,b) = (ft(:)-1)-1/(2*pi)*imag(detU2./detD2);
   omega3(:,b) = (ft(:)-1)-1/(2*pi)*imag(detU3./detD3);
   omega4(:,b) = (ft(:)-1)-1/(2*pi)*imag(detU4./detD4);
   
   phi2sec(:,b)        = -1/(2*pi)*imag(detU2sec./detD2);
   phi2sec_simple(:,b) = -1/(2*pi)*imag(vg(:,1)./vg(:,3));
   vt2g(:,b) = vg(:,3);
   phi3sec(:,b) = -1/(2*pi)*imag(detU3sec./detD3);
   phi4sec(:,b) = -1/(2*pi)*imag(detU4sec./detD4);
   STFT(:,b) = vg(:,1).*exp(1i*pi*(ft-1)');   
 end
 
 for b=1:nb
    for eta=1:neta
        if abs(STFT(eta,b))> gamma
            k = 1+round(omega(eta,b));
            if (k >= 1) && (k <= neta)
                % original reassignment
                SST(k,b) = SST(k,b) + STFT(eta,b);
            end
            %reassignment using new omega2
            k = 1+round(omega2(eta,b));
            if k>=1 && k<=neta
                % second-order reassignment: SST2
                SST2(k,b) = SST2(k,b) + STFT(eta,b);
            end
            %reassignment using new omega3
            k = 1+floor(omega3(eta,b));
            if k>=1 && k<=neta
             % third-order reassignment: SST3
             SST3(k,b) = SST3(k,b) + STFT(eta,b);
            end
            %reassignment using new omega4
            k = 1+floor(omega4(eta,b));
            if k>=1 && k<=neta
             % fourth-order reassignment: SST4
             SST4(k,b) = SST4(k,b) + STFT(eta,b);
            end
        end
    end
 end
end