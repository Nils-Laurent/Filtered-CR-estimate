function  [STFT,STFT_thresh,SST,SST2,SST3,SST4,omega,omega2,omega3,omega4,phi2sec,phi3sec,phi4sec] = sstn_det_simple_prec_der(s,sigma,Nfft,gamma,cas)

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

prec = 10^(-3);
L    = sigma*N;
l    = floor(L*sqrt(-log(prec)/pi))+1;
g    = amgauss(2*l+1,l+1,L);       

nb = N; 
neta = Nfft/2; 

n   = (0:2*l)'-l;
t0  = n/N;
t0  = t0(:);
a   = pi/sigma^2; 
gp  = -2*a*t0.*g; 


% Initialization
STFT = zeros(neta,nb);
STFT_thresh = zeros(neta,nb);
SST = zeros(neta,nb);
SST2 = zeros(neta,nb);
SST3 = zeros(neta,nb);
SST4 = zeros(neta,nb);
omega = zeros(neta,nb);
omega2 = zeros(neta,nb);
omega3 = zeros(neta,nb);
omega4 = zeros(neta,nb);

phi2sec = zeros(neta,nb);
phi3sec = zeros(neta,nb);
phi4sec = zeros(neta,nb);

%storing the different STFTs
vg  = zeros(neta,7);
vgp = zeros(neta,4);
bt = 1:N;       
ft = 1:neta; 

%if cas == 1,
    
 %% Computation of the different thresholds
 Vg = zeros(neta*N,7);
 Vgp = zeros(neta*N,4);

 for b=1:N
    time_inst = -min([l,b-1]):min([l,N-b]);
    
    % STFT, window x^i*g  
    for i = 0:6       
     y = fft(s(bt(b)+time_inst).*((time_inst'/N).^i).*g(l+time_inst+1),Nfft)/N;
     Vg((b-1)*neta+1:b*neta,i+1) = y(ft); 
    end 
    
    % STFT, window x^i*gp
    for i = 0:3
     y=fft(s(bt(b)+time_inst).*((time_inst'/N).^i).*gp(l+time_inst+1),Nfft)/N;   
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
    time_inst = -min([l,b-1]):min([l,N-b]);
    for i = 0:6  
     y = fft(s(bt(b)+time_inst).*((time_inst'/N).^i).*g(l+time_inst+1),Nfft)/N;
     vg(:,i+1) = y(ft);
%      if (cas == 1)
%       vg(:,i+1) = vg(:,i+1).*(abs(vg(:,i+1)) > 3*Vg_thresh(i+1));
%      end    
    end

    % STFT, window x^i*gp
    for i = 0:3
     y = fft(s(bt(b)+time_inst).*((time_inst'/N).^i).*gp(l+time_inst+1),Nfft)/N;
     vgp(:,i+1) = y(ft);
%      if (cas == 1)
%       vgp(:,i+1) = vgp(:,i+1).*(abs(vgp(:,i+1)) > 3*Vgp_thresh(i+1));
%      end 
    end

   %computation of the different determinant to compute the different
   %omegas  
   
   %determinant used in omega2
   
   detD2 =  vg(:,1).*vg(:,3)-vg(:,2).*vg(:,2);
   detU2 = -vg(:,1).*vg(:,2);
   
   %determinant for the computation of the modulation operator
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
   
   %estimation of second order modulation
   phi2sec(:,b) = -1/(2*pi)*imag(detU2sec./detD2);
   phi3sec(:,b) = -1/(2*pi)*imag(detU3sec./detD3);
   phi4sec(:,b) = -1/(2*pi)*imag(detU4sec./detD4);
   
   STFT(:,b) = vg(:,1).*(exp(2*1i*pi*min(l,b-1)*(ft-1)'/Nfft));   
 end  
 if cas == 1
  %gamma = 3*Vg_thresh(1);
  gamma = 2*Vg_thresh(1); 
 end
 
 STFT_thresh = STFT.*(abs(STFT) > gamma);
 for b=1:nb
    for eta=1:neta
        if abs(STFT(eta,b))> gamma
            k = 1+round(Nfft/N*omega(eta,b));
            if (k >= 1) && (k <= neta)
                % original reassignment
                SST(k,b) = SST(k,b) + STFT(eta,b);
            end
            %reassignment using new omega2
            k = 1+round(Nfft/N*omega2(eta,b));
            if k>=1 && k<=neta
                % second-order reassignment: SST2
                SST2(k,b) = SST2(k,b) + STFT(eta,b);
            end
            %reassignment using new omega3
            k = 1+floor(Nfft/N*omega3(eta,b));
            if k>=1 && k<=neta
             % third-order reassignment: SST3
             SST3(k,b) = SST3(k,b) + STFT(eta,b);
            end
            %reassignment using new omega4
            k = 1+floor(Nfft/N*omega4(eta,b));
            if k>=1 && k<=neta
             % fourth-order reassignment: SST4
             SST4(k,b) = SST4(k,b) + STFT(eta,b);
            end
        end
    end
 end
end