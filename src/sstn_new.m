function  [STFT,SST,SST2,SST3,SST4,omega,omega2,omega3,omega4,phi22p,phi23p,phi33p] = sstn_new(s,sigma,Nfft,gamma,cas)

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
%[1] D.-H. Pham and S. Meignen, “High-order synchrosqueezing transform for 
%multicomponent signals analysis - with an application to gravitational-wave signal,”
%IEEE Transac tions on Signal Processing, vol. 65, pp. 3168–3178, June 2017.

s = s(:);

N = length(s);
prec = 10^(-3);
L    = sigma*N;
l    = floor(L*sqrt(-log(prec)/pi))+1;
g    = amgauss(2*l+1,l+1,L);       

% Window definition

n   = (0:2*l)'-l;
t0  = n/N;
t0  = t0(:);
a   = pi/sigma^2; 
gp  = -2*a*t0.*g; 

nb = N; 
neta = Nfft; 

% Initialization
STFT = zeros(neta,nb);
SST = zeros(neta,nb);
SST2 = zeros(neta,nb);
SST3 = zeros(neta,nb);
SST4 = zeros(neta,nb);
omega = zeros(neta,nb);
tau2 = zeros(neta,nb);
tau3 = zeros(neta,nb);
tau4 = zeros(neta,nb);
omega2 = zeros(neta,nb);
omega3 = zeros(neta,nb);
omega4 = zeros(neta,nb);
phi22p = zeros(neta,nb);
phi23p = zeros(neta,nb);
phi33p = zeros(neta,nb);
phi24p = zeros(neta,nb);
phi34p = zeros(neta,nb);
phi44p = zeros(neta,nb);
vg = zeros(neta,7);
vgp = zeros(neta,5);
Y = zeros(neta,4,4);

%% Computes STFT and reassignment operators
bt = 1:N;       
ft = 1:Nfft; 
if cas == 1,
    
 %% Computation of the different thresholds
 Vg = zeros(neta*N,9);
 Vgp = zeros(neta*N,9);

 for b=1:N
    time_inst = -min([l,b-1]):min([l,N-b]);
    
    % STFT, window x^i*g  
    for i = 0:8       
     Vg((b-1)*neta+1:b*neta,i+1) = fft(s(bt(b)+time_inst).*((time_inst)'/N).^i.*g(l+time_inst+1),Nfft)/N;
    end

    % STFT, window x^i*gp
    for i = 0:8
     Vgp((b-1)*neta+1:b*neta,i+1) = fft(s(bt(b)+time_inst).*((time_inst)'/N).^i.*gp(l+time_inst+1),Nfft)/N;
    end
 end

 Vg_thresh  = zeros(1,7);
 Vgp_thresh = zeros(1,4);

 for i=1:9,
  Vg_thresh(i) = median(abs(real(Vg(:,i))))/0.6745;
 end

 for i=1:9,
  Vgp_thresh(i) = median(abs(real(Vgp(:,i))))/0.6745;
 end
end


for b=1:N
    time_inst = -min([l,b-1]):min([l,N-b]);
 	
    %% STFT, window x^i*g  
    for i = 0:8       
     vg(:,i+1) = fft(s(bt(b)+time_inst).*((time_inst)'/N).^i.*g(l+time_inst+1),Nfft)/N;
      if (cas == 1)
       vg(:,i+1) = vg(:,i+1).*(abs(vg(:,i+1)) > 3*Vg_thresh(i+1));
      end    
    end
    
    %% STFT, window x^i*gp
    for i = 0:8
     vgp(:,i+1) = fft(s(bt(b)+time_inst).*((time_inst)'/N).^i.*gp(l+time_inst+1),Nfft)/N;
      if (cas == 1)
       vgp(:,i+1) = vgp(:,i+1).*(abs(vgp(:,i+1)) > 3*Vgp_thresh(i+1));
      end 
    end
    
    %% second-order operator tau
    tau2(:,b) = -vg(:,2)./vg(:,1);
    %% third order operator tau
    tau3(:,b) = -vg(:,3)./vg(:,1);
    %% four order operator tau
    tau4(:,b) = -vg(:,4)./vg(:,1);
     
       
    %% Y expressions
    for i = 1:7
        for j = 1:7
            if i>=j
                Y(:,i,j) = vg(:,1).*vg(:,i+1) - vg(:,j).*vg(:,i-j+2);
            end
        end
    end
    
    %% W expressions
    W2 = -1/2/1i/pi*(vg(:,1).^2+vg(:,1).*vgp(:,2)-vg(:,2).*vgp(:,1));
    W3 = -1/2/1i/pi*(2*vg(:,1).*vg(:,2)+vg(:,1).*vgp(:,3)-vg(:,3).*vgp(:,1));
    W4 = -1/2/1i/pi*(2*vg(:,1).*vg(:,3)+2*vg(:,2).^2+vg(:,1).*vgp(:,4) - vg(:,4).*vgp(:,1)+vg(:,2).*vgp(:,3) - vg(:,3).*vgp(:,2));
    
    %% operator omega    
    %% adpted to the case Nfft different from N
    omega(:,b) = N/Nfft*(ft-1)'-real(vgp(:,1)/2/1i/pi./vg(:,1));        
    
    %% operator hat p: estimations of frequency modulation  
    %% SST2   
    phi22p(:,b) = W2(:)./Y(:,2,2);
    
    %% second order reassignment operator
    omega2(:,b) = omega(:,b) + real(phi22p(:,b).*tau2(:,b));
         
    %% SST3
    phi33p(:,b) = (W3(:).*Y(:,2,2)-W2(:).*Y(:,3,3))./(Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3));
    phi23p(:,b) = W2(:)./Y(:,2,2) - phi33p(:,b).*Y(:,3,2)./Y(:,2,2);
    
    omega3(:,b) = omega(:,b) + real(phi23p(:,b).*tau2(:,b))+ real(phi33p(:,b).*tau3(:,b));
       
    %% SST4    
    phi44p(:,b) =((Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3)).*W4-(W3.*Y(:,2,2)-W2.*Y(:,3,3)).*(Y(:,5,4)+Y(:,5,3)-Y(:,5,2))+(W3.*Y(:,3,2)-W2.*Y(:,4,3)).*(Y(:,4,4)+Y(:,4,3)-Y(:,4,2)))...
        ./((Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3)).*(Y(:,6,4)+Y(:,6,3)-Y(:,6,2))-(Y(:,5,3).*Y(:,2,2)-Y(:,4,2).*Y(:,3,3)).*(Y(:,5,4)+Y(:,5,3)-Y(:,5,2))+(Y(:,5,3).*Y(:,3,2)-Y(:,4,2).*Y(:,4,3)).*(Y(:,4,4)+Y(:,4,3)-Y(:,4,2)));
    phi34p(:,b) = (W3(:).*Y(:,2,2)-W2(:).*Y(:,3,3))./(Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3))-phi44p(:,b).*(Y(:,5,3).*Y(:,2,2)-Y(:,4,2).*Y(:,3,3))...
        ./(Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3));
    phi24p(:,b) = W2(:)./Y(:,2,2) - phi34p(:,b).*Y(:,3,2)./Y(:,2,2)...
        - phi44p(:,b).*Y(:,4,2)./Y(:,2,2);

    omega4(:,b) = omega(:,b) + real(phi24p(:,b).*tau2(:,b))+ ...
                    real(phi34p(:,b).*tau3(:,b))+ real(phi44p(:,b).*tau4(:,b));

   %% Storing STFT    
   STFT(:,b) = vg(:,1).*(exp(2*1i*pi*min(l,b-1)*(ft-1)'/Nfft));   
end

%% reassignment step
for b=1:nb
    for eta=1:neta
        if abs(STFT(eta,b))> gamma
            %if (abs(real(STFT(eta,b))) > gamma1) || (abs(imag(STFT(eta,b))) > gamma2)
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
            k = 1+round(Nfft/N*omega3(eta,b));
            if k>=1 && k<=neta
                % third-order reassignment: SST3
                SST3(k,b) = SST3(k,b) + STFT(eta,b);
            end
            %reassignment using new omega4
            k = 1+round(Nfft/N*omega4(eta,b));
            if k>=1 && k<=neta
                % fourth-order reassignment: SST4
                SST4(k,b) = SST4(k,b) + STFT(eta,b);
            end
        end
    end
end
end

