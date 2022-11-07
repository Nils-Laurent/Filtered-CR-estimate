function  [omega,omega_simple,omega2,omega2_simple] = omega_comput(s,s_wn,n,sigma,Nfft)


s = s(:);
s_wn = s_wn(:);
n = n(:);

N = length(s);

% Padding
sz=zeros(N,1);
sleft = flipud(conj(sz(2:N/2+1)));
sright = flipud(sz(end-N/2:end-1));
x = [sleft; s ; sright];
x_wn = [sleft; s_wn ; sright];
x_n = [sleft; n ; sright];


%clear xleft xright;

nb = N; 
neta = Nfft/2; 

% Window definition
 t = -0.5:1/N:0.5-1/N;
 t=t';
 g =  1/sigma*exp(-pi/sigma^2*(t.^2));
 a   = pi/sigma^2;
 gp  = -2*a*t.*g; 
 
 omega         = zeros(neta,nb);
 omega_simple  = zeros(neta,nb);
 omega2         = zeros(neta,nb);
 omega2_simple  = zeros(neta,nb);
 
 %storing the different STFTs
 vg  = zeros(neta,2);
 vg_wn = zeros(neta,2);
 vg_n = zeros(neta,2);
 
 bt = 1:N;       
 ft = 1:neta; 
 
 for b=1:N   
   
    % STFT, window x^i*g  
    for i = 0:2       
     y = fft(x(bt(b):bt(b)+N-1).*(t.^i).*g)/N;
     vg(:,i+1) = y(ft); 
     y = fft(x_wn(bt(b):bt(b)+N-1).*(t.^i).*g)/N;
     vg_wn(:,i+1) = y(ft); 
     y = fft(x_n(bt(b):bt(b)+N-1).*(t.^i).*g)/N;
     vg_n(:,i+1) = y(ft); 
    end
    for i = 0:0
     y = fft(x(bt(b):bt(b)+N-1).*(t.^i).*gp)/N;
     vgp = y(ft);
    end 
    %omega(:,b) = (ft(:)-1)-1/(2*pi)*imag(vgp./vg(:,1));
    omega(:,b) = (ft(:)-1)+1/(sigma^2)*imag(vg(:,2)./vg(:,1));
    omega_simple(:,b) = (ft(:)-1)+1/(sigma^2)*imag(vg_wn(:,2)./vg_wn(:,1))...
                             +1/(sigma^2)*imag(vg_n(:,2)./vg_wn(:,1));
    detD2 =  vg(:,1).*vg(:,3)-vg(:,2).*vg(:,2);
    detU2 = -vg(:,1).*vg(:,2);
    omega2(:,b) = (ft(:)-1)-1/(2*pi)*imag(detU2./detD2);
    omega2_simple(:,b) = (ft(:)-1)+ 1/(2*pi)*imag(vg_wn(:,2)./vg_wn(:,3)+ vg_n(:,2)./vg_wn(:,3));
                        %-1/(sigma^2)*imag(vg_n(:,1).*vg_wn(:,2))./(vg_wn(:,1).^2)...
                                               
 end
end