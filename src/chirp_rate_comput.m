function  [phi2sec,phi2sec_simple,phi2sec_simple2,val_car] = chirp_rate_comput(s,s_wn,n,sigma,Nfft)


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

phi2sec         = zeros(neta,nb);
phi2sec_simple  = zeros(neta,nb);
phi2sec_simple2 = zeros(neta,nb);
val_car         = zeros(neta,nb);

%storing the different STFTs
vg  = zeros(neta,7);
vgp = zeros(neta,4);
vg_wn = zeros(neta,4);
vgp_wn = zeros(neta,7);
vg_n = zeros(neta,4);
vgp_n = zeros(neta,7);
 
bt = 1:N;       
ft = 1:neta; 
 
 for b=1:N   
   
    % STFT, window x^i*g  
    for i = 0:6       
     y = fft(x(bt(b):bt(b)+N-1).*(t.^i).*g)/N;
     vg(:,i+1) = y(ft); 
     y = fft(x_wn(bt(b):bt(b)+N-1).*(t.^i).*g)/N;
     vg_wn(:,i+1) = y(ft); 
     y = fft(x_n(bt(b):bt(b)+N-1).*(t.^i).*g)/N;
     vg_n(:,i+1) = y(ft); 
    end 
    
    % STFT, window x^i*gp
    for i = 0:3
     y=fft(x(bt(b):bt(b)+N-1).*(t.^i).*gp)/N;   
     vgp(:,i+1) = y(ft);
     y=fft(x_wn(bt(b):bt(b)+N-1).*(t.^i).*gp)/N;   
     vgp_wn(:,i+1) = y(ft);
     y=fft(x_n(bt(b):bt(b)+N-1).*(t.^i).*gp)/N;   
     vgp_n(:,i+1) = y(ft);
    end

    detD2 =  vg(:,1).*vg(:,3)-vg(:,2).*vg(:,2);
    detU2sec =  vg(:,1).*vg(:,1);
     
    phi2sec(:,b)         = -1/(2*pi)*imag(detU2sec./detD2);
    phi2sec_simple(:,b)  = -1/(2*pi)*imag(vg(:,1)./vg(:,3));
    phi2sec_simple2(:,b) =-1/(2*pi)*imag(vg_wn(:,1)./vg_wn(:,3))+...
                           1/(2*pi)*imag(vg_n(:,3).*vg_wn(:,1)./(vg_wn(:,3).^2))...
                           -1/(2*pi)*imag(vg_n(:,1)./vg_wn(:,3));%...
      %+1/(2*pi)*imag(vg_n(:,1).*vg_n(:,3)./(vg_wn(:,3).^2)); 
    val_car(:,b) = vg_n(:,3)./vg_wn(:,3)- vg_n(:,1)./vg_wn(:,1);
 end
end