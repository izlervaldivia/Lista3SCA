clc
clear
ap= [   -0.0375   0.0375   
        0.0468  -0.0647 ]

bp = [0.5   0    
      0  0.625]

cp = [  1   0
         0   1 ]

dp = 0*ones(2,2)
%2.Polos y ceros
sys=ss(ap,bp,cp,dp);
figure; pzmap(sys)
pause
% 3. Calculo de las barreras de estabilidad y desempeño
    w=logspace(-1,3,400);
    a=50;b=250;c=100;
    x1=20*ones(1,a); x2=60*zeros(1,b); x3=-20*ones(1,c);
    xt=[x1 x2 x3];
    x4=[10*ones(1,80) 0*zeros(1,320)];
    % barreras de estabilidad

    semilogx(w, xt,'r')
    ylim([-60,60])
    grid;
    hold on;
    semilogx(w,x4,'b')
    title('Barreras de Estabilidad')
    xlabel('w(rad/s)')
    ylabel('dB')

    hold off;
    pause
% 4. Calculo de los valores singulares de la planta
    w=logspace(-2,2,100); 
    sv=sigma(sys,w);
    vsmax(1)=max(sv(1,:));    vsmin(1)=min(sv(1,:));
    vsmax(2)=max(sv(2,:));    vsmin(2)=min(sv(2,:));
    vsmax     %vector de valores singulares maximos
    vsmin     %vector de valores singulares min
pause  
% 5. Plotear valores singulares
% Valores singulares planta original
w = logspace(-2,3,100);
sv = sigma(ss(ap, bp, cp, dp),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Valores singulares de la planta original')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

% Planta con integrador
[ns nc] = size(bp);                      % ns = number of inputs;  nc = number of controls;   
a = [ ap             bp
      0*ones(nc,ns)    0*ones(nc,nc) ]

b = [ 0*ones(ns,nc)
      eye(nc)      ]

c = [ cp  0*ones(nc,nc) ]

d = 0*ones(nc,nc)
% Valores singulares de planta con integrador
sv = sigma(ss(a, b, c, d),w);
sv = 20*log10(sv);
semilogx(w, sv)
title('Valores singulares de la planta con integrador')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

% 6. LQR
q = c'*c;                                            % State Weighting Matrix
rho = 1e-3;                                         
r = rho*eye(nc);                                      % Control Weigthing Matrix
[~, ~, g] = care(a,b,q,r)               % Compute Control Gain Matrix G
pause

% 7. Observador Filtro Kalman
ll =  inv(cp*inv(-ap)*bp + dp);     % Choose ll and lh to match singular values at all frequencies
lh = -inv(ap)*bp*ll;
l = [lh 
     ll];                           % ll, lh - for low and high frequency loop shaping

pnint = eye(nc);                                    % Process Noise Intensity Matrix
mu = 0.01;                                         
mnint = mu*eye(nc)  ;                                            
[~, ~, g1] = care(a',c',l*l', mnint);                          
h = g1'

% 8. Compensador
ak = [ a-b*g-h*c  0*ones(ns+nc,nc)
       g          0*ones(nc,nc) ]

bk = [ h
       0*ones(nc,nc) ]

ck = [0*ones(nc, ns+nc) eye(nc,nc) ]
% Valores singulares del compensador
sv = sigma(ss(ak, bk, ck, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Compensator Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause


% 9.Sensibilidad S y T
% Open Loop Analysis

al = [ ap                     bp*ck
       0*ones(ns+nc+nc,ns)    ak    ]

bl = [ 0*ones(ns,nc)
       bk ]
    
cl = [ cp  0*ones(nc,ns+nc+nc) ]
sv = sigma(ss(al-bl*cl, bl, -cl, eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Sensitivity Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 

sv = sigma(ss(al-bl*cl, bl, cl, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Complementary Sensitivity Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause     
% 10. Step al sistema
% Closed Loop Analysis
[y,t] = step(ss(al-bl*cl, bl, cl, 0*eye(nc)));
subplot(2,2,1)
plot(t,y(:,1,1))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 response caused by input 1')

subplot(2,2,2)
plot(t,y(:,1,2))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 response caused by input 2')

subplot(2,2,3)
plot(t,y(:,2,1))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 2 response caused by input 1')

subplot(2,2,4)
plot(t,y(:,2,2))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 2 response caused by input 2')

