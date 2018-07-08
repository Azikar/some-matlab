A1=1;
A2=3;
A3=0.5;
f1=288;
f2=800;
f3=100;
th1=0;
th2=pi;
th3=pi/2;
M=48;
%Filtras=ZD;
T=1;
T1=0.005;
T2=0.1;
fd=8000;
K=4;



% 1:
x = fix(10*rand(1,8));
disp('pirmosios atsitiktines sekos reiksmes')
x1 = x(1:4)
disp('antrosios atsitiktines sekos reiksmes')
x2 = x(5:8)
disp('seku koreliacines funkcijos reiksmes')
korr = xcorr(x1, x2)

%2
z = 1/fd; % laiko ašies žingsnis
N = T1/z; % modeliuojamo signalo reiksmiu skaicius
Td=1/(fd/2);
t = (0:7999)*Td;
%t = (0:N-1)*z;
y1 = A1*sin(2*pi*f1*t + th1)+A2*sin(2*pi*f2*t + th2);
y2 = A1*sin(2*pi*f1*t + th1)+A3*sin(2*pi*f3*t + th3);
figure(13)
subplot(2,1,1)
stem(t(1:N), y1(1:N),'.'); box('off'); axis tight;
set(gca, 'fontsize', 12);
subplot(2,1,2)
stem(t(1:N), y2(1:N),'.'); box('off'); axis tight;
set(gca, 'fontsize', 12);

% 3

% dz = 1/fd; % laiko ašies žingsnis
% Nd = T/dz; % reiksmiu skaicius diskrecioje laiko asyje
% Td = (0:Nd-1)*dz; % diskrecios laiko asies reiksmes
Td=1/(fd/2);
korr = xcorr(y1, y2);
KORR = abs(fft(korr,15999));
df = 1/(15999*Td);
fasis = (0:8000-1)*df;
t = (-7999:7999)*Td;
Nf = 15999; % reiksmiu kiekis spektre
figure(14)
subplot(2,1,1)
plot(t, korr); box('off');
set(gca, 'fontsize', 12);
subplot(2,1,2)
stem(fasis(1:8000), KORR(1:8000),'.'); box('off');
set(gca, 'fontsize', 12);



%4###########################
ats1 = randn(1,8000);
korrats = xcorr(ats1, ats1);
figure(3)
plot(t, korrats); box('off');
set(gca, 'fontsize', 12)


%5#######################

ats = randn(1,8000);
f=100;
etal=sin(2*pi*f*t);
korrats = xcorr(y1, etal);
KORRATS = abs(fft(korrats,15999));
figure(4)
subplot(2,1,1)
plot(t, korrats(1:15999)); box('off');
set(gca, 'fontsize', 12);
subplot(2,1,2)
stem(fasis(1:8000), KORRATS(1:8000),'.'); box('off');
set(gca, 'fontsize', 12);



%6#########################

n = 0:8000-1;
t = n*Td;
N1 = 111; % reiksmiu kiekis tenkantis keturiems daznio f1 periodams
zings = zeros(1,8000);
zings(1:N1) = 1;
sinusas = A1*sin(2*pi*f1*t + th1);
etal = zings.*sinusas; % etaloninis signalas
atsp = zeros(1,8000);
% I T2 sekundziu telpa 800 diskreciu reiksmiu
atsp(800:8000) = etal(1:8000-800+1);
atsp = atsp + randn(1,8000);% atspindzio signalas su triuksmu
figure(5)
subplot(2,1,1)
stem(t(1:278), etal(1:278),'.'); box('off'); axis tight;
set(gca, 'fontsize', 12); title('Etaloninis signalas')
subplot(2,1,2)
stem(t,atsp,'.'); box('off'); axis tight;
set(gca, 'fontsize', 12); title('Atspindzio signalas')


%7############################ 12 individ
korr = xcorr(atsp, etal);
[a,b] = max(korr);
disp('atspindzio pasirodymo laikas')
(b-7999)*Td
figure(6)
t = (-7999:7999)*Td;
plot(korr);
korr
% %8#################################
% r12 = ifft(fft(atsp,8000).*conj(fft(etal,8000)));
% [a,b] = max(r12);
% disp('atspindzio pasirodymo laikas')
% (b-7999)*Td
% figure(7)
% t = (0:7999)*Td;
% plot(t,r12,'.-',t(b),r12(b),'or'); box('off'); axis tight;
% set(gca, 'fontsize', 12);
% 
% %9#####################################
% clc
% disp(' xcorr() Koreliacines funkcijos skaiciavimo trukme')
% tic
% for i = 1:500
%  korr = xcorr(atsp, etal);
% end
% toc
% disp('Greitas algoritmas. Koreliacines funkcijos skaiciavimo trukme')
% tic
% for i = 1:500
%  r12 = ifft(fft(atsp,8000).*conj(fft(etal,8000)));
% end
% toc
% %10##############################
% 
% df = 1;
% fasis = (0:4000-1)*df;
% Y = abs(fft(y1,8000));
% y = y1+y2;
% figure(8)
% subplot(2,1,1)
% stem(t(1:40), y(1:40)); box('off');
% set(gca, 'fontsize', 12);
% subplot(2,1,2)
% stem(fasis(1:4000),Y(1:4000),'.'); box off
% set(gca,'FontSize',12);
% 
% %11##############################
% figure(9)
% b = fir1(48, 550/8000,'low');
% [h, fn, s] = freqz(b, 1, 8000, 8000);
% s.plot = 'both';
% s.xunits = 'hz';
% s.yunits = 'linear';
% freqzplot(h, fn, s);
% 
% 
% %12################################
% %-- vel – velinimo reiksmes
% %-- dazn - dazniu asies reiksmes
% [vel, dazn] = grpdelay(b, 1, 8000, 8000);
% clc
% disp('filtro velinimas')
% vel_max = fix(max(vel))*1.25e-004
% figure(10)
% plot(dazn,vel); box('off'); axis tight;
% set(gca, 'fontsize', 12);
% 
% 
% 
% %13#################################
% 
% figure(11)
% stem(t(1:49),b,'.'); box('off'); axis tight;
% set(gca, 'fontsize', 12);
% 
% 
% %14###############################
% 
% % - Tf – filtro reakcijos i vienetiniimpulsa trukme
% % - dz_ch - filtro daznine charakteristika
% Tf = 0.006125;
% df = 1/Tf;
% Nf = fix(8000/df);
% dz_ch = abs(fft(b,Nf));
% fasis = (0:Nf-1)*df;
% figure(12)
% subplot(2,1,1)
% plot(fasis(1:25), dz_ch(1:25),'*-');
% box('off'); set(gca, 'fontsize', 12);
% faz_ch = unwrap(angle(fft(b,Nf)))*(180/pi);
% subplot(2,1,2)
% plot(fasis(1:25), faz_ch(1:25),'*-');
% box('off'); set(gca, 'fontsize', 12);
% 
% %15#########################################
% 
% vimp = zeros(1, 8000);
% vimp(1) = 1;
% VIMP = abs(fft(vimp,8000));
% df = 1/T;
% fasis = (0:8000-1)*df;
% t = (0:8000-1)*Td;
% figure(13)
% subplot(2,1,1)
% stem(t(1:40), vimp(1:40),'.'); box('off');
% set(gca, 'fontsize', 12);
% subplot(2,1,2)
% plot(fasis, VIMP); box('off'); set(gca, 'fontsize', 12);
% 
% 
% %16#################################
% filt_vimp = filter(b,1, vimp);
% FILT_IMP = abs(fft(filt_vimp, 8000));
% figure(14)
% subplot(2,1,1)
% stem(t(1:50),filt_vimp(1:50),'.-');box('off');
% axis tight;
% set(gca, 'fontsize', 12);
% subplot(2,1,2)
% plot(fasis(1:4000),FILT_IMP(1:4000));box('off');
% axis tight;
% set(gca, 'fontsize', 12);
% 
% 
% %17#########################
% 
% ats = randn(1, 8000);
% ATS = abs(fft(ats,8000));
% df = 1/T;
% fasis = (0:8000-1)*df;
% t = (0:8000-1)*Td;
% figure(15)
% subplot(2,1,1)
% plot(t(1:8000), ats(1:8000)); box('off'); axis tight;
% set(gca, 'fontsize', 12); title('Triuksmo signalas')
% subplot(2,1,2)
% plot(fasis, ATS); box('off'); title('Triuksmo spektras')
% set(gca, 'fontsize', 12)
% 
% 
% 
% %18#################################################
% ats = randn(1, 8000);
% ATS = abs(fft(ats,8000));
% df = 1/T;
% fasis = (0:8000-1)*df;
% t = (0:8000-1)*Td;
% filt_vimp = filter(b,1, ats);
% filt_vimp1 = filter(b,1, ATS);
% 
% figure(16)
% subplot(2,1,1)
% plot(t(1:8000), filt_vimp(1:8000)); box('off'); axis tight;
% set(gca, 'fontsize', 12); title('Triuksmo signalas')
% subplot(2,1,2)
% plot(fasis, filt_vimp1); box('off'); title('Triuksmo spektras')
% set(gca, 'fontsize', 12)
% 
% %19#########################################
% filt_y = filter(b,1, y);
% FILT_Y = abs(fft(filt_y, 8000));
% figure(19)
% subplot(2,1,1)
% plot(t(1:80),filt_y(1:80));box('off'); axis tight;
% set(gca, 'fontsize', 12); title('Filtruotas sinusu sumos signalas')
% subplot(2,1,2)
% stem(fasis(1:4000),FILT_Y(1:4000),'.');box('off'); axis tight;
% set(gca, 'fontsize', 12); title('Filtruoto sinusu sumos signalo spektras')
% 
