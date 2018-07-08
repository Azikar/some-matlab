A1=6;
A2=4;
A3=8;
f1=821;
f2=247;
f3=748;
th1=pi/2;
th2=pi;
th3=0;
M=55;
%Filtras=ZD;
T=0.9912;
T1=0.0884;
T2=0.12;
fd=5831;
K=5;



% 1:
x = fix(10*rand(1,8));
disp('pirmosios atsitiktines sekos reiksmes')
x1 = x(1:4)
disp('antrosios atsitiktines sekos reiksmes')
x2 = x(5:8)
disp('seku koreliacines funkcijos reiksmes')
korr = xcorr(x1, x2)


korel=zeros(1,length(x1)+length(x2)-1)
ilg=size(korel,2);


for j=1:ilg
    arg = (j-length(x1)); 

   if(arg < 0)
       neigiamas = 1;
       limitas = length(x1) + arg;
   else
       neigiamas = 0;
       limitas = length(x1) - arg;
   end
    

   for n=1:limitas
        if (neigiamas==0)
              korel(j) = korel(j) + x1(arg + n) * x2(n);
            
        else
        
           
             korel(j) = korel(j) + x1(n) * x2(n - arg);
        end
       
   end
end
korel
%n=length(x1)+length(x2)-1
%fftshift(ifft(fft(x1,korelilg).*conj(fft(x2,korelilg))))

%2
z = 1/fd; % laiko ašies žingsnis
N = T1/z; % modeliuojamo signalo reiksmiu skaicius
Td=1/(fd/2);
t = (0:fd/2-1)*Td;
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

ik=N;
% 3

% dz = 1/fd; % laiko ašies žingsnis
% Nd = T/dz; % reiksmiu skaicius diskrecioje laiko asyje
% Td = (0:Nd-1)*dz; % diskrecios laiko asies reiksmes
Td=1/(fd/2);
korr = xcorr(y1, y2);
KORR = abs(fft(korr,fd-1));
df = 1/(floor(fd-1)*Td);
fasis = (0:floor(fd/2)-1)*df;

t = (-2914:(fd/2)-1)*Td;
Nf = fd-1; % reiksmiu kiekis spektre
figure(14)
subplot(2,1,1)
plot(t, korr); box('off');
set(gca, 'fontsize', 12);
subplot(2,1,2)
stem(fasis(1:fd/2), KORR(1:fd/2),'.'); box('off');
set(gca, 'fontsize', 12);



%4###########################
ats1 = randn(1,floor(fd/2));
korrats = xcorr(ats1, ats1);
figure(3)
plot(t, korrats); box('off');
set(gca, 'fontsize', 12)


%5#######################

ats = randn(1,floor(fd/2));
korrats = xcorr(y1, ats);
KORRATS = abs(fft(korrats,fd-1));
figure(4)
subplot(2,1,1)
plot(t, korrats); box('off');
set(gca, 'fontsize', 12);
subplot(2,1,2)
stem(fasis(1:fd/2), KORRATS(1:fd/2),'.'); box('off');
set(gca, 'fontsize', 12);



%6#########################

n = 0:fd-1;
t = n*Td;
N1 = (fd/f1)*5; % reiksmiu kiekis tenkantis keturiems daznio f1 periodams
zings = zeros(1,floor(fd));
zings(1:N1) = 1;
sinusas = A1*sin(2*pi*f1*t + th1);
etal = zings.*sinusas; % etaloninis signalas
atsp = zeros(1,floor(fd/2));
% I T2 sekundziu telpa 800 diskreciu reiksmiu
eh=T2*fd
atsp(eh:fd) = etal(1:fd-eh+1);
atsp = atsp + randn(1,floor(fd));% atspindzio signalas su triuksmu
figure(5)
subplot(2,1,1)
stem(t(1:278), etal(1:278),'.'); box('off'); axis tight;
set(gca, 'fontsize', 12); title('Etaloninis signalas')
subplot(2,1,2)
stem(t,atsp,'.'); box('off'); axis tight;
set(gca, 'fontsize', 12); title('Atspindzio signalas')


%7############################ 12 individ
korr = xcorr(atsp, etal);
korel=zeros(1,length(atsp)+length(etal)-1)
ilg=size(korel,2);
for j=1:ilg
    arg = (j-length(atsp)); 

   if(arg < 0)
       neigiamas = 1;
       limitas = length(atsp) + arg;
   else
       neigiamas = 0;

       
       limitas = length(atsp) - arg;
   end
    

   for n=1:limitas
        if (neigiamas==0)
              korel(j) = korel(j) + atsp(arg + n) * etal(n);
            
        else
        
           
             korel(j) = korel(j) + atsp(n) * etal(n - arg);
        end
       
   end
end
disp('korel')
korel

[a,b] = max(korr);

disp('atspindzio pasirodymo laikas')
(b-(fd-1))*Td
figure(6)
% t = (floor(-2914):floor(fd/2-1))*Td;
% plot(t,korr,'.-',t(b),korr(b),'or'); 
% figure(30)
% [c,d]=max(korel);
% t = (floor(-2914):floor(fd/2-1))*Td;

korel(d);
% plot(t,korel,'.-',t(d),korel(d),'or'); 
%8#################################
r12 = ifft(fft(atsp,fd).*conj(fft(etal,fd)));
[a,b] = max(r12);
disp('atspindzio pasirodymo laikas')
(b-(fd-1))*Td
figure(7)
t = (0:fd-1)*Td;
plot(t,r12,'.-',t(b),r12(b),'or'); box('off'); axis tight;
set(gca, 'fontsize', 12);

%9#####################################
clc
disp(' xcorr() Koreliacines funkcijos skaiciavimo trukme')
tic
for i = 1:500
 korr = xcorr(atsp, etal);
end
toc
disp('Greitas algoritmas. Koreliacines funkcijos skaiciavimo trukme')
tic
for i = 1:500
 r12 = ifft(fft(atsp,floor(fd)).*conj(fft(etal,floor(fd))));
end
toc
%10##############################

df = 1;
fasis = (0:floor(fd/4)-1)*df;
Y = abs(fft(y1,floor(fd/2)));
y = y1+y2;
figure(8)
subplot(2,1,1)
stem(t(1:20), y(1:20)); box('off');
set(gca, 'fontsize', 12);
subplot(2,1,2)
stem(fasis(1:floor(fd/4)),Y(1:floor(fd/4)),'.'); box off
set(gca,'FontSize',12);

%11##############################
figure(9)
b = fir1(55, 550/floor(fd/2),'low');
[h, fn, s] = freqz(b, 1, Nf, floor(fd/2));
s.plot = 'both';
s.xunits = 'hz';
s.yunits = 'linear';
freqzplot(h, fn, s);


%12################################
%-- vel – velinimo reiksmes
%-- dazn - dazniu asies reiksmes
[vel, dazn] = grpdelay(b, 1, Nf, floor(fd/2));
clc
disp('filtro velinimas')
vel_max = fix(max(vel))*1.25e-004
figure(10)
plot(dazn,vel); box('off'); axis tight;
set(gca, 'fontsize', 12);



%13#################################

figure(11)
stem(t(1:56),b,'.'); box('off'); axis tight;
set(gca, 'fontsize', 12);


%14###############################

% - Tf – filtro reakcijos i vienetiniimpulsa trukme
% - dz_ch - filtro daznine charakteristika
Tf = 0.006125;
df = 1/Tf;
Nf = fix((fd/2)/df);
dz_ch = abs(fft(b,Nf));
fasis = (0:Nf-1)*df;
figure(12)
subplot(2,1,1)
plot(fasis(1:9), dz_ch(1:9),'*-');
box('off'); set(gca, 'fontsize', 12);
faz_ch = unwrap(angle(fft(b,Nf)))*(180/pi);
subplot(2,1,2)
plot(fasis(1:9), flip(faz_ch(1:9)),'*-');
box('off'); set(gca, 'fontsize', 12);

%15#########################################

vimp = zeros(1, floor(fd/2));
vimp(1) = 1;
VIMP = abs(fft(vimp,fd/2));
df = 1/T;
fasis = (0:fd/2-1)*df;
t = (0:fd/2-1)*Td;
figure(13)
subplot(2,1,1)
stem(t(1:40), vimp(1:40),'.'); box('off');
set(gca, 'fontsize', 12);
subplot(2,1,2)
plot(fasis, VIMP); box('off'); set(gca, 'fontsize', 12);


%16#################################
filt_vimp = filter(b,1, vimp);
FILT_IMP = abs(fft(filt_vimp, fd/2));
figure(14)
subplot(2,1,1)
stem(t(1:50),filt_vimp(1:50),'.-');box('off');
axis tight;
set(gca, 'fontsize', 12);
subplot(2,1,2)
plot(fasis(1:2000),FILT_IMP(1:2000));box('off');
axis tight;
set(gca, 'fontsize', 12);


%17#########################

ats = randn(1, floor(fd/2));
ATS = abs(fft(ats,floor(fd/2)));
df = 1/T;
fasis = (0:floor(fd/2)-1)*df;
t = (0:floor(fd/2)-1)*Td;
figure(15)
subplot(2,1,1)
plot(t(1:floor(fd/2)), ats(1:floor(fd/2))); box('off'); axis tight;
set(gca, 'fontsize', 12); title('Triuksmo signalas')
subplot(2,1,2)
plot(fasis, ATS); box('off'); title('Triuksmo spektras')
set(gca, 'fontsize', 12)



%18#################################################
ats = randn(1, floor(fd/2));
ATS = abs(fft(ats,fd/2));
df = 1/T;
fasis = (0:floor(fd/2)-1)*df;
t = (0:floor(fd/2)-1)*Td;
filt_vimp = filter(b,1, ats);
filt_vimp1 = filter(b,1, ATS);

figure(16)
subplot(2,1,1)
plot(t(1:fd/2), filt_vimp(1:fd/2)); box('off'); axis tight;
set(gca, 'fontsize', 12); title('Triuksmo signalas')
subplot(2,1,2)
plot(fasis, filt_vimp1); box('off'); title('Triuksmo spektras')
set(gca, 'fontsize', 12)

%19#########################################
filt_y = filter(b,1, y);
FILT_Y = abs(fft(filt_y, floor(fd/2)));
figure(19)
subplot(2,1,1)
plot(t(1:ik),filt_y(1:ik));box('off'); axis tight;
set(gca, 'fontsize', 12); title('Filtruotas sinusu sumos signalas')
subplot(2,1,2)
stem(fasis(1:floor(fd/4)),FILT_Y(1:floor(fd/4)),'.');box('off'); axis tight;
set(gca, 'fontsize', 12); title('Filtruoto sinusu sumos signalo spektras')
