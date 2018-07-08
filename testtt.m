fileID = fopen('C:\Users\merp\Downloads\Duom_2.dat','r');
formatSpec = '%f';
[A,N] = fscanf(fileID,formatSpec);
kernel = ones(41971, 1);
smoothedSignal = conv(A, kernel, 'same');
Ysp = abs(fft(A));
clc
figure(1)
plot(A,'.')

N
fd=24000  %disk daznis
Td=1/fd    %zingsnis laiko asy
T=N*Td     %signalo trukme?
f=1/T  %dazniu asies zingsnis 
m=fd/f   % reiksmiu sk spektre
fasis = (0:1:m-1)*f;

Ysp = abs(fft(A)/(N/2));
% 
% plot(A,'.')
C=zeros(1,fix(m/2));
for i= 1:size(Ysp,1)
    if Ysp(i)>0.05
        C(i)=Ysp(i); 
    end
end
figure(15)

stem(fasis(1:fix(N/2)),Ysp(1:fix(N/2))); box off
set(gca,'FontSize',12); axis tight;
B=unique(C,'stable');
d=zeros(20986,2);
figure(16)

stem(fasis(1:fix(m/2)),C(1:fix(m/2))); box off
set(gca,'FontSize',12); axis tight;
B=unique(C,'stable');
d=zeros(20986,2);

peek=zeros(1,fix((m/2)/100));
for i=1:fix((m/2)/100)
    a=C(i*100:(i+1)*100);
%     [ind, val]=max(a);
%     peek(i,1)=ind;
%     peek(i,2)=val;
    peek(i)=max(a);
end
peek;
 unique(peek,'stable');
 peek(peek==0)=[];
 peek
dazn=zeros(1,fix((m/2)/100));

for j=1:4
    for i=1:fix(m/2)
    
        if round(C(i),4)==round(peek(j),4)
            dazn(j)=fasis(i);
        
         end
    end
end
 dazn(dazn==0)=[];
 dazn

Td=1/fd;
N=T/Td;
n = 0:N;
t = n*Td;
T=N*Td;
df = 1/T;
f0=dazn(4)





fasis = (0:(fd*2)-1)*df;
%ats= rand(1,41971);


etal=sin(2*pi*f0*t);
 

k=xcorr(A,etal);
figure(19)
max(k)
 
et= abs(fft(etal));
plot(et)
figure(20)
plot(k);
KORRATS = abs(fft(k,fd*2));
figure(21)
stem(fasis(1:(fd*2)-1),KORRATS(1:(fd*2)-1));
% figure(16)
% stem(C(1:41000),(1:41000)'.'); box off
% set(gca,'FontSize',12); axis tight;
% 
% figure(17)
% stem(C(1:41971),(1:41971)'.'); box off
% set(gca,'FontSize',12); axis tight;
% for i=1 :20986
%     skait=0;
%     for j=1:41971
%         if B(i)==C(j)
%             skait=skait+1;
%         end
%          d(i,1)=B(i);
%          d(i,2)=skait;
%     end
% end



%  for g= 1:41971
%      if C(g)== 1.913865419221690e+03
%          C(g)
%      end
%  end


% C = unique(C)
% clc

% figure(17)
% stem(C(1:6),(1:6),'.'); box off
% set(gca,'FontSize',12); axis tight;
%  figure(15);
% % hist(A);
%  size(C,2)
% 8*10^4


