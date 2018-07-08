a=16;
f=120;
Fm=3840;
Tm=1/Fm;
T=1/f;
N=T/Tm;
t=(0:N-1)*Tm;
y=5*sin(2*pi*f*(t));
y1=16*sin(2*pi*f*(t));
sum=y+y;

fd=300
N=1750
td=1/fd   %laiko asies zingsnis
T=N*td   %viso laiko
df=1/T   %dazniu asies zingsnis
(0:2)*td
(0:2)*df
