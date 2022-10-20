clc
clear all
syms TL TD T zL zD arfa t_t t0 Q R r n M v r1 rs k i j Kn zyc k kl T P seita d shuchu z0 Dkthis1 Dkthis2 Dkthis3 kongxilv np quzhedu daoshu1 pengzhuangzuai1 pengzhuangzuai2 xifuzuai zuai ratio Dkother Dkothereffective;

zL=0;%计算深度起点，浅
zD=-1;%计算深度终点，深
R=8.314;%普适气体常数
M=18;%水的分子量
t0=1e-13;%振动时间
Q=54000;%使用黄土的首层饱和吸附热值进行计算

r=(0.036*(10^-3))/2;%孔隙半径三维%%%体积中值孔径%%%
kongxilv=0.343484;%(rs-r1)/rs;压汞测出来
n=((kongxilv)^(2/3))/((r^2)*pi);

np=n*pi;
daoshu1=0/((np^0.5)*(kongxilv^(2/3)));
quzhedu=2.049;%压汞测出来
v=(8*R*T/(pi*M))^0.5;%分子运动速度
k=1.3806505*10^(-23);%玻尔兹曼常数
seita=2.7*10^-10;%水分子直径
P=1*10^1
TL=19;%月壤浅层位置温度
TD=19;%月壤深层位置温度
for i=1:361
if TL<360
TD=TD+1;
TL=TL+1;
T=TL;
kl=(TL-TD)/(zL-zD);
arfa=exp(-0.0788-(9.6928*10^-10)*T^3.9159);%%%凝结系数，拟合得来


t_t=t0*exp(Q/(R*T));%停留时间

v=(8*R*T/(pi*M))^0.5;%分子运动速度

xifuzuai=1+(arfa*t_t*v)/(2*((1+(daoshu1)^2)^0.5)*(((kongxilv^(2/3))/np)^0.5));%吸附阻碍系数
pengzhuangzuai1=(1+4*np/(kongxilv^(2/3)))^0.5;

ratio=xifuzuai/pengzhuangzuai1;
zuai=xifuzuai*pengzhuangzuai1;%总阻碍系数
Dkthis11=(v*(n^0.5)*r*2)/(pengzhuangzuai1);

Dkthis2=(v*r*2)/(zuai);

Dkother=2/3*r*v;
Dkothereffective=kongxilv/quzhedu*Dkother;

shuchu=[Q T pengzhuangzuai1 xifuzuai ratio zuai Dkthis11 Dkthis2 Dkother Dkothereffective];
dlmwrite('Loess.txt',shuchu,'delimiter',',','-append','newline','pc','precision',10)
end
end
