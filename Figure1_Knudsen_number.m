clc
clear all
syms TL TD T zL zD z arfa t_t t0 Q R r n M v r1 rs k i j Kn zyc k kl T P seita d shuchu fai;%下表面温度，上表面温度，温度变量，上表面位置，下表面位置，深度变量，凝结系数计算公式，水分子停留时间，振动时间，挥发焓吸附热，普适气体常数常数，孔隙半径，孔隙通道数量，水的分子量，分子运动速度，月壤干密度，月壤颗粒密度。

zL=0;%计算深度起点，浅
zD=-0.1;%计算深度终点，深
% r1=1390*(-1*z)^0.056;%月壤密度
R=8.314;%普适气体常数
t0=1e-13;%振动时间
M=18;%水的分子量
Q=50000;%挥发焓


d=0.21*10^-3;%2*r;

fai=0.46
k=1.3806505*10^(-23);%玻尔兹曼常数
seita=2.7*10^-10;%水分子直径
P=1*10^-9

TL=39;%月壤浅层位置温度
TD=39;%月壤深层位置温度
r=d/2
n=((fai)^(2/3))/(pi*r^2);%孔隙数.......√((1-139/310 z^0.056 )^(2/3)/nπ)
for i=1:361
if TL<360
TD=TD+1;
TL=TL+1;


zyc=(k.*TD)/(P*pi*seita^2*sqrt(2));%自由程
Kn=zyc/d;%克努森数
Kn=vpa(Kn);%计算结果




shuchu=[TD n Kn]
dlmwrite('Kn1_D0_0.21_0.46.txt',shuchu,'delimiter',',','-append','newline','pc','precision',10)
end
end

