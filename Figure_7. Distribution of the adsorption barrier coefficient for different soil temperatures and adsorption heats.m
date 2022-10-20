clc
clear all
syms TL TD T zL zD z arfa t_t t0 Q R r n M v r1 rs k i j;%下表面温度，上表面温度，温度变量，上表面位置，下表面位置，深度变量，凝结系数计算公式，水分子停留时间，振动时间，挥发焓吸附热，普适气体常数常数，孔隙半径，孔隙通道数量，水的分子量，分子运动速度，月壤干密度，月壤颗粒密度。


zL=0;%计算深度起点，浅
zD=-0.1;%计算深度终点，深

fai=0.23;%孔隙率
r=(0.21*10^-3)/2;%假设的孔隙半径
n=((fai)^(2/3))/(pi*r^2);%孔隙数.......√((1-139/310 z^0.056 )^(2/3)/nπ)
R=8.314;%普适气体常数
t0=1e-13;%振动时间
M=18;%水的分子量
Q=0;%挥发焓
for j=1:310
if Q<150000

TL=19;%月壤浅层位置温度
TD=19;%月壤深层位置温度
for i=1:361
if TL<360
TD=TD+1;
TL=TL+1;
k=(TL-TD)/(zL-zD);
T=TL+k*z;
arfa=exp(-0.0788-(9.6928*10^-10)*T^3.9159);%%%凝结系数，拟合得来
arfa2=int(arfa,z,zD,zL);%凝结系数的第二层积分

t_t=t0*exp(Q/(R*T));%停留时间
tt2=int(t_t,z,zD,zL);%水分子停留时间的积分



v=(8*R*T/(pi*M))^0.5;%分子运动速度
v2=int(v,z,zD,zL);%水分子运动速度的积分


rr=int(r,z,zD,zL);%孔隙半径积分.......∫_0^(z_0)▒√((1-139/310 z^0.056 )^(2/3)/nπ) dz
%%%%%符号积分法计算


fc=1+(arfa2*tt2*v2/(rr)^2)/(zL-zD);%吸附阻碍系数

fc2=(1+(4*(zL-zD)^2)/(rr)^2)^0.5;%碰撞阻碍系数

fc3=fc*fc2;%总阻碍系数
TD=vpa(TD);
TL=vpa(TL);
Q=vpa(Q);
fc=vpa(fc);
fc2=vpa(fc2);
fc3=vpa(fc3);
%计算结果

shuchu=[Q n r fai TD TL fc fc2 fc3]
dlmwrite('QK_2.1_.txt',shuchu,'delimiter',',','-append','newline','pc','precision',10)
end
end
Q=Q+500;
end 
end
