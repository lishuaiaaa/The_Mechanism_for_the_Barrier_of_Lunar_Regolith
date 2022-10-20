clc
clear all
syms TL TD T zL zD z arfa t_t t0 Q R r n M v r1 rs k fai;

zL=0;%计算深度起点，浅
zD=-0.1;%计算深度终点，深

fai=0.23;%孔隙率
r=(0.21*10^-3)/2;%假设的孔隙半径
n=((fai)^(2/3))/(pi*r^2);%孔隙数.

R=8.314;%普适气体常数
t0=1e-13;%振动时间
M=18;%水的分子量
Q=30000;%挥发焓
TL=19;%月壤浅层位置温度
TD=19;%月壤深层位置温度
for i=1:3601
if TL<360
TD=TD+0.5;
TL=TL+0.5;
k=(TL-TD)/(zL-zD);
T=TL+k*z;
arfa=exp(-0.0788-(9.6928*10^-10)*T^3.9159);%%%凝结系数，拟合得来
arfa2=int(arfa,z,zD,zL);%凝结系数的第二层积分

t_t=t0*exp(Q/(R*T));%停留时间
tt2=int(t_t,z,zD,zL);%水分子停留时间的积分



v=(8*R*T/(pi*M))^0.5;%分子运动速度
v2=int(v,z,zD,zL);%水分子运动速度的积分


rr=int(r,z,zD,zL);%孔隙半径积分.......∫_0^(z_0)▒√((1-139/310 z^0.056 )^(2/3)/nπ) dz


fc=1+(arfa2*tt2*v2/(rr)^2)/(zL-zD);%吸附阻碍系数

fc2=(1+(4*(zL-zD)^2)/(rr)^2)^0.5;%碰撞阻碍系数

fc3=fc*fc2;%总阻碍系数
TD=vpa(TD);
TL=vpa(TL);
fc=vpa(fc);
fc2=vpa(fc2);
fc3=vpa(fc3);

shuchu=[Q n r fai TD TL fc fc2 fc3]
dlmwrite('Q_30000_fai_0.23_D_0.21.txt',shuchu,'delimiter',',','-append','newline','pc','precision',10)
end
end

