    clc
clear all
syms TL TD T zL zD z arfa t_t t0 Q R r n M v r1 rs k;%下表面温度，上表面温度，温度变量，上表面位置，下表面位置，深度变量，凝结系数计算公式，水分子停留时间，振动时间，挥发焓吸附热，普适气体常数常数，孔隙半径，孔隙通道数量，水的分子量，分子运动速度，月壤干密度，月壤颗粒密度。


zL=0;%计算深度起点，浅
zD=-0.1;%计算深度终点，深

fai=0.23;%孔隙率
r=(0.21*10^-3)/2;%假设的孔隙半径
n=((fai)^(2/3))/(pi*r^2);%孔隙数.......√((1-139/310 z^0.056 )^(2/3)/nπ)
R=8.314;%普适气体常数
t0=1e-13;%振动时间
M=18;%水的分子量
Q=0;%挥发焓
TL=40;%月壤浅层位置温度
TD=40;%月壤深层位置温度
for i=1:3601
if Q<150000
Q=Q+200;
k=(TL-TD)/(zL-zD);
T=TL+k*z;
arfa=exp(-0.0788-(9.6928*10^-10)*T^3.9159);%%%凝结系数，拟合得来
arfa2=int(arfa,z,zD,zL);%凝结系数的第二层积分

t_t=t0*exp(Q/(R*T));%停留时间
tt2=int(t_t,z,zD,zL);%水分子停留时间的积分



v=(8*R*T/(pi*M))^0.5;%分子运动速度
v2=int(v,z,zD,zL);%水分子运动速度的积分

% r=(((rs-r1)/rs)^(1/3))/(n*pi)^0.5;%孔隙半径.......√((1-139/310 z^0.056 )^(2/3)/nπ)
rr=int(r,z,zD,zL);%孔隙半径积分.......∫_0^(z_0)▒√((1-139/310 z^0.056 )^(2/3)/nπ) dz
%%%%%符号积分法计算


fc=1+(arfa2*tt2*v2/(rr)^2)/(zL-zD);%吸附阻碍系数

fc2=(1+(4*(zL-zD)^2)/(rr)^2)^0.5;%碰撞阻碍系数

fc3=fc*fc2;%总阻碍系数
%计算结果

Q=vpa(Q);
fc=vpa(fc);
fc2=vpa(fc2);
fc3=vpa(fc3);

shuchu=[Q n r fai TD TL fc fc2 fc3]
dlmwrite('40K_2.1_200.txt',shuchu,'delimiter',',','-append','newline','pc','precision',10)
end
end


%%%%%第二区计算


Q=0;%挥发焓
TL=120;%月壤浅层位置温度
TD=120;%月壤深层位置温度
for i=1:3601
if Q<150000
Q=Q+200;
k=(TL-TD)/(zL-zD);
T=TL+k*z;
arfa=exp(-0.0788-(9.6928*10^-10)*T^3.9159);%%%凝结系数，拟合得来
arfa2=int(arfa,z,zD,zL);%凝结系数的第二层积分

t_t=t0*exp(Q/(R*T));%停留时间
tt2=int(t_t,z,zD,zL);%水分子停留时间的积分



v=(8*R*T/(pi*M))^0.5;%分子运动速度
v2=int(v,z,zD,zL);%水分子运动速度的积分

% r=(((rs-r1)/rs)^(1/3))/(n*pi)^0.5;%孔隙半径.......√((1-139/310 z^0.056 )^(2/3)/nπ)
rr=int(r,z,zD,zL);%孔隙半径积分.......∫_0^(z_0)▒√((1-139/310 z^0.056 )^(2/3)/nπ) dz
%%%%%符号积分法计算


fc=1+(arfa2*tt2*v2/(rr)^2)/(zL-zD);%吸附阻碍系数

fc2=(1+(4*(zL-zD)^2)/(rr)^2)^0.5;%碰撞阻碍系数

fc3=fc*fc2;%总阻碍系数
%计算结果

Q=vpa(Q);
fc=vpa(fc);
fc2=vpa(fc2);
fc3=vpa(fc3);

shuchu=[Q n r fai TD TL fc fc2 fc3]
dlmwrite('120K2.1_200.txt',shuchu,'delimiter',',','-append','newline','pc','precision',10)
end
end



%%%%%第三区计算


Q=0;%挥发焓
TL=200;%月壤浅层位置温度
TD=200;%月壤深层位置温度
for i=1:3601
if Q<150000
Q=Q+200;
k=(TL-TD)/(zL-zD);
T=TL+k*z;
arfa=exp(-0.0788-(9.6928*10^-10)*T^3.9159);%%%凝结系数，拟合得来
arfa2=int(arfa,z,zD,zL);%凝结系数的第二层积分
%%%%
t_t=t0*exp(Q/(R*T));%停留时间
tt2=int(t_t,z,zD,zL);%水分子停留时间的积分



v=(8*R*T/(pi*M))^0.5;%分子运动速度
v2=int(v,z,zD,zL);%水分子运动速度的积分

% r=(((rs-r1)/rs)^(1/3))/(n*pi)^0.5;%孔隙半径.......√((1-139/310 z^0.056 )^(2/3)/nπ)
rr=int(r,z,zD,zL);%孔隙半径积分.......∫_0^(z_0)▒√((1-139/310 z^0.056 )^(2/3)/nπ) dz
%%%%%符号积分法计算


fc=1+(arfa2*tt2*v2/(rr)^2)/(zL-zD);%吸附阻碍系数

fc2=(1+(4*(zL-zD)^2)/(rr)^2)^0.5;%碰撞阻碍系数

fc3=fc*fc2;%总阻碍系数
%计算结果

Q=vpa(Q);
fc=vpa(fc);
fc2=vpa(fc2);
fc3=vpa(fc3);

shuchu=[Q n r fai TD TL fc fc2 fc3]
dlmwrite('200K_2.1_200.txt',shuchu,'delimiter',',','-append','newline','pc','precision',10)
end
end
