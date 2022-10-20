clc
clear all
syms TL TD T zL zD z arfa t_t t0 Q R r n M v r1 rs k plunarH20 Psvsurface Psvdepth Jflu shuchu porosity;%下表面温度，上表面温度，温度变量，上表面位置，下表面位置，深度变量，凝结系数计算公式，水分子停留时间，振动时间，挥发焓吸附热，普适气体常数常数，孔隙半径，孔隙通道数量，水的分子量，分子运动速度，月壤干密度，月壤颗粒密度。
zL=0;%计算深度起点，浅
zD=-5;%计算深度终点，深


n=1.0*10^7;%孔隙通道数量


R=8.314;%普适气体常数
t0=1e-13;%振动时间
M=18;%水的分子量
plunar=1*10^-9;%月表环境压强
plunarH20=4.2*10^-14;%月表水分子的压强
Q=40000;%挥发焓
TD=210;
TL=250;
for j=1:1000
    if n<1.1*10^8
        
        
        porosity=0.92;%孔隙率;
        for i=1:21
            if porosity>0.0570
                porosity=porosity/2;
                
                k=(TL-TD)/(zL-zD);
                T=TL+k*z;
                arfa=exp(-0.0788-(9.6928*10^-10)*T^3.9159);%%%凝结系数，拟合得来
                arfa2=int(arfa,z,zD,zL);%凝结系数的第二层积分
                
                t_t=t0*exp(Q/(R*T));%停留时间
                tt2=int(t_t,z,zD,zL);%水分子停留时间的积分
                
                v=(8*R*T/(pi*M))^0.5;%分子运动速度
                v2=int(v,z,zD,zL);%水分子运动速度的积分
                
                r=((porosity)^(1/3))/((n*pi)^0.5);%孔隙半径.......√((1-139/310 z^0.056 )^(2/3)/nπ)
                rr=int(r,z,zD,zL);%孔隙半径积分.......∫_0^(z_0)▒√((1-139/310 z^0.056 )^(2/3)/nπ) dz
               

                
                fc=1+(arfa2*tt2*v2)/(((rr)^2)*(zL-zD));%吸附阻碍系数
                
                fc2=(1+(4*(zL-zD)^2)/(rr)^2)^0.5;%碰撞阻碍系数
                
                fc3=fc*fc2;%总阻碍系数
                fc=vpa(fc);
                fc2=vpa(fc2);
                fc3=vpa(fc3);

                Dkthis=(v2/(zL-zD)*(n^0.5)*r*2)/fc3;
                Dkthis=vpa(Dkthis);
                Psvsurface=(611*exp((-51058/(8.314))*((1/TL)-(1/273.16))));%计算点深度位置的饱和蒸汽压
                if Psvsurface<plunarH20
                    Psvsurface=Psvsurface;
                else
                    Psvsurface=plunarH20;
                end
                
                
                
                Psvdepth=(611*exp((-51058/(8.314))*((1/TD)-(1/273.16))));%计算点深度位置的饱和蒸汽压
                
                if Psvdepth<plunar
                    Psvdepth=Psvdepth;
                else
                    Psvdepth=plunar;
                end
                
                Jflu=(-1*Dkthis*M/R*(((Psvdepth/TD)-(Psvsurface/TL))/(zL-zD)))*365*24*3600;%菲克第一定律的扩散量计算（g）
                
                Jflu=vpa(Jflu);
                shuchu=[porosity n Q TD TL fc fc2 fc3 Psvsurface Psvdepth Jflu];
                dlmwrite('Water_loss_rate.txt',shuchu,'delimiter',',','-append','newline','pc','precision',10)
            end
        end
        n=n+0.01*10^7;
    end
end

