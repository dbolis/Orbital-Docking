%[w_t_out;rho_t_out;drho_t_out;w_c_out;rho_c_out;drho_c_out;q_c_out]



%plot: rho: x, y, z. drho: dx, dy, dz. qr: q0, q1, q2, q3. wr: w1, w2 ,w3. 

% x=[x1, x2, x3, x4, x5 ,x6, x7, x8, x9, x10...
%    x11, x12, x13, x14, x15 ,x16, x17, x18, x19, x20...
%    x21, x22, x23, x24, x25 ,x26, x28, x29...
%    x31, x32, x33, x34, x35 ,x36, x37, x38, x39, x40...
%    x41, x42, x43, x44, x45 ,x46, x47, x48, x49, x50...
%    x51, x52, x53, x54, x55 ,x56, x57, x58, x59, x60...
%    x61, x62, x63, x64, x65 ,x66, x67, x68, x69, x70...
%    x71, x72, x73, x74, x75 ,x76, x77, x78, x79, x80...
%    x81, x82, x83, x84, x85 ,x86, x87, x89, x90...
%    x91, x92, x93, x94, x95 ,x96, x97, x98, x99, x100]

% x=[x27, x30, x88]
t=transpose([1:1:751]);

X=[]
Y=[]
Z=[]
dx=[]
dy=[]
dz=[]
q0=[]
q1=[]
q2=[]
q3=[]
wx=[]
wy=[]
wz=[]


for i=0:1:96

timeAll=size(t);
% error=transpose(x(1:timeAll(1),12:18));
% derror=transpose(x(1:timeAll(1),19:25));
q_t=transpose(x(1:timeAll(1),5+i*34:8+i*34));
w_c=transpose(x(1:timeAll(1),9+i*34:11+i*34));
w_t=transpose(x(1:timeAll(1),2+i*34:4+i*34));
rho_c=transpose(x(1:timeAll(1),16+i*34:18+i*34));
drho_c=transpose(x(1:timeAll(1),19+i*34:21+i*34));
rho_t=transpose(x(1:timeAll(1),22+i*34:24+i*34));
drho_t=transpose(x(1:timeAll(1),25+i*34:27+i*34));
rho=rho_c-rho_t
drho=drho_c-drho_t
% q_r=error(4:7,:)+q_d.*ones(4,timeAll(1));
% dq_r=derror(4:7,:);
% rho_out=error(1:3,:)+rho_d.*ones(3,timeAll(1));
q_c=transpose(x(1:timeAll(1),12+i*34:15+i*34));
q_r=quatProdMat(quatRecipMat(q_t),q_c);

X=[X; rho(1,:)]
Y=[Y; rho(2,:)]
Z=[Z; rho(3,:)]
dx=[dx; drho(1,:)]
dy=[dy; drho(2,:)]
dz=[dz; drho(3,:)]
q0=[q0; q_r(1,:)]
q1=[q1; q_r(2,:)]
q2=[q2; q_r(3,:)]
q3=[q3; q_r(4,:)]

for j=1:1:size(t)

    Ar=(q_r(1,j)^2-transpose(q_r(2:4,j))*q_r(2:4,j))*eye(3)+2*(q_r(2:4,j))*transpose(q_r(2:4,j))-2*q_r(1,j)*[0, -q_r(4,j), q_r(3,j); %EQ 10 Liu
                                                                                            q_r(4,j), 0, -q_r(2,j);
                                                                                            -q_r(3,j), q_r(2,j), 0]; 
    w_r(1:3,j)=w_c(1:3,j)-Ar*w_t(1:3,j);
    
end
    
wx=[wx; w_r(1,:)]
wy=[wy; w_r(2,:)]
wz=[wz; w_r(3,:)]

end 

fontSize=25
axsize=17


figure;
plot(t,X);
ax=gca;
ax.FontSize=axsize;
title("Relative Distance - X Component","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("Displacement [m]","FontSize",fontSize)


figure;
plot(t,Y);
ax=gca;
ax.FontSize=axsize;
title("Relative Distance - Y Component","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("Displacement [m]","FontSize",fontSize)

figure;
plot(t,Z);
ax=gca;
ax.FontSize=axsize;
title("Relative Distance - Z Component","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("Displacement [m]","FontSize",fontSize)

figure;
plot(t,dx);
ax=gca;
ax.FontSize=axsize;
title("Relative Velocity - V_x Component","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("Velocity [m/s]","FontSize",fontSize)

figure;
plot(t,dy);
ax=gca;
ax.FontSize=axsize;
title("Relative Velocity - V_y Component","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("Velocity [m/s]","FontSize",fontSize)

figure;
plot(t,dz);
ax=gca;
ax.FontSize=axsize;
title("Relative Velocity - V_z Component","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("Velocity [m/s]","FontSize",fontSize)

figure;
plot(t,q0);
ax=gca;
ax.FontSize=axsize;
title("Relative Quaternion - q_0 Component","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("q_0","FontSize",fontSize)

figure;
plot(t,q1);
ax=gca;
ax.FontSize=axsize;
title("Relative Quaternion - q_1 Component","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("q_1","FontSize",fontSize)

figure;
plot(t,q2);
ax=gca;
ax.FontSize=axsize;
title("Relative Quaternion - q_2 Component","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("q_2","FontSize",fontSize)

figure;
plot(t,q3);
ax=gca;
ax.FontSize=axsize;
title("Relative Quaternion - q_3 Component","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("q_3","FontSize",fontSize)

figure;
plot(t,wx);
ax=gca;
ax.FontSize=axsize;
title("Relative Angular Rate - \omega_x Component","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize);
ylabel("Angular Velocity [rad/s]","FontSize",fontSize)

figure;
plot(t,wy);
ax=gca;
ax.FontSize=axsize;
title("Relative Angular Rate - \omega_y Component","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("Angular Velocity [rad/s]","FontSize",fontSize)


figure;
plot(t,wz);
ax=gca;
ax.FontSize=axsize;
title("Relative Angular Rate - \omega_z Component","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("Angular Velocity [rad/s]","FontSize",fontSize)

%%%%%%%%%%%%%%%%%

% figure;
% plot(t,X);
% ax=gca;
% ax.FontSize=axsize;
% title("Relative Distance - X Component","FontSize",fontSize);
% xlabel("Time [s]","FontSize",fontSize)
% ylabel("Displacement [m]","FontSize",fontSize)
% 
% figure;
% plot(t,[Y;Z]);
% ax=gca;
% ax.FontSize=axsize;
% title("Relative Distance - Y and Z Components","FontSize",fontSize);
% xlabel("Time [s]","FontSize",fontSize)
% ylabel("Displacement [m]","FontSize",fontSize)
% 
% figure;
% plot(t,[dx;dy;dz]);
% ax=gca;
% ax.FontSize=axsize;
% title("Relative Velocity - V_x, V_y, and V_z Components","FontSize",fontSize);
% xlabel("Time [s]","FontSize",fontSize)
% ylabel("Velocity [m/s]","FontSize",fontSize)
% 
% figure;
% plot(t,q0);
% ax=gca;
% ax.FontSize=axsize;
% title("Relative Quaternion - q0 Component","FontSize",fontSize);
% xlabel("Time [s]","FontSize",fontSize)
% ylabel("q_0","FontSize",fontSize)
% 
% figure;
% plot(t,[q1;q2;q3]);
% ax=gca;
% ax.FontSize=axsize;
% title("Relative Quaternion - q_1, q_2, and q_3 Components","FontSize",fontSize);
% xlabel("Time [s]","FontSize",fontSize)
% ylabel("q_1, q_2, q_3","FontSize",fontSize)
% 
% figure;
% plot(t,[wx;wy;wz]);
% ax=gca;
% ax.FontSize=axsize;
% title("Relative Angular Rate - \omega_x, \omega_y, and \omega_z Components","FontSize",fontSize);
% xlabel("Time [s]","FontSize",fontSize)
% ylabel("Angular Velocity [rad/s]","FontSize",fontSize)