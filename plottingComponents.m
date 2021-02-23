%[w_t_out;rho_t_out;drho_t_out;w_c_out;rho_c_out;drho_c_out;q_c_out]



%plot: rho: x, y, z. drho: dx, dy, dz. qr: q0, q1, q2, q3. wr: w1, w2 ,w3. 

x=[x1, x2, x3, x4, x5 ,x6, x7, x8, x9, x10]
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


for i=0:1:9

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

figure;
plot(t,X);
title("Relative Distance - X Component");
xlabel("time [s]")
ylabel("displacement [m]")


figure;
plot(t,Y);
title("Relative Distance - Y Component");
xlabel("time [s]")
ylabel("displacement [m]")

figure;
plot(t,Z);
title("Relative Distance - Z Component");
xlabel("time [s]")
ylabel("displacement [m]")

figure;
plot(t,dx);
title("Relative Velocity - V_x Component");
xlabel("time [s]")
ylabel("velocity [m/s]")

figure;
plot(t,dy);
title("Relative Velocity - V_y Component");
xlabel("time [s]")
ylabel("velocity [m/s]")

figure;
plot(t,dz);
title("Relative Velocity - V_z Component");
xlabel("time [s]")
ylabel("velocity [m/s]")

figure;
plot(t,q0);
title("Relative Quaternion - q_0 Component");
xlabel("time [s]")
ylabel("q_0")

figure;
plot(t,q1);
title("Relative Quaternion - q_1 Component");
xlabel("time [s]")
ylabel("q_1")

figure;
plot(t,q2);
title("Relative Quaternion - q_2 Component");
xlabel("time [s]")
ylabel("q_2")

figure;
plot(t,q3);
title("Relative Quaternion - q_3 Component");
xlabel("time [s]")
ylabel("q_3")

figure;
plot(t,wx);
title("Relative Angular Rate - omega_x Component");
xlabel("time [s]")
ylabel("angular velocity [rad/s]")

figure;
plot(t,wy);
title("Relative Angular Rate - omega_y Component");
xlabel("time [s]")
ylabel("angular velocity [rad/s]")


figure;
plot(t,wz);
title("Relative Angular Rate - omega_z Component");
xlabel("time [s]")
ylabel("angular velocity [rad/s]")

