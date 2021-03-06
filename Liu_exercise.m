% 5 standard lines of starting code ************** clear % clear variables
clc % clear command window
echo off
echo on % make sure command echo is on
format compact % command window formatting
close all % close all graphic windows

% ************************************************ %

global count
count=0 % used in ODE integration to count iterations


mu=3.986e14 % gravitational parameter [km^3/sec^2]

a_c= 7178160 % chaser semimajor axis [m]
a_t= 7178160 % target semimajor axis [m]
a_v= 7178160 % virtual semimajor axis [m]

period_c=2*pi*sqrt(a_c^3/mu) % chaser orbital period [s]
period_t=2*pi*sqrt(a_t^3/mu) % target orbital period [s]
period_v=2*pi*sqrt(a_v^3/mu) % virtual orbital period [s]

T_c=period_c % chaser period [s]
T_t=period_t % target period [s]
T_v=period_v % virtual period [s]
e_c=0 % chaser eccentricity [-]
e_t=0 % target eccentricity [-]
e_v=0 % virtual eccentricity [-]
n_c=(2*pi)/T_c % average chaser ang velocity [rad/s]
n_t=(2*pi)/T_t % average target ang velocity [rad/s]

p_c=a_c*(1-e_c^2) % chaser parameter [m]
p_t=a_t*(1-e_t^2) % target parameter [m]
p_v=a_v*(1-e_v^2) % target parameter [m]

i_v=0.2 % inclination of virtual orbit
omega_v=0.1 % RAAN of virtual orbit
Re=6.3781e6 % Earth radius
J2=1.083e-3 % J2 perturbation 
omega_E=7.2921159e-5 % ang velocity of earth around its axis
AerS=pi*4^2 % need to add for both chaser and target
c_D=2.2 % drag coefficient



sidelength_c=1.5 % chasesr cube side length
mc=100 % chaser mass [kg]
mt=200 % target mass [kg]

Ic=[(1/6)*mc*sidelength_c^2, 0, 0; % chaser inertial matrix [kg*m^2]
    0, (1/6)*mc*sidelength_c^2, 0;
    0, 0, (1/6)*mc*sidelength_c^2]
It=[50, 0, 0; % target inertial matrix [kg*m^2]
    0, 70, 0;
    0, 0, 100]


w_c_0=[0;0;0] % initial ang velocity chaser [rad/s]
w_t_0=[0.05;-0.06;0.02] % initial ang velocity target [rad/s]

q_c_0=[0.8;0.3464;0.3464;0.3464] % initial chaser quaternion
q_t_0=[1;0;0;0]% initial target quaternion


Ac_0=(q_c_0(1)^2-transpose(q_c_0(2:4))*q_c_0(2:4))*eye(3)+2*(q_c_0(2:4))*transpose(q_c_0(2:4))-2*q_c_0(1)*[0, -q_c_0(4), q_c_0(3);
                                                                                            q_c_0(4), 0, -q_c_0(2);
                                                                                            -q_c_0(3), q_c_0(2), 0] % chaser quaternion rotation mat
                  
At_0=(q_t_0(1)^2-transpose(q_t_0(2:4))*q_t_0(2:4))*eye(3)+2*(q_t_0(2:4))*transpose(q_t_0(2:4))-2*q_t_0(1)*[0, -q_t_0(4), q_t_0(3);
                                                                                            q_t_0(4), 0, -q_t_0(2);
                                                                                            -q_t_0(3), q_t_0(2), 0] % target quaternion rotation mat
                                                                                        
Act_0=Ac_0*At_0 % error quaternion rotion mat

w_r_0=w_c_0-Act_0*w_t_0 % relative ang velocity 

w_r_0tilde = [0, -w_r_0(3), w_r_0(2);
             w_r_0(3), 0, -w_r_0(1);
              -w_r_0(2), w_r_0(1), 0] 

q_r_0=quatProd(quatRecip(q_t_0),q_c_0) % error quaternion
dq_r_0=0.5*[0, -transpose(w_r_0);
            w_r_0, -w_r_0tilde]*q_r_0 % derror quaternion

dq_r_0vec=[0.05;0.3;-0.1] % error quaternion vector component (override previous error quaternion
        
T_0=[0, -q_r_0(4), q_r_0(3); q_r_0(4), 0, -q_r_0(2);-q_r_0(3), q_r_0(2), 0] + q_r_0(1)*eye(3) % transform from dqr to wr

w_e_0=2*inv(T_0)*dq_r_0vec % relative ang velocity 

w_c_0=w_e_0+Act_0*w_t_0 % chaser ang velocity 

theta_c_0=0 % initial chaser angle [rad]
theta_t_0=0 % initial target angle [rad]
theta_v_0=0 % initial virtual angle [rad]
dtheta_c_0=(n_c*(1+e_c*cos(theta_c_0))^2)/((1-e_c^2)^(3/2)) % initial chaser ang velocity [rad/s]
dtheta_t_0=(n_t*(1+e_t*cos(theta_t_0))^2)/((1-e_t^2)^(3/2)) % initial target ang velocity [rad/s]
r_c_0=p_c/(1+e_c*cos(theta_c_0)) %initial chaser distance from earth CoM [m]
r_t_0=p_t/(1+e_t*cos(theta_t_0)) %initial target distance from earth CoM [m]

rho_0=[0; 0; 0] % initial relative position
drho_0=[0; 0; 0] % initial relative velocity

rho_c_0=[1; 2; 3] % initial chaser relative position
drho_c_0=[0; 0; 0] % initial chaser relative velocity

rho_t_0=[2; 3; 4] % initial target relative position (if at 0;0;0 get imaginary number)
drho_t_0=[0; 0; 0] % initial target relative velocity


beta=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1] % sliding parameters Liu
p=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1] % sliding parameters Liu
q=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1] % sliding parameters Liu
eta=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1] % sliding parameters Liu

rho_d=[2; 0; 0] % desired final relative position
drho_d=[0; 0; 0] % desired final relative velocity

q_d=[1;0;0;0] % desired final relatitve atttitude


tsim_noActu = 20 % simulation no actuation
tstep = 1 % integration time step
options = 0
X0Li=[theta_t_0;w_t_0;q_t_0;w_c_0;q_c_0;rho_c_0;drho_c_0;rho_t_0;drho_t_0]
numb=20
% inits=inits2
% X0Li=[theta_t_0;inits(1:3,numb);q_t_0;inits(10:12,numb);inits(19:22,numb);inits(13:15,numb);inits(16:18,numb);inits(4:6,numb);inits(7:9,numb)]

tsim_Actu = 20 % simulation with actuation


[pulses,out]=thrusterIntegration(X0Li,options,tsim_Actu,e_v,p_v,i_v,omega_v,Re,J2,omega_E,AerS,c_D,e_c,n_c,e_t,n_t,Ic,It,mc,mt,mu,beta,p,q,eta,q_d,rho_d,drho_d) % with actuation

%% Plotting 

%%% no actuation

% [t,x] = ode45(@integrationLi,0:tstep:tsim_noAct,X0Li,options,e_v,p_v,i_v,omega_v,Re,J2,omega_E,AerS,c_D,e_c,n_c,e_t,n_t,Ic,It,mc,mt,mu,beta,p,q,eta,q_d,rho_d,drho_d)
% t=transpose([1:1:tsim_noAct+1]) % no actuation 

%%% actuation

t=transpose([1:1:tsim_Actu+1]) % actuation 
x=out
timeAll=size(t);
w_t=transpose(x(1:timeAll(1),2:4));
q_t=transpose(x(1:timeAll(1),5:8));
w_c=transpose(x(1:timeAll(1),9:11));
q_c=transpose(x(1:timeAll(1),12:15));
rho_c=transpose(x(1:timeAll(1),16:18));
drho_c=transpose(x(1:timeAll(1),19:21));
rho_t=transpose(x(1:timeAll(1),22:24));
drho_t=transpose(x(1:timeAll(1),25:27));
rho=rho_c-rho_t
drho=drho_c-drho_t
q_r=quatProdMat(quatRecipMat(q_t),q_c);

for i=1:1:size(t)

Ar=(q_r(1,i)^2-transpose(q_r(2:4,i))*q_r(2:4,i))*eye(3)+2*(q_r(2:4,i))*transpose(q_r(2:4,i))-2*q_r(1,i)*[0, -q_r(4,i), q_r(3,i); 
                                                                                            q_r(4,i), 0, -q_r(2,i);
                                                                                            -q_r(3,i), q_r(2,i), 0]; 
                                                                                        
w_r=w_c(1:3,i)-Ar*w_t(1:3,i);
w_rOut(1:3,i)=w_r;

end


fontSize=25


figure;
plot(t,x(:,1));
title("theta_v");
xlabel("time [s]")
ylabel("angle [rad]")
legend;

figure;
plot(t,x(:,2:4));
title("w_t");
xlabel("time [s]")
ylabel("angular velocity [rad/s]")
legend;

figure;
plot(t,x(:,5:8));
title("q_t");
xlabel("time [s]")
legend;

figure;
plot(t,q_c);
title("q_c");
xlabel("time [s]")
legend;

figure;
plot(t,w_c);
title("w_c");
xlabel("time [s]")
ylabel("velocity [m/s]")
legend;

figure;
plot(t,q_r);
title("q_r");
xlabel("time [s]")
ylabel("quaternions")
legend("q0","q1","q2","q3");

figure;
plot(t,w_rOut);
title("w_r");
xlabel("time [s]")
ylabel("angular velocity [rad/s]")
legend("w_x","w_y","w_z");

figure;
plot(t,x(:,16:18));
ax=gca;
ax.FontSize=17;
title("\rho_c","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("Displacement [m]","FontSize",fontSize)
legend("X","Y","Z")

figure;
plot(t,x(:,19:21));
title("drho_c");
xlabel("time [s]")
ylabel("velocity [m/s]")
legend;

figure;
plot(t,x(:,22:24));
ax=gca;
ax.FontSize=17;
title("\rho_t","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("Displacement [m]","FontSize",fontSize)
legend("X","Y","Z")

figure;
plot(t,x(:,19:21));
title("drho_c");
xlabel("time [s]")
ylabel("velocity [m/s]")
legend;

figure;
plot(t,rho);
title("rho");
xlabel("time [s]")
ylabel("displacement [m]")
legend("x","y","z");

figure;
plot(t,drho);
title("drho");
xlabel("time [s]")
ylabel("velocity [m/s]")
legend("v_x","v_y","v_z");

figure;
plot(t,vecnorm(rho_c));
title("norm(c)");
xlabel("time [s]")
ylabel("velocity [m/s]")
legend;

figure;
plot(t,vecnorm(rho_t));
title("norm(t)");
xlabel("time [s]")
ylabel("velocity [m/s]")
legend;

figure;
plot(t,x(:,29:31));
ax=gca;
ax.FontSize=17;
title("Commanded Force","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("Thrust [N]","FontSize",fontSize)
legend("F_x", "F_y", "F_z","FontSize",17);

figure;
plot(t,x(:,32:34));
ax=gca;
ax.FontSize=17;
title("Commanded Torque","FontSize",fontSize);
xlabel("Time [s]","FontSize",fontSize)
ylabel("Torque [Nm]","FontSize",fontSize)
legend("T_x", "T_y", "T_z","FontSize",17);

% figure;
% spy(pulses(1:12,100:150));
% ax=gca;
% ax.FontSize=17;
% title("Thruster Pulses","FontSize",fontSize);
% xlabel("0.1-second timesteps","FontSize",fontSize)
% ylabel("Thruster #","FontSize",fontSize)

% figure;
% plot(t,logicMat2);
% title("Thrusters");
% xlabel("time [s]")
% ylabel("Force [m]")
% legend;


