% 5 standard lines of starting code ************** clear % clear variables
clc % clear command window
echo off
echo on % make sure command echo is on
format compact % command window formatting
close all % close all graphic windows

% ************************************************ %

global count w_rout TOUT Logic outt
count=0
w_rout=[]
TOUT=[]
Logic=[]
outt=[]

period=2*pi*sqrt(7178160^3/3.986e14)

T_c=period % chaser period [s]
T_t=period % target period [s]
T_v=period % virtual period [s]
e_c=0 % chaser eccentricity [-]
e_t=0 % target eccentricity [-]
e_v=0 % virtual eccentricity [-]
n_c=(2*pi)/T_c % average chaser ang velocity [rad/s]
n_t=(2*pi)/T_t % average target ang velocity [rad/s]

i_v=0.2 % inclination of virtual orbit
omega_v=0.1 % RAAN of virtual orbit
Re=6.3781e6 % Earth radius
J2=1.083e-3 % J2 perturbation 
omega_E=7.2921159e-5 % ang velocity of earth around its axis
AerS=pi*4^2 % need to add for both chaser and target
c_D=2.2 % drag coefficient



mu=3.986e14 % gravitational parameter [km^3/sec^2]
a_c= 7178160 % chaser semimajor axis [m]
a_t= 7178160 % target semimajor axis [m]
a_v= 7178160 % virtual semimajor axis [m]

p_c=a_c*(1-e_c^2) % chaser parameter [m]
p_t=a_t*(1-e_t^2) % target parameter [m]
p_v=a_v*(1-e_v^2) % target parameter [m]

Ic=[40, -12, 20; % chaser inertial matrix [kg*m^2]
    -12, 80, 10;
    20, 10, 50]
It=[50, 0, 0; % target inertial matrix [kg*m^2]
    0, 70, 0;
    0, 0, 100]
mc=20 % chaser mass [kg]


w_c_0=[0;0;0] % initial ang velocity chaser [rad/s]
w_t_0=[0.05;-0.06;0.02] % initial ang velocity target [rad/s]

q_c_0=[0.8;0.3464;0.3464;0.3464] % initial chaser quaternion
q_t_0=[1;0;0;0]% initial target quaternion


Ac_0=(q_c_0(1)^2-transpose(q_c_0(2:4))*q_c_0(2:4))*eye(3)+2*(q_c_0(2:4))*transpose(q_c_0(2:4))-2*q_c_0(1)*[0, -q_c_0(4), q_c_0(3);
                                                                                            q_c_0(4), 0, -q_c_0(2);
                                                                                            -q_c_0(3), q_c_0(2), 0] % Equation 9 Liu
                       
At_0=(q_t_0(1)^2-transpose(q_t_0(2:4))*q_t_0(2:4))*eye(3)+2*(q_t_0(2:4))*transpose(q_t_0(2:4))-2*q_t_0(1)*[0, -q_t_0(4), q_t_0(3);
                                                                                            q_t_0(4), 0, -q_t_0(2);
                                                                                            -q_t_0(3), q_t_0(2), 0] % Equation 9 Liu 
                                                                                        
Act_0=Ac_0*At_0 % Equation 10 Liu

w_r_0=w_c_0-Act_0*w_t_0 % Equation 11 Liu

w_r_0tilde = [0, -w_r_0(3), w_r_0(2);
             w_r_0(3), 0, -w_r_0(1);
              -w_r_0(2), w_r_0(1), 0] 

q_r_0=quatProd(quatRecip(q_t_0),q_c_0) % Equation 13 Liu
dq_r_0=0.5*[0, -transpose(w_r_0);
            w_r_0, -w_r_0tilde]*q_r_0 % Equation 14 Liu

dq_r_0vec=[0.05;0.3;-0.1]
        
T_0=[0, -q_r_0(4), q_r_0(3); q_r_0(4), 0, -q_r_0(2);-q_r_0(3), q_r_0(2), 0] + q_r_0(1)*eye(3)

w_e_0=2*inv(T_0)*dq_r_0vec

w_c_0=w_e_0+Act_0*w_t_0

theta_c_0=0 % initial chaser angle [rad]
theta_t_0=0 % initial target angle [rad]
theta_v_0=0 % initial virtual angle [rad]
dtheta_c_0=(n_c*(1+e_c*cos(theta_c_0))^2)/((1-e_c^2)^(3/2)) % initial chaser ang velocity [rad/s]
dtheta_t_0=(n_t*(1+e_t*cos(theta_t_0))^2)/((1-e_t^2)^(3/2)) % initial target ang velocity [rad/s]
r_c_0=p_c/(1+e_c*cos(theta_c_0)) %initial chaser distance from earth CoM [m]
r_t_0=p_t/(1+e_t*cos(theta_t_0)) %initial target distance from earth CoM [m]

rho_0=[0; 0; 0] % initial relative position
drho_0=[0; 0; 0] % initial relative velocity

rho_c_0=[-20; 12; -7] % initial chaser relative position
drho_c_0=[0.5; -0.7; 1] % initial chaser relative velocity

rho_t_0=[1; 0; 0] % initial target relative position (if at 0;0;0 get imaginary number)
drho_t_0=[0; 0; 0] % initial target relative velocity


beta=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1] % sliding parameters Liu
p=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1] % sliding parameters Liu
q=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1] % sliding parameters Liu
eta=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1] % sliding parameters Liu

rho_d=[2; 0; 0] % desired final relative position
drho_d=[0; 0; 0] % desired final relative velocity

q_d=[1;0;0;0] % desired final relatitve atttitude

error_0=[rho_0-rho_d;q_r_0-q_d] % initial error [position; attitude]
derror_0=[drho_0-drho_d; dq_r_0] % initial derror [velocity; angular velocity]

tsim = 50
tstep = 1
options = 0
X0=[theta_c_0; dtheta_c_0; theta_t_0; dtheta_t_0; rho_0; drho_0; w_c_0; w_t_0;q_c_0;q_t_0;w_r_0;q_r_0;dq_r_0; error_0; derror_0]
X0min=[theta_t_0; rho_0; drho_0;w_t_0;q_t_0;w_r_0;q_r_0;dq_r_0]
X0Error=[theta_t_0;w_t_0;q_t_0;w_r_0;error_0;derror_0]
X0Li=[theta_t_0;w_t_0;q_t_0;w_c_0;q_c_0;rho_c_0;drho_c_0;rho_t_0;drho_t_0]%;q_r_0;dq_r_0]
% [t,x] = ode45(@integrationLiu,0:tstep:tsim,X0,options,e_c,n_c,e_t,n_t,p_c,p_t,Ic,It,mc,mu,beta,p,q,eta)
% [t,x] = ode45(@integrationLiuMin,0:tstep:tsim,X0min,options,e_c,n_c,e_t,n_t,p_c,p_t,Ic,It,mc,mu,beta,p,q,eta)
% [t,x] = ode45(@integrationLiuError,0:tstep:tsim,X0Error,options,e_c,n_c,e_t,n_t,p_c,p_t,Ic,It,mc,mu,beta,p,q,eta,q_d,rho_d,drho_d)

% for i =1:1:size(t)
%     theta_c=transpose(x(i,1))
%     theta_t=transpose(x(i,3))
%     
%     r_c(i,1)=p_c/(1+e_c*cos(theta_c))
%     r_t(i,1)=p_t/(1+e_t*cos(theta_t))
%     
%     
%     
% end


% figure
% polarplot(x(:,1),r_c)
% hold on
% polarplot(x(:,3),r_t,"g")
% 
% figure
% plot(t,x(:,11:13))
% title("omega_c")
% 
% figure
% plot(t,x(:,14:16))
% title("omega_t")
% 
% figure
% plot(t,x(:,17:20))
% title("quat_c")
% 
% figure
% plot(t,x(:,21:24))
% title("quat_t")
% 
% figure
% plot(t,x(:,25:27))
% title("omega_r")
% 
% figure
% plot(t,x(:,28:31))
% title("quat_r")
% 
% figure
% plot(t,x(:,36:42))
% title("Error")


%% Plotting integrationMin
% 
% 
% for i=1:1:size(t)
%     theta_t=transpose(x(i,1))
%     rho=transpose(x(i,2:4))
%     q_t=transpose(x(i,11:14))
%     q_r=transpose(x(i,18:21))
%     w_r=transpose(x(i,15:17))
%     w_t=transpose(x(i,8:10))
%     
%     r_t(i,1)=p_t/(1+e_t*cos(theta_t))
%     r_cVec= [0; 0; r_t(i,1)]+rho
%     r_c(i,1)= norm(r_cVec)
%     
% 
%     
%     
%     q_c=quatProd(q_t,q_r)
%     
%     Ac=(q_c(1)^2-transpose(q_c(2:4))*q_c(2:4))*eye(3)+2*transpose(q_c(2:4))*q_c(2:4)-2*q_c(1)*[0, -q_c(4), q_c(3);
%                                                                                             q_c(4), 0, -q_c(2);
%                                                                                            -q_c(3), q_c(2), 0]                   
%     At=(q_t(1)^2-transpose(q_t(2:4))*q_t(2:4))*eye(3)+2*transpose(q_t(2:4))*q_t(2:4)-2*q_t(1)*[0, -q_t(4), q_t(3);
%                                                                                             q_t(4), 0, -q_t(2);
%                                                                                            -q_t(3), q_t(2), 0]
%     Act=Ac*At
%     
%     w_c(i,1:3)=w_r+Act*w_t
%     rho_out(i,1:3)=rho
% 
%     q_c_out1(i,1:4)=q_c
%     q_r_out1(i,1:4)=q_r
%     q_t_out1(i,1:4)=q_t
% 
% end
% 
% 
% 
% figure
% plot(t,[r_t,r_c])
% title("r_c, rt")
% legend("r_c","r_t")
% 
% figure
% plot(t,q_c_out1)
% title("q_c")
% 
% figure
% plot(t,w_c)
% title("w_c")
% 
% figure
% plot(t,x(:,8:10))
% title("w_t")
% 
% figure
% plot(t,x(:,11:14))
% title("q_t")
% 
% figure
% plot(t,x(:,15:17))
% title("w_r")
% 
% 
% figure
% plot(t,x(:,18:21))
% title("q_r")
% 
% figure
% plot(t,x(:,22:25))
% title("dq_r")
% 
% figure
% plot(t,rho_out)
% title("rho")






%% Plotting integrationError
% tic
% [t,x] = ode45(@integrationLiuError,0:tstep:tsim,X0Error,options,e_c,n_c,e_t,n_t,p_c,p_t,Ic,It,mc,mu,beta,p,q,eta,q_d,rho_d,drho_d)
% timeAll=size(t)
% error=transpose(x(1:timeAll(1),12:18))
% derror=transpose(x(1:timeAll(1),19:25))
% q_t=transpose(x(1:timeAll(1),5:8))
% w_r=transpose(x(1:timeAll(1),9:11))
% w_t=transpose(x(1:timeAll(1),2:4))
% % q_r=transpose(x(1:timeAll(1),26:29))
% q_r=error(4:7,:)+q_d.*ones(4,timeAll(1))
% dq_r=derror(4:7,:)
% rho_out=error(1:3,:)+rho_d.*ones(3,timeAll(1))
% q_c=quatProdMat(q_t,q_r)
% 
% tic
% timeAll=size(t)
% % error=transpose(zeros(timeAll(1),7))
% % derror=transpose(zeros(timeAll(1),7))
% % q_t=transpose(zeros(timeAll(1),4))
% % w_r=transpose(zeros(timeAll(1),3))
% % w_t=transpose(zeros(timeAll(1),3))
% % q_r_out2=zeros(timeAll(1),4)
% % q_t_out2=zeros(timeAll(1),4)
% % q_c_out2=zeros(timeAll(1),4)
% % dq_r=zeros(timeAll(1),4)
% % rho_out=zeros(timeAll(1),3)
% w_c=ones(3,timeAll(1))
% toc
% 
% 
% for i=1:1:size(t)
% %     
% %    
% %     error=transpose(x(i,12:18))
% %     derror=transpose(x(i,19:25))
% %     q_t=transpose(x(i,5:8))
% %     w_r=transpose(x(i,9:11))
% %     w_t=transpose(x(i,2:4))
% %     
% %     
% %     q_r=error(4:7)+q_d
% %     dq_r(i,1:4)=derror(4:7)
% %     rho_out(i,1:3)=error(1:3)+rho_d
% % 
% % 
% %     q_c=quatProd(q_t,q_r)
% 
%     Ac=(q_c(1,i)^2-transpose(q_c(2:4,i))*q_c(2:4,i))*eye(3)+2*(q_c(2:4,i))*transpose(q_c(2:4,i))-2*q_c(1,i)*[0, -q_c(4,i), q_c(3,i);
%                                                                                             q_c(4,i), 0, -q_c(2,i);
%                                                                                            -q_c(3,i), q_c(2,i), 0]                   
%     At=(q_t(1,i)^2-transpose(q_t(2:4,i))*q_t(2:4,i))*eye(3)+2*(q_t(2:4,i))*transpose(q_t(2:4,i))-2*q_t(1,i)*[0, -q_t(4,i), q_t(3,i);
%                                                                                             q_t(4,i), 0, -q_t(2,i);
%                                                                                            -q_t(3,i), q_t(2,i), 0]
%     Act=Ac*At
% %     
%     w_c(1:3,i)=w_r(:,i)+Act*w_t(:,i)
%     
%     
% %     
% %     
% %     q_r_out2(1:4,i)=q_r
% %     q_t_out2(1:4,i)=q_t
% %     q_c_out2(1:4,i)=q_c
% end
% time2=toc
% 
% 
% 
% 
% 
% figure
% plot(t,x(:,1))
% title("theta")
% legend
% 
% figure
% plot(t,x(:,2:4))
% title("w_t")
% legend
% 
% figure
% plot(t,x(:,5:8))
% title("q_t")
% legend
% 
% figure
% plot(t,x(:,9:11))
% title("w_r")
% legend
% 
% figure
% plot(t,x(:,12:14))
% title("position error")
% legend
% 
% figure
% plot(t,x(:,15:18))
% title("quaternian error")
% legend
% 
% figure
% plot(t,x(:,19:21))
% title("position derivative error")
% legend
% 
% figure
% plot(t,x(:,22:25))
% title("quaternion derivative error")
% legend
% 
% figure
% plot(t,q_c)
% title("q_c")
% legend
% 
% figure
% plot(t,w_c)
% title("w_c")
% legend
% 
% figure
% plot(t,q_r)
% title("q_r")
% legend
% 
% figure
% plot(t,dq_r)
% title("dq_r")
% legend
% 
% figure
% plot(t,rho_out)
% title("rho")
% legend

%% Plotting integrationLi
tic
[t,x] = ode45(@integrationLi,0:tstep:tsim,X0Li,options,e_v,p_v,i_v,omega_v,Re,J2,omega_E,AerS,c_D,e_c,n_c,e_t,n_t,Ic,It,mc,mu,beta,p,q,eta,q_d,rho_d,drho_d)
timeAll=size(t);
% error=transpose(x(1:timeAll(1),12:18));
% derror=transpose(x(1:timeAll(1),19:25));
q_t=transpose(x(1:timeAll(1),5:8));
w_c=transpose(x(1:timeAll(1),9:11));
w_t=transpose(x(1:timeAll(1),2:4));
rho_c=transpose(x(1:timeAll(1),16:18));
drho_c=transpose(x(1:timeAll(1),19:21));
rho_t=transpose(x(1:timeAll(1),22:24));
drho_t=transpose(x(1:timeAll(1),25:27));
rho=rho_c-rho_t
drho=drho_c-drho_t
% q_r=error(4:7,:)+q_d.*ones(4,timeAll(1));
% dq_r=derror(4:7,:);
% rho_out=error(1:3,:)+rho_d.*ones(3,timeAll(1));
q_c=transpose(x(1:timeAll(1),12:15));
% q_r=transpose(x(1:timeAll(1),16:19));

% tic;
% timeAll=size(t);
% error=transpose(zeros(timeAll(1),7));
% derror=transpose(zeros(timeAll(1),7));
% q_t=transpose(zeros(timeAll(1),4));
% w_r=transpose(zeros(timeAll(1),3));
% w_t=transpose(zeros(timeAll(1),3));
% q_r_out2=zeros(timeAll(1),4);
% q_t_out2=zeros(timeAll(1),4);
% q_c_out2=zeros(timeAll(1),4);
% dq_r=zeros(timeAll(1),4);
% rho_out=zeros(timeAll(1),3);
% w_c=ones(3,timeAll(1));
% toc
q_r=quatProdMat(quatRecipMat(q_t),q_c);
logicMat=[]
for i=1:1:size(t)
%     
%    
%     error=transpose(x(i,12:18))
%     derror=transpose(x(i,19:25))
%     q_t=transpose(x(i,5:8))
%     w_r=transpose(x(i,9:11))
%     w_t=transpose(x(i,2:4))
%     
%     
%     q_r=error(4:7)+q_d
%     dq_r(i,1:4)=derror(4:7)
%     rho_out(i,1:3)=error(1:3)+rho_d
% 
% 
%     q_c=quatProd(q_t,q_r)
% 
%     Ac=(q_c(1,i)^2-transpose(q_c(2:4,i))*q_c(2:4,i))*eye(3)+2*(q_c(2:4,i))*transpose(q_c(2:4,i))-2*q_c(1,i)*[0, -q_c(4,i), q_c(3,i);
%                                                                                             q_c(4,i), 0, -q_c(2,i);
%                                                                                            -q_c(3,i), q_c(2,i), 0]                   
%     At=(q_t(1,i)^2-transpose(q_t(2:4,i))*q_t(2:4,i))*eye(3)+2*(q_t(2:4,i))*transpose(q_t(2:4,i))-2*q_t(1,i)*[0, -q_t(4,i), q_t(3,i);
%                                                                                             q_t(4,i), 0, -q_t(2,i);
%                                                                                            -q_t(3,i), q_t(2,i), 0]
%     Ar=(q_r(1,i)^2-transpose(q_r(2:4,i))*q_r(2:4,i))*eye(3)+2*(q_r(2:4,i))*transpose(q_r(2:4,i))-2*q_r(1,i)*[0, -q_r(4,i), q_r(3,i) %EQ 10 Liu
%                                                                                             q_r(4,i), 0, -q_r(2,i)
%                                                                                             -q_r(3,i), q_r(2,i), 0]                                                                                                       
%                                                                                        
% %     Act=Ac*At
% %     Act=Ar
% %     
%     w_r(1:3,i)=w_c(:,i)-Ar*w_t(:,i);
% % %     
%     
%     q_r_out2(1:4,i)=q_r
%     q_t_out2(1:4,i)=q_t
%     q_c_out2(1:4,i)=q_c

    
%     Eulersq_r(1:3,i)=quat2euler(q_r_test(1,i),q_r_test(2,i),q_r_test(3,i),q_r_test(4,i))
%     Eulersq_c(1:3,i)=quat2euler(q_c(1,i),q_c(2,i),q_c(3,i),q_c(4,i))
%     Eulersq_t(1:3,i)=quat2euler(q_t(1,i),q_t(2,i),q_t(3,i),q_t(4,i))

% Relative position
Ch=0.5
rho_obs=rho_d
d_0=1.9
delta_0=2
alpha=0.5
etaPos=5e-4


if norm(rho(1:3,i)-rho_obs)>=delta_0
    k=0
else
    k=0
end

gradUatt = Ch*(rho(1:3,i)-rho_d)/sqrt((norm(rho(1:3,i)-rho_d))^2+1)

gradUrep = 4*k*(rho(1:3,i)-rho_obs)*((norm(rho(1:3,i)-rho_obs))^2-delta_0^2)/((norm(rho(1:3,i)-rho_obs))^2-d_0^2)^3

sRel=drho(1:3,i)+alpha*(gradUatt+gradUrep)
sRelout(1:3,i)=sRel

if abs(sRel(1))>=etaPos
    s1Pos=tanh(sRel(1))
else
    s1Pos=sRel(1)/etaPos
end

if abs(sRel(2))>=etaPos
    s2Pos=tanh(sRel(2))
else
    s2Pos=sRel(2)/etaPos
end

if abs(sRel(3))>=etaPos
    s3Pos=tanh(sRel(3))
else
    s3Pos=sRel(3)/etaPos
end

satPos=[s1Pos; s2Pos; s3Pos]


grad2Uatt=Ch*(((norm(rho(1:3,i)-rho_d))^2+1)*eye(3)-(rho(1:3,i)-rho_d)*transpose(rho(1:3,i)-rho_d))/((norm(rho(1:3,i)-rho_d))^2+1)^(3/2)





gamma1=0.02 ;
gamma2=0.005;
gamma3=0.05;
pdiag=[0.5, 0, 0; 0, 0.5, 0; 0, 0, 0.5];
mdiag=[0.3, 0, 0; 0, 0.3, 0; 0, 0, 0.3];
% spos=drho; % s position control, excluding APF term  Eq 20 Li
% satt=derror(2:4); % s att control, excluding APF term Eq 36 
Fclog=-mc*alpha*grad2Uatt*drho(1:3,i)-(gamma1+mc*gamma3)*satPos-pdiag*satPos;
Fc(1:3,i)=-mc*alpha*grad2Uatt*drho(1:3,i)-(gamma1+mc*gamma3)*satPos-pdiag*satPos;  %Eq 27 Li excluding APF terms

Ar=(q_r(1,i)^2-transpose(q_r(2:4,i))*q_r(2:4,i))*eye(3)+2*(q_r(2:4,i))*transpose(q_r(2:4,i))-2*q_r(1,i)*[0, -q_r(4,i), q_r(3,i); %EQ 10 Liu
                                                                                            q_r(4,i), 0, -q_r(2,i);
                                                                                            -q_r(3,i), q_r(2,i), 0]; 

w_r=w_c(1:3,i)-Ar*w_t(1:3,i);
w_rOut(1:3,i)=w_r;

w_r_tilde = [0, -w_r(3), w_r(2); %EQ 12 Liu
             w_r(3), 0, -w_r(1);
              -w_r(2), w_r(1), 0];
          
G_r=[0, -transpose(w_r);
    w_r, -w_r_tilde];

dq_r=0.5*G_r*q_r(1:4,i);


gamma1=0.02 ;
gamma2=0.005;
gamma3=0.05;
pdiag=[0.5, 0, 0; 0, 0.5, 0; 0, 0, 0.5];
mdiag=[0.3, 0, 0; 0, 0.3, 0; 0, 0, 0.3];

Chrel=0.3;
etaAtt=5e-4
beta=0.5;
gradUquat=Chrel*(q_r(2:4,i)-[0;0;0])/sqrt((norm(q_r(2:4,i)-[0;0;0]))^2+1);

sAtt=dq_r(2:4)+beta*gradUquat
sAttout(1:3,i)=sAtt

if abs(sAtt(1))>=etaAtt
    s1Att=sign(sAtt(1))
else
    s1Att=sAtt(1)/etaAtt
    
end

if abs(sAtt(2))>=etaAtt
    s2Att=sign(sAtt(2))
else
    s2Att=sAtt(2)/etaAtt
end

if abs(sAtt(3))>=etaAtt
    s3Att=sign(sAtt(3))
else
    s3Att=sAtt(3)/etaAtt
end

satAtt=[s1Att; s2Att; s3Att]

grad2Uquat=Chrel*(((norm(q_r(2:4,i)-[0;0;0]))^2+1)*eye(3)-(q_r(2:4,i)-[0;0;0])*transpose(q_r(2:4,i)-[0;0;0]))/((norm(q_r(2:4,i)-[0;0;0]))^2+1)^(3/2);

T=[0, -q_r(4,i), q_r(3,i); q_r(4,i), 0, -q_r(2,i);-q_r(3,i), q_r(2,i), 0] + q_r(1,i)*eye(3);
dT=[0, -dq_r(4), dq_r(3); dq_r(4), 0, -dq_r(2);-dq_r(3), dq_r(2), 0] + dq_r(1)*eye(3);
P=inv(T);
dP=inv(dT);
Istar=transpose(P)*Ic*P;

IcPdq_r=Ic*P*dq_r(2:4);
Pdq_r=P*dq_r(2:4);
Arw_t=Ar*w_t(1:3,i);
twoPdq_r=2*P*dq_r(2:4);
w_t_tilde = [0, -w_t(3,i), w_t(2,i); %EQ 8 Liu
             w_t(3,i), 0, -w_t(1,i);
              -w_t(2,i), w_t(1,i), 0];
dw_t=It\(-w_t_tilde*It*w_t(1:3,i));

Cstar=-Istar/dP*P-2*transpose(P)*[0, -IcPdq_r(3), IcPdq_r(2); IcPdq_r(3), 0, -IcPdq_r(1); -IcPdq_r(2), IcPdq_r(1), 0]*P;
Nstar=transpose(P)*([0, -Pdq_r(3), Pdq_r(2); Pdq_r(3), 0, -Pdq_r(1); -Pdq_r(2), Pdq_r(1), 0]*...
    Ic*Ar*w_t(1:3,i))+transpose(P)*([0, -Arw_t(3), Arw_t(2); Arw_t(3), 0, -Arw_t(1);-Arw_t(2), Arw_t(1), 0]*Ic*P*dq_r(2:4))+...
    0.5*transpose(P)*([0, -Arw_t(3), Arw_t(2); Arw_t(3), 0, -Arw_t(1); -Arw_t(2), Arw_t(1), 0]*Ic*Ar*w_t(1:3,i))-...
    0.5*transpose(P)*Ic*([0, -twoPdq_r(3), twoPdq_r(2); twoPdq_r(3), 0, -twoPdq_r(1);-twoPdq_r(2), twoPdq_r(1), 0]*Ar*w_t(1:3,i)-Ar*dw_t);

Tc=Nstar-Cstar*beta*gradUquat-Istar*beta*grad2Uquat*dq_r(2:4)-gamma2*satAtt-mdiag*satAtt;
Tcprimelog = 2*transpose(T)*Tc
Tcprime(1:3,i) = 2*transpose(T)*Tc;
Tcout(1:3,i)=Tc

logicMat=[logicMat, thrusterLogic(transpose(Fclog),transpose(Tcprimelog),q_c(1:4,i), x(i,1), i_v, omega_v)]
end
time2=toc




figure;
plot(t,x(:,1));
title("theta_v");
xlabel("time [s]")
ylabel("angle [rad]")
legend;


% figure;
% plot(t,x(:,2:4));
% title("w_t");
% xlabel("time [s]")
% ylabel("angular velocity [rad/s]")
% legend;
% 
% figure;
% plot(t,x(:,5:8));
% title("q_t");
% xlabel("time [s]")
% legend;

% figure;
% plot(t,w_rOut);
% title("w_r");
% xlabel("time [s]")
% ylabel("angular velocity [rad/s]")
% legend;

figure;
plot(t,x(:,16:18));
title("rho_c");
xlabel("time [s]")
ylabel("displacement [m]")
legend;

figure;
plot(t,x(:,19:21));
title("drho_c");
xlabel("time [s]")
ylabel("velocity [m/s]")
legend;

figure;
plot(t,x(:,22:24));
title("rho_t");
xlabel("time [s]")
ylabel("displacement [m]")
legend;

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
legend;

figure;
plot(t,drho);
title("drho");
xlabel("time [s]")
ylabel("velocity [m/s]")
legend;

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
plot(t,Fc);
title("Fc");
xlabel("time [s]")
ylabel("thrust [N]")
legend;

% figure;
% plot(t,sRelout);
% title("sRelout");
% legend;
% figure
% plot(t,x(:,19:21))
% title("position derivative error")
% legend
% 
% figure
% plot(t,x(:,22:25))
% title("quaternion derivative error")
% legend
% 
% figure;
% plot(t,q_c);
% title("q_c");
% xlabel("time [s]")
% legend;
% 
% figure;
% plot(t,w_c);
% title("w_c");
% xlabel("time [s]")
% ylabel("velocity [m/s]")
% legend;

% figure
% plot(t,q_r)
% title("q_r")
% legend
% 
figure;
plot(t,q_r);
title("q_r");
xlabel("time [s]")
legend;

figure;
plot(t,Tcprime);
title("Tcprime");
xlabel("time [s]")
ylabel("Torque [Nm]")
legend;
% 
% figure;
% plot(t,sAttout);
% title("sAttout");
% legend;
% figure
% plot(t,Eulersq_r)
% title("Eulersq_r")
% legend
% 
% figure
% plot(t,Eulersq_c)
% title("Eulersq_c")
% legend
% 
% figure
% plot(t,Eulersq_t)
% title("Eulersq_t")
% legend

% 
% figure
% plot(t,dq_r)
% title("dq_r")
% legend

% figure
% plot(t,rho_out)
% title("rho")
% legend