% 6 standard lines of starting code ************** clear % clear variables
clc % clear command window
echo off
echo on % make sure command echo is on
format compact % command window formatting
close all % close all graphic windows

% ************************************************ %
 
Jc = [1200, 100, -200; 
    100, 2200, 300; 
    -200, 300, 3100] % Inertial Matrix
q_0 = [0.1; -0.2; 0.5] % initial Euler params
q0_0 = sqrt(1-((q_0(1))^2 + (q_0(2))^2+ (q_0(3))^2)) %initial base eueler pram calculated using identity
w_0 = [-0.001; 0.003; 0.005] % initial angular velocity
p=42164
r=42164
mu=3.986e14
h=sqrt(p*mu)
etaDot=2*pi/86164 %h/r^2

w_c = [0; -etaDot; 0] % commanded angular rate
q_f = [0; 0; 0] % commanded eueler params
q0_f = 1 % commanded base reuler param
w_tilde_0 = [0, -w_0(3), w_0(2); 
             w_0(3), 0, -w_0(1);
              -w_0(2), w_0(1), 0]


  
          
omega_n= 0.05
zeta=1
c1=2*omega_n^2
c2=zeta/omega_n

A = inv(c1*eye(3))
B = c2*eye(3)

bodyFrame=rotMat2quat(0,1,0,0,0,-1,-1,0,0)
q0c_0 = bodyFrame(1) % initial commanded attitude find from intitial conditions in ECI and LVLH frame
qc_0 = bodyFrame(2:4) % initial commanded attitude find from intitial conditions in ECI and LVLH frame

qc_tilde_0 = [0, -qc_0(3), qc_0(2);
            qc_0(3), 0, -qc_0(1);
            -qc_0(2), qc_0(1), 0]


qe_0=[q0c_0, transpose(qc_0); -qc_0, q0c_0*eye(3)-qc_tilde_0]*[q0_0;q_0] % initial error for all 4 qauternians 

tsim = 2000
tstep = 5
options = 0
X0=[w_0;q0_0;q_0;q0c_0;qc_0]
[t,x] = ode45(@integrationTrackGeoSat,0:tstep:tsim,X0,options,A,B,Jc,w_c,qe_0(1))

for i =1: 1:size(t)
   w=transpose(x(i,1:3))
   q0=transpose(x(i,4))
   q=transpose(x(i,5:7))
   q0c=transpose(x(i,8))
   qc=transpose(x(i,9:11))
    

    w_tilde = [0, -w(3), w(2);
                 w(3), 0, -w(1);
                  -w(2), w(1), 0]
    we=w-w_c

    qc_tilde = [0, -qc(3), qc(2);
                qc(3), 0, -qc(1);
                -qc(2), qc(1), 0]

    qe=[q0c, transpose(qc); -qc, q0c*eye(3)-qc_tilde]*[q0;q]
    qe_out(i,1:4)=qe
    Tc(i,1:3)=w_tilde*Jc*w-Jc*(A\B)*we-sign(qe_0(1))*Jc*(A\qe(2:4))
    Euler(i,1:3)=quat2euler(q0,q(1),q(2),q(3))
    
    
end 

plot(t,Tc)
legend("Tc1","Tc2","Tc3")
title("Control Torque")
figure
plot(t,x(:,1:3))
legend("omega1","omega2","omega3")
title("Angular Rate")
figure
plot(t,Euler)
legend("theta","phi","psi")
title("Euler Angles")
figure
plot(t,x(:,4:7))
legend("q0","q1","q2","q3")
title("Euler Parameters")
figure
plot(t,x(:,8:11))
legend("q0c","q1c","q2c","q3c")
title("Commanded Euler Parameters")
figure
plot(t,(x(:,1:3)-repmat(transpose(w_c),length(x),1))*360/(2*pi)) %deg
legend("omega1err","omega2err","omega3err")
title("Angular Rate Error (deg)")
figure
plot(t,qe_out)
legend("q0err","q1err","q2err","q3err")
title("Euler Parameters error")
