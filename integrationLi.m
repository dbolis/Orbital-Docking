function der = integrationLi(t,x,e_v,p_v,i_v,omega_v,Re,J2,omega_E,AerS,c_D,e_c,n_c,e_t,n_t,Ic,It,mc,mu,beta,p,q,eta,q_d,rho_d,drho_d)
theta_v=x(1); % True anamoly 
w_t=x(2:4); % target ang velocity
q_t=x(5:8); % target attitude
w_c=x(9:11);% relative angular velocity
q_c=x(12:15)
rho_c=x(16:18)
drho_c=x(19:21)
rho_t=x(22:24)
drho_t=x(25:27)

rho=rho_c-rho_t
drho=drho_c-drho_t

% error=x(12:18); % error
% derror=x(19:25); % velocity error
% q_r=x(16:19)
% dq_r=x(20:23)
% q_r=error(4:7)+q_d; % relative attitude 
% dq_r=derror(4:7); % relative ang rate
% rho=error(1:3)+rho_d; % relative position
% drho=derror(1:3)+drho_d; % relative velocity


global count w_rout TOUT Logic outt logicMat2
count=count+1



%% Dynamics with Perterbances
dtheta_v=sqrt(mu/p_v^3)*(1+e_v*cos(theta_v))^2
d2theta_v=-2*(mu/p_v^3)*(e_v*sin(theta_v))*(1+e_v*cos(theta_v))^3

%%% Chaser

r_v=p_v/(1+e_v*cos(theta_v)) % virtual radius from CoM
r_cVec_v=[r_v;0;0]+rho_c % chaser radius vector in virtual frame
v_cVec_v=[sqrt(mu/p_v)*e_v*sin(theta_v)+drho_c(1)-dtheta_v*rho_c(2); sqrt(mu/p_v)*(1+e_v*cos(theta_v))+drho_c(2)+dtheta_v*rho_c(1); drho_c(3)] % chaser velocity vector in virtual frame

r_cVec_ECI=transpose(r_cVec_v)*ECI2LVLH313(theta_v,i_v,omega_v) % chaser radius vector in ECI
v_cVec_ECI=transpose(v_cVec_v)*ECI2LVLH313(theta_v,i_v,omega_v) % chaser velocity vector in ECI


h_cVec=cross(r_cVec_ECI,v_cVec_ECI) % specific angular momentum vector
h_cVers=h_cVec/norm(h_cVec) % specific angular momentum versor

r_cVers=r_cVec_ECI/norm(r_cVec_ECI) % chaser radius versor in ECI
theta_cVers=cross(h_cVers,r_cVers) % theta versor

phi_c=asin(r_cVers(3)) % latitude
clambda_c=r_cVers(1)/cos(phi_c) 
slambda_c=r_cVers(2)/cos(phi_c) 

lambda_c=2*atan(slambda_c/(1+clambda_c)) % absolute longitutude

szeta_c=theta_cVers(3)/cos(phi_c) 
czeta_c=h_cVers(3)/cos(phi_c)

zeta_c=2*atan(szeta_c/(1+czeta_c)) % heading angle

%%% J2
fJ2_c_c=3*Re^2*J2*mu/(norm(r_cVec_ECI))^4*[(3*sin(phi_c)^2-1)/2, -sin(phi_c)*cos(phi_c)*sin(zeta_c), -sin(phi_c)*cos(phi_c)*cos(zeta_c)] % J2 disturbance on chaser in chaser frame

fJ2_c_v=fJ2_c_c*ECI2LVLH123(zeta_c, phi_c, lambda_c)*transpose(ECI2LVLH313(theta_v,i_v,omega_v)) % J2 disturbance on chaser in virtual frame

%%% Drag

a_c=mu/((2*mu/norm(r_cVec_ECI))-norm(v_cVec_ECI)^2)
e_c=sqrt(1-(norm(h_cVec)^2/(mu*a_c)))
v_radial_c=dot(v_cVec_ECI,r_cVec_ECI)/norm(r_cVec_ECI)
p_c=a_c*(1-e_c^2)

cThetaStar_c=1/e_c*(p_c/norm(r_cVec_ECI)-1)
sThetaStar_c=v_radial_c/e_c*sqrt(p_c/mu)

theta_c=2*atan(sThetaStar_c/(1+cThetaStar_c))


v_r_c=sqrt(mu/p_c)*e_c*sin(theta_c)
v_e_c=sqrt(mu/p_c)*(1+e_c*cos(theta_c))*cos(zeta_c)
v_n_c=sqrt(mu/p_c)*(1+e_c*cos(theta_c))*sin(zeta_c)

v_r_cVec=[v_r_c, v_e_c-omega_E*norm(r_cVec_ECI)*cos(phi_c) v_n_c]*transpose([1 0, 0; 0 cos(zeta_c) -sin(zeta_c); 0 sin(zeta_c) cos(zeta_c)])

atmDensity_c=1.225*exp(-(norm(r_cVec_ECI)-Re)/10332.6) %% apply better fit.. negligible at 800 km altitude

fDrag_c_c=-0.5*c_D*AerS/mc*atmDensity_c*norm(v_r_cVec)*v_r_cVec % in LVLH chaser frame
fDrag_c_v=fDrag_c_c*ECI2LVLH123(zeta_c, phi_c, lambda_c)*transpose(ECI2LVLH313(theta_v,i_v,omega_v)) % in LVLH virtual frame

%%% Solar Radiation Pressure
theta_sol_0=0 % theta_sol at t=0
theta_sol=theta_sol_0+(2*pi)*(t-0)/31558149.5 % sun position in ECI frame
ecl_obliq=deg2rad(23.45) % earth obliquity 

i_c=acos(h_cVers(3))
somega_c=h_cVers(1)/sin(i_c)
comega_c=-h_cVers(2)/sin(i_c)
omega_c=2*atan(somega_c/(1+comega_c))

r_sunVec_c=1.495978e11*[cos(theta_sol), cos(ecl_obliq)*sin(theta_sol), sin(ecl_obliq)*sin(theta_sol)]*transpose(ECI2LVLH313(theta_c,i_c,omega_c))

nu=1 % shadow function 
P_sr=4.56e-6 % solar radiation pressure
c_R_c=1.5 % radiation pressure coefficient
fSolar_c_c=-nu*P_sr*AerS*c_R_c/mc*r_sunVec_c/norm(r_sunVec_c)
fSolar_c_v=fSolar_c_c*ECI2LVLH123(zeta_c, phi_c, lambda_c)*transpose(ECI2LVLH313(theta_v,i_v,omega_v))

%%% target


r_tVec_v=[r_v;0;0]+rho_t  % target radius vector in virtual frame
v_tVec_v=[sqrt(mu/p_v)*e_v*sin(theta_v)+drho_t(1)-dtheta_v*rho_t(2); sqrt(mu/p_v)*(1+e_v*cos(theta_v))+drho_t(2)+dtheta_v*rho_t(1); drho_t(3)] % target velocity vector in virtual frame

r_tVec_ECI=transpose(r_tVec_v)*ECI2LVLH313(theta_v,i_v,omega_v) % target radius vector in ECI
v_tVec_ECI=transpose(v_tVec_v)*ECI2LVLH313(theta_v,i_v,omega_v) % target velocity vector in ECI


h_tVec=cross(r_tVec_ECI,v_tVec_ECI) % specific angular momentum vector
h_tVers=h_tVec/norm(h_tVec) % specific angular momentum versor

r_tVers=r_tVec_ECI/norm(r_tVec_ECI) % target radius versor in ECI
theta_tVers=cross(h_tVers,r_tVers) % theta versor 

phi_t=asin(r_tVers(3)) % latitude
clambda_t=r_tVers(1)/cos(phi_t)
slambda_t=r_tVers(2)/cos(phi_t)

lambda_t=2*atan(slambda_t/(1+clambda_t)) % absolute longitutude

szeta_t=theta_tVers(3)/cos(phi_t)
czeta_t=h_tVers(3)/cos(phi_t)

zeta_t=2*atan(szeta_t/(1+czeta_t)) % heading angle

%%J2
fJ2_t_t=3*Re^2*J2*mu/(norm(r_tVec_ECI))^4*[(3*sin(phi_t)^2-1)/2, -sin(phi_t)*cos(phi_t)*sin(zeta_t), -sin(phi_t)*cos(phi_t)*cos(zeta_t)] % J2 disturbance on target in target frame

fJ2_t_v=fJ2_t_t*ECI2LVLH123(zeta_t, phi_t, lambda_t)*transpose(ECI2LVLH313(theta_v,i_v,omega_v)) % J2 disturbance on target in virtual frame

%%% Drag
a_t=mu/((2*mu/norm(r_tVec_ECI))-norm(v_tVec_ECI)^2)
e_t=sqrt(1-(norm(h_tVec)^2/(mu*a_t)))
v_radial_t=dot(v_tVec_ECI,r_tVec_ECI)/norm(r_tVec_ECI)
p_t=a_t*(1-e_t^2)

cThetaStar_t=1/e_t*(p_t/norm(r_tVec_ECI)-1)
sThetaStar_t=v_radial_t/e_t*sqrt(p_t/mu)

theta_t=2*atan(sThetaStar_t/(1+cThetaStar_t))

v_r_t=sqrt(mu/p_t)*e_t*sin(theta_t)
v_e_t=sqrt(mu/p_t)*(1+e_t*cos(theta_t))*cos(zeta_t)
v_n_t=sqrt(mu/p_t)*(1+e_t*cos(theta_t))*sin(zeta_t)

v_r_tVec=[v_r_t, v_e_t-omega_E*norm(r_tVec_ECI)*cos(phi_t) v_n_t]*transpose([1 0, 0; 0 cos(zeta_t) -sin(zeta_t); 0 sin(zeta_t) cos(zeta_t)])

atmDensity_t=1.225*exp(-(norm(r_tVec_ECI)-Re)/10332.6)

fDrag_t_t=-0.5*c_D*AerS/mc*atmDensity_t*norm(v_r_tVec)*v_r_tVec
fDrag_t_v=fDrag_t_t*ECI2LVLH123(zeta_t, phi_t, lambda_t)*transpose(ECI2LVLH313(theta_v,i_v,omega_v))
%outt=[outt, [fDrag_t_t;fDrag_t_v;v_r_tVec]] 
% outt=[outt, e_t]

%%% Solar Radiation Pressure
i_t=acos(h_tVers(3))
somega_t=h_tVers(1)/sin(i_t)
comega_t=-h_tVers(2)/sin(i_t)
omega_t=2*atan(somega_t/(1+comega_t))

r_sunVec_t=1.495978e11*[cos(theta_sol), cos(ecl_obliq)*sin(theta_sol), sin(ecl_obliq)*sin(theta_sol)]*transpose(ECI2LVLH313(theta_t,i_t,omega_t))

% nu=1 % shadow function 
% P_sr=4.56e-6 % solar radiation pressure
c_R_t=1.5 % radiation pressure coefficient
fSolar_t_t=-nu*P_sr*AerS*c_R_t/mc*r_sunVec_t/norm(r_sunVec_t) 
fSolar_t_v=fSolar_t_t*ECI2LVLH123(zeta_t, phi_t, lambda_t)*transpose(ECI2LVLH313(theta_v,i_v,omega_v))

d_c=transpose(fJ2_c_v)+transpose(fDrag_c_v)+transpose(fSolar_c_v)
d_t=transpose(fJ2_t_v)+transpose(fDrag_t_v)+transpose(fSolar_t_v)

%% Relative Orbit Dynamics

% dtheta_t=(n_t*(1+e_t*cos(theta_t))^2)/((1-e_t^2)^(3/2)); % target orbital velocity
% d2theta_t=(-2*n_t^2*e_t*sin(theta_t)*(1+e_t*cos(theta_t))^3)/((1-e_t^2)^3); % target orbital acceleration

% dtheta_t=sqrt(mu/p_t^3)*(1+e_t*cos(theta_t))^2
% d2theta_t=-2*(mu/p_t^3)*(e_t*sin(theta_t))*(1+e_t*cos(theta_t))^3
% 
% 
% r_t=p_t/(1+e_t*cos(theta_t)); % radius from Earth CoM to target
% r_cVec= [r_t;0;0]+rho; % radius from Earth CoM to chaser VECTOR Pontani coordinates
% % r_cVec= [0;0;r_t]+rho % radius from Earth CoM to chaser VECTOR Liu coordinates
% r_c= norm(r_cVec); % % radius from Earth CoM to chaser scalar



%% Liu coordinates
% A1=[dtheta_t^2, 0, d2theta_t;  % Matrices EQ 4 Liu
%     0, 0, 0;
%     -d2theta_t, 0, dtheta_t^2]
% A2= [0, 0, 2*dtheta_t;
%     0, 0, 0;
%     -2*dtheta_t, 0, 0]
%    
% A3= [-mu*rho(1)/r_c^3;
%     -mu*rho(2)/r_c^3;
%     -mu/r_t^2-mu*(rho(3)-r_t)/r_c^3]


%% Pontani Coordinates

%%% chaser 
X_c=rho_c(1);
Y_c=rho_c(2);
Z_c=rho_c(3);

dX_c=drho_c(1);
dY_c=drho_c(2);
dZ_c=drho_c(3);

Q_c=((norm(rho_c))^2+2*transpose([r_v;0;0])*rho_c)/r_v^2 % EQ 13 Pontani

d2X_c= (mu/r_v^2)*(Q_c*(2+Q_c+(1+Q_c)^0.5))/((1+Q_c)^(3/2)*((1+Q_c)^0.5+1))-mu*X_c/((norm(r_cVec_v))^3)+d2theta_v*Y_c+2*dtheta_v*dY_c+dtheta_v^2*X_c; %EQ 16 Pontani
d2Y_c= -mu*Y_c/((norm(r_cVec_v))^3)-d2theta_v*X_c-2*dtheta_v*dX_c+dtheta_v^2*Y_c; %EQ 17 Pontani
d2Z_c= -mu*Z_c/((norm(r_cVec_v))^3); %EQ 18 Pontani

d2rho_c=[d2X_c;d2Y_c;d2Z_c]; %d2rho 

g_c=mc*(-mu*r_cVec_v/(norm(r_cVec_v)^3)+mu*r_tVec_v/(norm(r_tVec_v)^3)) %% need to rewrite!

%%% target

X_t=rho_t(1);
Y_t=rho_t(2);
Z_t=rho_t(3);

dX_t=drho_t(1);
dY_t=drho_t(2);
dZ_t=drho_t(3);

Q_t=((norm(rho_t))^2+2*transpose([r_v;0;0])*rho_t)/r_v^2; % EQ 13 Pontani

d2X_t= (mu/r_v^2)*(Q_t*(2+Q_t+(1+Q_t)^0.5))/((1+Q_t)^(3/2)*((1+Q_t)^0.5+1))-mu*X_t/(norm(r_tVec_v)^3)+d2theta_v*Y_t+2*dtheta_v*dY_t+dtheta_v^2*X_t; %EQ 16 Pontani
d2Y_t= -mu*Y_t/(norm(r_tVec_v)^3)-d2theta_v*X_t-2*dtheta_v*dX_t+dtheta_v^2*Y_t; %EQ 17 Pontani
d2Z_t= -mu*Z_t/(norm(r_tVec_v)^3); %EQ 18 Pontani

d2rho_t=[d2X_t;d2Y_t;d2Z_t]; %d2rho 




%% Relative Attitude Dynamics


 
w_t_tilde = [0, -w_t(3), w_t(2); %EQ 8 Liu
             w_t(3), 0, -w_t(1);
              -w_t(2), w_t(1), 0];
           
G_t=[0, -transpose(w_t); %EQ 8 Liu
    w_t, -w_t_tilde];

q_r=quatProd(quatRecip(q_t),q_c); %EQ 13 Liu solved for q_c


Ac=(q_c(1)^2-transpose(q_c(2:4))*q_c(2:4))*eye(3)+2*(q_c(2:4))*transpose(q_c(2:4))-2*q_c(1)*[0, -q_c(4), q_c(3); %EQ 10 Liu
                                                                                            q_c(4), 0, -q_c(2);
                                                                                           -q_c(3), q_c(2), 0];  
                                                                                       
At=(q_t(1)^2-transpose(q_t(2:4))*q_t(2:4))*eye(3)+2*(q_t(2:4))*transpose(q_t(2:4))-2*q_t(1)*[0, -q_t(4), q_t(3); %EQ 10 Liu
                                                                                            q_t(4), 0, -q_t(2);
                                                                                            -q_t(3), q_t(2), 0];
                                                                                        
Ar=(q_r(1)^2-transpose(q_r(2:4))*q_r(2:4))*eye(3)+2*(q_r(2:4))*transpose(q_r(2:4))-2*q_r(1)*[0, -q_r(4), q_r(3); %EQ 10 Liu
                                                                                            q_r(4), 0, -q_r(2);
                                                                                            -q_r(3), q_r(2), 0];                                                                                        
                                                                             
Act=Ac*At; %EQ 10 Liu



w_r=w_c-Ar*w_t; %EQ 11 Liu solved for w_c
% w_rout=[w_rout w_r];
w_c_tilde = [0, -w_c(3), w_c(2); %EQ 12 Liu
             w_c(3), 0, -w_c(1);
              -w_c(2), w_c(1), 0];

% 
w_r_tilde = [0, -w_r(3), w_r(2); %EQ 12 Liu
             w_r(3), 0, -w_r(1);
              -w_r(2), w_r(1), 0];

% 
% dw_t=It\(-w_t_tilde*It*w_t); %EQ 6 Liu
% tc=Ic\(-w_c_tilde*Ic*w_c)-Ar*dw_t+w_r_tilde*w_c;
%                                                 
% tc_tilde = [0, -tc(3), tc(2); %EQ 12 Liu
%              tc(3), 0, -tc(1);
%               -tc(2), tc(1), 0];
% 
% 
G_r=[0, -transpose(w_r);
    w_r, -w_r_tilde]

G_c=[0, -transpose(w_c);
    w_c, -w_c_tilde]

dq_r=0.5*G_r*q_r

% B1=[0, -transpose(tc); %EQ 16 Liu
%     tc, -tc_tilde];
% B2=[0, -transpose(w_r); %EQ 16 Liu
%     w_r, -w_r_tilde];
% B3=[-transpose(q_r(2:4)); %EQ 16 Liu
%     q_r(1)*eye(3)+[0, -q_r(4), q_r(3);
%                    q_r(4), 0, -q_r(2);
%                    -q_r(3), q_r(2), 0]];
%                
% 
% f=[A1*rho+A2*drho+A3;
%    0.5*B1*q_r+0.5*B2*dq_r] %% Eq 16 Liu with Liu coordinates

% f=[d2rho;   
%    0.5*B1*q_r+0.5*B2*dq_r]; %% Eq 16 Liu with pontani coordinates


%% LI Control

% Relative position
Ch=0.5
rho_obs=rho_d
d_0=1.9
delta_0=2
alpha=0.5
etaPos=5e-4


if norm(rho-rho_obs)>=delta_0
    k=0
else
    k=0
end

gradUatt = Ch*(rho-rho_d)/sqrt((norm(rho-rho_d))^2+1)

gradUrep = 4*k*(rho-rho_obs)*((norm(rho-rho_obs))^2-delta_0^2)/((norm(rho-rho_obs))^2-d_0^2)^3

sRel=drho+alpha*(gradUatt+gradUrep)

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

grad2Uatt=Ch*(((norm(rho-rho_d))^2+1)*eye(3)-(rho-rho_d)*transpose(rho-rho_d))/((norm(rho-rho_d))^2+1)^(3/2)





gamma1=0.02 ;
gamma2=0.005;
gamma3=0.05;
pdiag=[0.5, 0, 0; 0, 0.5, 0; 0, 0, 0.5];
mdiag=[0.5, 0, 0; 0, 0.5, 0; 0, 0, 0.5];
% spos=drho; % s position control, excluding APF term  Eq 20 Li
% satt=derror(2:4); % s att control, excluding APF term Eq 36 

Fc=-mc*alpha*grad2Uatt*drho-(gamma1+mc*gamma3)*satPos-pdiag*satPos;  %Eq 27 Li excluding APF terms


% Relative attitude 
Chrel=0.5
etaAtt=5e-4
beta=0.6
gradUquat=Chrel*(q_r(2:4)-[0;0;0])/sqrt((norm(q_r(2:4)-[0;0;0]))^2+1)

sAtt=dq_r(2:4)+beta*gradUquat

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
% satAtt=[tanh(sAtt(1)); tanh(sAtt(2)); tanh(sAtt(3))]

grad2Uquat=Chrel*(((norm(q_r(2:4)-[0;0;0]))^2+1)*eye(3)-(q_r(2:4)-[0;0;0])*transpose(q_r(2:4)-[0;0;0]))/((norm(q_r(2:4)-[0;0;0]))^2+1)^(3/2)

T=[0, -q_r(4), q_r(3); q_r(4), 0, -q_r(2);-q_r(3), q_r(2), 0] + q_r(1)*eye(3)
dT=[0, -dq_r(4), dq_r(3); dq_r(4), 0, -dq_r(2);-dq_r(3), dq_r(2), 0] + dq_r(1)*eye(3)
P=inv(T)
dP=inv(dT)
Istar=transpose(P)*Ic*P

IcPdq_r=Ic*P*dq_r(2:4)
Pdq_r=P*dq_r(2:4)
Arw_t=Ar*w_t
twoPdq_r=2*P*dq_r(2:4)
dw_t=It\(-w_t_tilde*It*w_t)

Cstar=-Istar/dP*P-2*transpose(P)*[0, -IcPdq_r(3), IcPdq_r(2); IcPdq_r(3), 0, -IcPdq_r(1); -IcPdq_r(2), IcPdq_r(1), 0]*P


Nstar=transpose(P)*([0, -Pdq_r(3), Pdq_r(2); Pdq_r(3), 0, -Pdq_r(1); -Pdq_r(2), Pdq_r(1), 0]*Ic*Ar*w_t)...
    +transpose(P)*([0, -Arw_t(3), Arw_t(2); Arw_t(3), 0, -Arw_t(1);-Arw_t(2), Arw_t(1), 0]*Ic*P*dq_r(2:4))+...
    0.5*transpose(P)*([0, -Arw_t(3), Arw_t(2); Arw_t(3), 0, -Arw_t(1); -Arw_t(2), Arw_t(1), 0]*Ic*Ar*w_t)-...
    0.5*transpose(P)*Ic*([0, -twoPdq_r(3), twoPdq_r(2); twoPdq_r(3), 0, -twoPdq_r(1);-twoPdq_r(2), twoPdq_r(1), 0]*Ar*w_t-Ar*dw_t)

Tc=Nstar-Cstar*beta*gradUquat-Istar*beta*grad2Uquat*dq_r(2:4)-gamma2*satAtt-mdiag*satAtt
Tcprime = 2*transpose(T)*Tc
% TOUT=[TOUT, Tcprime];
% TEST=0.5*transpose(P)*Tcprime


                            

% Tc=-gamma2*sign(satt)-mdiag*sign(satt); %Eq 37 excluding APF terms
%  
% T=[0, -q_r(4), q_r(3);
%    q_r(4), 0, -q_r(2);
%     -q_r(3), q_r(2), 0] + q_r(1)*eye(3);
% P=inv(T);
% 
% Jstar=transpose(P)*Ic*P;
% 
% B=[(1/mc)*eye(3), zeros(3); % B matrix.. attitude term is wrong
%     zeros(4,3), B3/Ic] ;
% u=[Fc;Tc];
% 


 logicMat2=[logicMat2, thrusterLogic2(Fc,Tcprime,q_c,theta_v, i_v, omega_v)]


%% Derivatives
9000000000
t
der(1,1)=sqrt(mu/p_v^3)*(1+e_v*cos(theta_v))^2 % derivative True anamoly 
der(2:4,1)=It\(-w_t_tilde*It*w_t) % derivative target ang velocity
der(5:8,1)=0.5*G_t*q_t % derivative target attitude
% if t<150
%     der(9:11,1)=Ic\(-w_c_tilde*Ic*w_c)+Ic\Tcprime % derivative relative angular velocity 
% else
%     der(9:11,1)=Ic\(-w_c_tilde*Ic*w_c) % derivative relative angular velocity 
% end
der(9:11,1)=Ic\(-w_c_tilde*Ic*w_c)+Ic\Tcprime % derivative relative angular velocity
der(12:15,1)=0.5*G_c*q_c
der(16:18,1)=drho_c
% if t<95
%     der(19:21,1)= g/mc%+Fc/mc %% simulation slow after convergence.. attitude converge takes longer though..
% else
%     der(19:21,1)= g/mc
% end
der(19:21,1)=d_c+Fc/mc+g_c/mc
der(22:24,1)=drho_t
der(25:27,1)=d2rho_t+d_t
% der(12:18,1)=derror % derivative error
% der(19:25,1)=f%+B*u % derivative velocity error
% der(16:19,1)=dq_r
% der(20:23,1)=0.5*B1*q_r+0.5*B2*dq_r
