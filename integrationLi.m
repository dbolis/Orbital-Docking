function der = integrationLi(t,x,e_v,p_v,e_c,n_c,e_t,n_t,p_c,p_t,Ic,It,mc,mu,beta,p,q,eta,q_d,rho_d,drho_d)
theta_t=x(1); % True anamoly 
w_t=x(2:4); % target ang velocity
q_t=x(5:8); % target attitude
w_c=x(9:11);% relative angular velocity
q_c=x(12:15)
rho=x(16:18)
drho=x(19:21)
theta_v=x(22)
rho_c=x(23:25)
drho_c=x(26:28)


% error=x(12:18); % error
% derror=x(19:25); % velocity error
% q_r=x(16:19)
% dq_r=x(20:23)
% q_r=error(4:7)+q_d; % relative attitude 
% dq_r=derror(4:7); % relative ang rate
% rho=error(1:3)+rho_d; % relative position
% drho=derror(1:3)+drho_d; % relative velocity


global count w_rout TOUT
count=count+1



%% Dynamics with Perterbances
dtheta_v=sqrt(mu/p_v^3)*(1+e_v*cos(theta_v))^2
d2theta_v=-2*(mu/p_v^3)*(e_v*sin(theta_v))*(1+e_v*cos(theta_v))^3

%%% Chaser

r_v=p_v/(1+e_v*cos(theta_v)) %% virtual radius from CoM
r_cVec=[r_v;0;0]+rho_c %% in virtual coordinates
v_cVec=[sqrt(mu/p_v)*e_v*sin(theta_v)+drho_c(1)-dtheta_v*rho_c(2); sqrt(mu/p_v)*(1+e_v*cos(theta_v))+drho_c(2)+dtheta_v*rho_c(1); drho_c(3)] 




%% Relative Orbit Dynamics

% dtheta_t=(n_t*(1+e_t*cos(theta_t))^2)/((1-e_t^2)^(3/2)); % target orbital velocity
% d2theta_t=(-2*n_t^2*e_t*sin(theta_t)*(1+e_t*cos(theta_t))^3)/((1-e_t^2)^3); % target orbital acceleration

dtheta_t=sqrt(mu/p_t^3)*(1+e_t*cos(theta_t))^2
d2theta_t=-2*(mu/p_t^3)*(e_t*sin(theta_t))*(1+e_t*cos(theta_t))^3


r_t=p_t/(1+e_t*cos(theta_t)); % radius from Earth CoM to target
r_cVec= [r_t;0;0]+rho; % radius from Earth CoM to chaser VECTOR Pontani coordinates
% r_cVec= [0;0;r_t]+rho % radius from Earth CoM to chaser VECTOR Liu coordinates
r_c= norm(r_cVec); % % radius from Earth CoM to chaser scalar



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

X=rho(1);
Y=rho(2);
Z=rho(3);

dX=drho(1);
dY=drho(2);
dZ=drho(3);

Q=((norm(rho))^2+2*transpose([r_t;0;0])*rho)/r_t^2; % EQ 13 Pontani

d2X= (mu/r_t^2)*(Q*(2+Q+(1+Q)^0.5))/((1+Q)^(3/2)*((1+Q)^0.5+1))-mu*X/(r_c^3)+d2theta_t*Y+2*dtheta_t*dY+dtheta_t^2*X; %EQ 16 Pontani
d2Y= -mu*Y/r_c^3-d2theta_t*X-2*dtheta_t*dX+dtheta_t^2*Y; %EQ 17 Pontani
d2Z= -mu*Z/r_c^3; %EQ 18 Pontani

d2rho=[d2X;d2Y;d2Z]; %d2rho 

g=mc*(-mu*r_cVec/r_c^3+mu*[r_t;0;0]/r_t^3)

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
mdiag=[0.3, 0, 0; 0, 0.3, 0; 0, 0, 0.3];
% spos=drho; % s position control, excluding APF term  Eq 20 Li
% satt=derror(2:4); % s att control, excluding APF term Eq 36 

Fc=-mc*alpha*grad2Uatt*drho-(gamma1+mc*gamma3)*satPos-pdiag*satPos;  %Eq 27 Li excluding APF terms


% Relative attitude 
Chrel=0.3
etaAtt=5e-4
beta=0.5
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



%% Derivatives

t
der(1,1)=sqrt(mu/p_t^3)*(1+e_t*cos(theta_t))^2 % derivative True anamoly 
der(2:4,1)=It\(-w_t_tilde*It*w_t) % derivative target ang velocity
der(5:8,1)=0.5*G_t*q_t % derivative target attitude
% if t<150
%     der(9:11,1)=Ic\(-w_c_tilde*Ic*w_c)+Ic\Tcprime % derivative relative angular velocity 
% else
%     der(9:11,1)=Ic\(-w_c_tilde*Ic*w_c) % derivative relative angular velocity 
% end
der(9:11,1)=Ic\(-w_c_tilde*Ic*w_c)%+Ic\Tcprime % derivative relative angular velocity
der(12:15,1)=0.5*G_c*q_c
der(16:18,1)=drho
if t<95
    der(19:21,1)= g/mc%+Fc/mc %% simulation slow after convergence.. attitude converge takes longer though..
else
    der(19:21,1)= g/mc
end
der(22,1)=dtheta_v
der(23:25,1)=0 %placeholder
der(26:28,1)=0 %placeholder
% der(12:18,1)=derror % derivative error
% der(19:25,1)=f%+B*u % derivative velocity error
% der(16:19,1)=dq_r
% der(20:23,1)=0.5*B1*q_r+0.5*B2*dq_r
