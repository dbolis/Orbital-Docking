theta_v=3
mu=3.986e14
p_v=7178169
e_v=0
i_v=0.2
omega_v=0.1
rho_c=[-20;12;-7]
drho_c=[0.5;-0.7;1]
Re=6.3781e6
J2=1.083e-3

dtheta_v=sqrt(mu/p_v^3)*(1+e_v*cos(theta_v))^2
d2theta_v=-2*(mu/p_v^3)*(e_v*sin(theta_v))*(1+e_v*cos(theta_v))^3

%%% Chaser

r_v=p_v/(1+e_v*cos(theta_v)) %% virtual radius from CoM
r_cVec_v=[r_v;0;0]+rho_c %% in virtual coordinates
v_cVec_v=[sqrt(mu/p_v)*e_v*sin(theta_v)+drho_c(1)-dtheta_v*rho_c(2); sqrt(mu/p_v)*(1+e_v*cos(theta_v))+drho_c(2)+dtheta_v*rho_c(1); drho_c(3)] %% in virtual coordinates

r_cVec_ECI=transpose(r_cVec_v)*ECI2LVLH313(theta_v,i_v,omega_v)
v_cVec_ECI=transpose(v_cVec_v)*ECI2LVLH313(theta_v,i_v,omega_v)


h_cVec=cross(r_cVec_ECI,v_cVec_ECI)
h_cVers=h_cVec/norm(h_cVec)

r_cVers=r_cVec_ECI/norm(r_cVec_ECI)
theta_cVers=cross(h_cVers,r_cVers)

phi_c=asin(r_cVers(3))
clambda_c=r_cVers(1)/cos(phi_c)
slambda_c=r_cVers(2)/cos(phi_c)

lambda_c=2*atan(slambda_c/(1+clambda_c))

szeta_c=theta_cVers(3)/cos(phi_c)
czeta_c=h_cVers(3)/cos(phi_c)

zeta_c=2*atan(szeta_c/(1+czeta_c))

fJ2_c_c=3*Re^2*J2*mu/(norm(r_cVec_ECI))^4*[(3*sin(phi_c)^2-1)/2, -sin(phi_c)*cos(phi_c)*sin(zeta_c), -sin(phi_c)*cos(phi_c)*cos(zeta_c)]

fJ2_c_v=fJ2_c_c*ECI2LVLH123(zeta_c, phi_c, lambda_c)*transpose(ECI2LVLH313(theta_v,i_v,omega_v))


