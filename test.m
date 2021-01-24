% theta_v=3
% mu=3.986e14
% p_v=7178169
% e_v=0
% i_v=0.2
% omega_v=0.1
% rho_c=[0;0;0]
% drho_c=[0.58;-0.7;1]
% rho_t=[-20;12;-7]
% drho_t=[0.5;-0.7;1]
% Re=6.3781e6
% J2=1.083e-3
% omega_E=7.2921159e-5 % ang velocity of earth around its axis
% c_D=2.2 % drag coefficient
% AerS=pi*4^2 % need to add for both chaser and target
% mc=20 % chaser mass [kg]
% 
% dtheta_v=sqrt(mu/p_v^3)*(1+e_v*cos(theta_v))^2
% d2theta_v=-2*(mu/p_v^3)*(e_v*sin(theta_v))*(1+e_v*cos(theta_v))^3
% 
% %%% Chaser
% 
% r_v=p_v/(1+e_v*cos(theta_v)) % virtual radius from CoM
% r_cVec_v=[r_v;0;0]+rho_c % chaser radius vector in virtual frame
% v_cVec_v=[sqrt(mu/p_v)*e_v*sin(theta_v)+drho_c(1)-dtheta_v*rho_c(2); sqrt(mu/p_v)*(1+e_v*cos(theta_v))+drho_c(2)+dtheta_v*rho_c(1); drho_c(3)] % chaser velocity vector in virtual frame
% 
% r_cVec_ECI=transpose(r_cVec_v)*ECI2LVLH313(theta_v,i_v,omega_v) % chaser radius vector in ECI
% v_cVec_ECI=transpose(v_cVec_v)*ECI2LVLH313(theta_v,i_v,omega_v) % chaser velocity vector in ECI
% 
% 
% h_cVec=cross(r_cVec_ECI,v_cVec_ECI) % specific angular momentum vector
% h_cVers=h_cVec/norm(h_cVec) % specific angular momentum versor
% 
% r_cVers=r_cVec_ECI/norm(r_cVec_ECI) % chaser radius versor in ECI
% theta_cVers=cross(h_cVers,r_cVers) % theta versor
% 
% phi_c=asin(r_cVers(3)) % latitude
% clambda_c=r_cVers(1)/cos(phi_c) 
% slambda_c=r_cVers(2)/cos(phi_c) 
% 
% lambda_c=2*atan(slambda_c/(1+clambda_c)) % absolute longitutude
% 
% szeta_c=theta_cVers(3)/cos(phi_c) 
% czeta_c=h_cVers(3)/cos(phi_c)
% 
% zeta_c=2*atan(szeta_c/(1+czeta_c)) % heading angle
% 
% fJ2_c_c=3*Re^2*J2*mu/(norm(r_cVec_ECI))^4*[(3*sin(phi_c)^2-1)/2, -sin(phi_c)*cos(phi_c)*sin(zeta_c), -sin(phi_c)*cos(phi_c)*cos(zeta_c)] % J2 disturbance on chaser in chaser frame
% 
% fJ2_c_v=fJ2_c_c*ECI2LVLH123(zeta_c, phi_c, lambda_c)*transpose(ECI2LVLH313(theta_v,i_v,omega_v)) % J2 disturbance on chaser in virtual frame
% 
% 
% %%% Drag
% 
% a_c=mu/((2*mu/norm(r_cVec_ECI))-norm(v_cVec_ECI)^2)
% e_c=sqrt(1-(norm(h_cVec)^2/(mu*a_c)))
% v_radial_c=dot(v_cVec_ECI,r_cVec_ECI)/norm(r_cVec_ECI)
% p_c=a_c*(1-e_c^2)
% 
% cThetaStar_c=1/e_c*(p_c/norm(r_cVec_ECI)-1)
% sThetaStar_c=v_radial_c/e_c*sqrt(p_c/mu)
% 
% theta_c=2*atan(sThetaStar_c/(1+cThetaStar_c))
% 
% 
% v_r_c=sqrt(mu/p_c)*e_c*sin(theta_c)
% v_e_c=sqrt(mu/p_c)*(1+e_c*cos(theta_c))*cos(zeta_c)
% v_n_c=sqrt(mu/p_c)*(1+e_c*cos(theta_c))*sin(zeta_c)
% 
% v_r_cVec=[v_r_c, v_e_c-omega_E*norm(r_cVec_ECI)*cos(phi_c) v_n_c]*transpose([1 0, 0; 0 cos(zeta_c) -sin(zeta_c); 0 sin(zeta_c) cos(zeta_c)])
% 
% atmDensity_c=1.225*exp(-(norm(r_cVec_ECI)-Re)/10332.6) %% apply better fit.. negligible at 800 km altitude
% 
% fDrag_c_c=-0.5*c_D*AerS/mc*atmDensity_c*norm(v_r_cVec)*v_r_cVec % in LVLH chaser frame
% fDrag_c_v=fDrag_c_c*ECI2LVLH123(zeta_c, phi_c, lambda_c)*transpose(ECI2LVLH313(theta_v,i_v,omega_v)) % in LVLH virtual frame
% 
% 
% %%% target
% 
% r_tVec_v=[r_v;0;0]+rho_t  % target radius vector in virtual frame
% v_tVec_v=[sqrt(mu/p_v)*e_v*sin(theta_v)+drho_t(1)-dtheta_v*rho_t(2); sqrt(mu/p_v)*(1+e_v*cos(theta_v))+drho_t(2)+dtheta_v*rho_t(1); drho_t(3)] % target velocity vector in virtual frame
% 
% r_tVec_ECI=transpose(r_tVec_v)*ECI2LVLH313(theta_v,i_v,omega_v) % target radius vector in ECI
% v_tVec_ECI=transpose(v_tVec_v)*ECI2LVLH313(theta_v,i_v,omega_v) % target velocity vector in ECI
% 
% 
% h_tVec=cross(r_tVec_ECI,v_tVec_ECI) % specific angular momentum vector
% h_tVers=h_tVec/norm(h_tVec) % specific angular momentum versor
% 
% r_tVers=r_tVec_ECI/norm(r_tVec_ECI) % target radius versor in ECI
% theta_tVers=cross(h_tVers,r_tVers) % theta versor 
% 
% phi_t=asin(r_tVers(3)) % latitude
% clambda_t=r_tVers(1)/cos(phi_t)
% slambda_t=r_tVers(2)/cos(phi_t)
% 
% lambda_t=2*atan(slambda_t/(1+clambda_t)) % absolute longitutude
% 
% szeta_t=theta_tVers(3)/cos(phi_t)
% czeta_t=h_tVers(3)/cos(phi_t)
% 
% zeta_t=2*atan(szeta_t/(1+czeta_t)) % heading angle
% 
% %%J2
% fJ2_t_t=3*Re^2*J2*mu/(norm(r_tVec_ECI))^4*[(3*sin(phi_t)^2-1)/2, -sin(phi_t)*cos(phi_t)*sin(zeta_t), -sin(phi_t)*cos(phi_t)*cos(zeta_t)] % J2 disturbance on target in target frame
% 
% fJ2_t_v=fJ2_t_t*ECI2LVLH123(zeta_t, phi_t, lambda_t)*transpose(ECI2LVLH313(theta_v,i_v,omega_v)) % J2 disturbance on target in virtual frame
% 
% %%% Drag
% a_t=mu/((2*mu/norm(r_tVec_ECI))-norm(v_tVec_ECI)^2)
% e_t=sqrt(1-(norm(h_tVec)^2/(mu*a_t)))
% v_radial_t=dot(v_tVec_ECI,r_tVec_ECI)/norm(r_tVec_ECI)
% p_t=a_t*(1-e_t^2)
% 
% cThetaStar_t=1/e_t*(p_t/norm(r_tVec_ECI)-1)
% sThetaStar_t=v_radial_t/e_t*sqrt(p_t/mu)
% 
% theta_t=2*atan(sThetaStar_t/(1+cThetaStar_t))
% 
% v_r_t=sqrt(mu/p_t)*e_t*sin(theta_t)
% v_e_t=sqrt(mu/p_t)*(1+e_t*cos(theta_t))*cos(zeta_t)
% v_n_t=sqrt(mu/p_t)*(1+e_t*cos(theta_t))*sin(zeta_t)
% 
% v_r_tVec=[v_r_t, v_e_t-omega_E*norm(r_tVec_ECI)*cos(phi_t) v_n_t]*transpose([1 0, 0; 0 cos(zeta_t) -sin(zeta_t); 0 sin(zeta_t) cos(zeta_t)])
% 
% atmDensity_t=1.225*exp(-(norm(r_tVec_ECI)-Re)/10332.6)
% 
% fDrag_t_t=-0.5*c_D*AerS/mc*atmDensity_t*norm(v_r_tVec)*v_r_tVec
% fDrag_t_v=fDrag_t_t*ECI2LVLH123(zeta_t, phi_t, lambda_t)*transpose(ECI2LVLH313(theta_v,i_v,omega_v))

clear all;

syms Rx Ry Rz Fx Fy Fz Tx Ty Tz;

A=[-1,-1,1,1,0,0,0,0,0,0,0,0;...
    0,0,0,0,-1,-1,1,1,0,0,0,0;...
    0,0,0,0,0,0,0,0,-1,-1,1,1;...
    0,0,0,0,Ry,-Ry,Ry,-Ry,0,0,0,0;...
    0,0,0,0,0,0,0,0,Rz,-Rz,Rz,-Rz;...
    Rx,-Rx,Rx,-Rx,0,0,0,0,0,0,0,0];

b=[Fx; Fy; Fz; Tx; Ty; Tz]

transpose(A)*(A*transpose(A))^(-1)

