function out=genRandInit(w_t,rho_t,drho_t,w_c,rho_c,drho_c)

num=2 % number of values generated

% standard deviations

w_tSD=deg2rad(5)
rho_tSD=10
drho_tSD=1
w_cSD=deg2rad(20)
rho_cSD=10
drho_cSD=1

% generate norm distributed values with inputted means and starndard devs

w_t_init=ones(3,num).*w_t
w_t_out=normrnd(w_t_init, w_tSD)

rho_t_init=ones(3,num).*rho_t
rho_t_out=normrnd(rho_t_init, rho_tSD)

drho_t_init=ones(3,num).*drho_t
drho_t_out=normrnd(drho_t_init, drho_tSD)

w_c_init=ones(3,num).*w_c
w_c_out=normrnd(w_c_init, w_cSD)

rho_c_init=ones(3,num).*rho_c
rho_c_out=normrnd(rho_c_init, rho_cSD)

drho_c_init=ones(3,num).*drho_c
drho_c_out=normrnd(drho_c_init, drho_cSD)

% generate random quaternion

PHI=pi*rand(1,num)
zeta=2*pi*rand(1,num)
r=2*rand(1,num)-ones(1,num)
phi=asin(r)

a1=cos(zeta).*cos(phi)
a2=cos(phi).*sin(zeta)
a3=sin(phi)

q_c_out=[cos(PHI/2); a1.*sin(PHI/2); a2.*sin(PHI/2); a3.*sin(PHI/2)]


out=[w_t_out;rho_t_out;drho_t_out;w_c_out;rho_c_out;drho_c_out;q_c_out]