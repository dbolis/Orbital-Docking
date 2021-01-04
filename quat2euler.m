function Euler = quat2euler(q0,q1,q2,q3)
R_IB=[q0^2+q1^2-q2^2-q3^2, 2*(q1*q2+q0*q3), 2*(q1*q3-q0*q2);
     2*(q1*q2-q0*q3), q0^2+q2^2-q1^2-q3^2, 2*(q2*q3+q0*q1);
     2*(q1*q3+q0*q2), 2*(q2*q3-q0*q1), q0^2+q3^2-q1^2-q2^2]
      
theta = acos(R_IB(3,3))

sphi = R_IB(1,3)/sin(theta)
cphi = R_IB(2,3)/sin(theta)
phi = 2*atan(sphi/(cphi+1))

spsi = R_IB(3,1)/sin(theta)
cpsi = -R_IB(3,2)/sin(theta)
psi=2*atan(spsi/(cpsi+1))


Euler = [theta,phi,psi]