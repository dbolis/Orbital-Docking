function quat = euler2quat(theta,phi,psi)
R_IB = [cos(phi)*cos(psi)-sin(phi)*cos(theta)*sin(psi), sin(phi)*cos(theta)*cos(psi)+cos(phi)*sin(psi), sin(phi)*sin(theta);
      -cos(phi)*cos(theta)*sin(psi)-sin(phi)*cos(psi), cos(phi)*cos(theta)*cos(psi)-sin(phi)*sin(psi), cos(phi)*sin(theta);
      sin(theta)*sin(psi), -sin(theta)*cos(psi), cos(theta)]
  
  q0 = 0.5*sqrt(R_IB(1,1)+R_IB(2,2)+R_IB(3,3)+1)
  q1 = (R_IB(2,3)-R_IB(3,2))/(4*q0)
  q2 = (R_IB(3,1)-R_IB(1,3))/(4*q0)
  q3 = (R_IB(1,2)-R_IB(2,1))/(4*q0)
  
  quat=[q0,q1,q2,q3]