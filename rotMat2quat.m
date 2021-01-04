function quat = rotMat2quat(a,b,c,d,e,f,g,h,i)
R_IB = [a,b,c;
        d,e,f;
        g,h,i]
  
  q0 = 0.5*sqrt(R_IB(1,1)+R_IB(2,2)+R_IB(3,3)+1)
  q1 = (R_IB(2,3)-R_IB(3,2))/(4*q0)
  q2 = (R_IB(3,1)-R_IB(1,3))/(4*q0)
  q3 = (R_IB(1,2)-R_IB(2,1))/(4*q0)
  
  quat=[q0;q1;q2;q3]