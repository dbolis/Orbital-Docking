function out=thruster2dynamics(thrusterVec,q_c,theta, i, omega)

Rx=5
Ry=5
Rz=5



t1x=thrusterVec(1)
t2x=thrusterVec(2)
t3x=thrusterVec(3)
t4x=thrusterVec(4)
t1y=thrusterVec(5)
t2y=thrusterVec(6)
t3y=thrusterVec(7)
t4y=thrusterVec(8)
t1z=thrusterVec(9)
t2z=thrusterVec(10)
t3z=thrusterVec(11)
t4z=thrusterVec(12)

Fx=-t1x-t2x+t3x+t4x
Fy=-t1y-t2y+t3y+t4y
Fz=-t1z-t2z+t3z+t4z
Tx=Ry*t1y-Ry*t2y+Ry*t3y-Ry*t4y
Ty=Rz*t1z-Rz*t2z+Rz*t3z-Rz*t4z
Tz=Rx*t1x-Rx*t2x+Rx*t3x-Rx*t4x

force_c=[Fx, Fy, Fz]
torque_c=[Tx, Ty, Tz]


force_v=force_c*((q_c(1)^2-transpose(q_c(2:4))*q_c(2:4))*eye(3)+2*(q_c(2:4))*transpose(q_c(2:4))-2*q_c(1)*[0, -q_c(4), q_c(3)
                                                                                            q_c(4), 0, -q_c(2)
                                                                                            -q_c(3), q_c(2), 0])*transpose(ECI2LVLH313(theta, i, omega))
                               
out=[transpose(force_v); transpose(torque_c)]                                           