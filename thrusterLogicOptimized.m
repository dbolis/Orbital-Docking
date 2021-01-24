function out=thrusterLogicOptimized(force,torque,q_c,theta, i, omega)

force_c=transpose(force)*ECI2LVLH313(theta, i, omega)*transpose((q_c(1)^2-transpose(q_c(2:4))*q_c(2:4))*eye(3)+2*(q_c(2:4))*transpose(q_c(2:4))-2*q_c(1)*[0, -q_c(4), q_c(3)
                                                                                            q_c(4), 0, -q_c(2)
                                                                                            -q_c(3), q_c(2), 0])

                                                                                        
Rx=1
Ry=1
Rz=1
                                                                                     
thrustName=["1x";
        "2x";
        "3x";
        "4x";
        "1y";
        "2y";
        "3y";
        "4y";
        "1z";
        "2z";
        "3z";
        "4z"]
    
    
A=[-1,-1,1,1,0,0,0,0,0,0,0,0;
    0,0,0,0,-1,-1,1,1,0,0,0,0;
    0,0,0,0,0,0,0,0,-1,-1,1,1;
    0,0,0,0,Ry,-Ry,Ry,-Ry,0,0,0,0;
    0,0,0,0,0,0,0,0,Rz,-Rz,Rz,-Rz;
    Rx,-Rx,Rx,-Rx,0,0,0,0,0,0,0,0];    
    
b=[transpose(force_c);torque]
    
thrustVals=transpose(A)*(A*transpose(A))^(-1)*b
            
                                                                                        
out=[thrustName, thrustVals]
% out=thrustVals