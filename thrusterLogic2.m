function out=thrusterLogic2(force,torque,q_c,theta, i, omega)

force_c=transpose(force)*ECI2LVLH313(theta, i, omega)*transpose((q_c(1)^2-transpose(q_c(2:4))*q_c(2:4))*eye(3)+2*(q_c(2:4))*transpose(q_c(2:4))-2*q_c(1)*[0, -q_c(4), q_c(3)
                                                                                            q_c(4), 0, -q_c(2)
                                                                                            -q_c(3), q_c(2), 0])

                                                                                        
Rx=5
Ry=5
Rz=5
                                                                                     
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
    
thrustVals=[0;
            0;
            0;
            0;
            0;
            0;
            0;
            0;
            0;
            0;
            0;
            0]

  
        
        
            
                                                                                        
                                                                                        
A=[-1,-1,0,0,0,0;
    0,0,-1,-1,0,0;
    0,0,0,0,-1,-1;
    0,0,Ry,-Ry,0,0;
    0,0,0,0,Rz,-Rz;
    Rx,-Rx,0,0,0,0];    
    
b=[transpose(force_c);torque]
    
y=A\b

if y(1)>0
    thrustVals(1)=y(1)
else
    thrustVals(4)=abs(y(1))
end

if y(2)>0
    thrustVals(2)=y(2)
else
    thrustVals(3)=abs(y(2))
end

if y(3)>0
    thrustVals(5)=y(3)
else
    thrustVals(8)=abs(y(3))
end

if y(4)>0
    thrustVals(6)=y(4)
else
    thrustVals(7)=abs(y(4))
end

if y(5)>0
    thrustVals(9)=y(5)
else
    thrustVals(12)=abs(y(5))
end

if y(6)>0
    thrustVals(10)=y(6)
else
    thrustVals(11)=abs(y(6))
end

% out=[thrustName, thrustVals]
out=thrustVals