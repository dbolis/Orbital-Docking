function out=thrusterLogic(force,torque,q_c,theta, i, omega)

force_c=force*((q_c(1)^2-transpose(q_c(2:4))*q_c(2:4))*eye(3)+2*(q_c(2:4))*transpose(q_c(2:4))-2*q_c(1)*[0, -q_c(4), q_c(3)
                                                                                            q_c(4), 0, -q_c(2)
                                                                                            -q_c(3), q_c(2), 0])*ECI2LVLH313(theta, i, omega) 

                                                                                        
rx=1
ry=1
rz=1
                                                                                     
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
            
                                                                                        
                                                                                        
% if force_c(1)>0
%     forceOut(1,:)=["1x", "2x", force_c(1)/2]
% else
%     forceOut(1,:)=["3x", "4x", force_c(1)/2]
% end
% 
% if force_c(2)>0
%     forceOut(2,:)=["1y", "2y", force_c(2)/2]
% else
%     forceOut(2,:)=["3y", "4y", force_c(2)/2]
% end
%          
% if force_c(3)>0
%     forceOut(3,:)=["1z", "2z", force_c(3)/2]
% else
%     forceOut(3,:)=["3z", "4z", force_c(3)/2]
% end

if force_c(1)>0
    thrustVals(3)= force_c(1)/2
    thrustVals(4)= force_c(1)/2
else
    thrustVals(1)= abs(force_c(1)/2)
    thrustVals(2)= abs(force_c(1)/2)
end

if force_c(2)>0
    thrustVals(7)= force_c(2)/2
    thrustVals(8)= force_c(2)/2
else
    thrustVals(5)= abs(force_c(2)/2)
    thrustVals(6)= abs(force_c(2)/2)
end
         
if force_c(3)>0
    thrustVals(11)= force_c(3)/2
    thrustVals(12)= force_c(3)/2
else
    thrustVals(9)= abs(force_c(3)/2)
    thrustVals(10)= abs(force_c(3)/2)
end

    
if torque(1)>0
    thrustVals(5)= thrustVals(5)+torque(1)/(2*rx) %y1
    thrustVals(7)= thrustVals(7)+torque(1)/(2*rx) %y3
else
    thrustVals(6)= thrustVals(6)+abs(torque(1)/(2*rx)) %y2
    thrustVals(8)= thrustVals(8)+abs(torque(1)/(2*rx)) %y4
end

if torque(2)>0
    thrustVals(9)= thrustVals(9)+torque(2)/(2*ry) %z1
    thrustVals(11)= thrustVals(11)+torque(2)/(2*ry) %z3
else
    thrustVals(10)= thrustVals(10)+abs(torque(2)/(2*ry)) %z2
    thrustVals(12)= thrustVals(12)+abs(torque(2)/(2*ry)) %z4
end
         
if torque(3)>0
    thrustVals(1)= thrustVals(1)+torque(3)/(2*rz)
    thrustVals(3)= thrustVals(3)+torque(3)/(2*rz)
else
    thrustVals(2)= thrustVals(2)+abs(torque(3)/(2*rz))
    thrustVals(4)= thrustVals(4)+abs(torque(3)/(2*rz))
end    
    

out=[thrustName, thrustVals]