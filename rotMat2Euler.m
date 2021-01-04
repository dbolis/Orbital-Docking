function Euler = rotMat2Euler(a,b,c,d,e,f,g,h,i)
R_IB = [a,b,c;
        d,e,f;
        g,h,i]
    
theta = acos(R_IB(3,3))

sphi = R_IB(1,3)/sin(theta)
cphi = R_IB(2,3)/sin(theta)
phi = 2*atan(sphi/(cphi+1))

spsi = R_IB(3,1)/sin(theta)
cpsi = -R_IB(3,2)/sin(theta)
psi=2*atan(spsi/(cpsi+1))


Euler = [theta,phi,psi]