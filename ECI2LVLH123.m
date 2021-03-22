function out=ECI2LVLH123(zeta, phi, lambda)
% 321 rotation matrix
out=[cos(phi)*cos(lambda), cos(phi)*sin(lambda), sin(phi);
    -sin(zeta)*sin(phi)*cos(lambda)-cos(zeta)*sin(lambda), -sin(zeta)*sin(phi)*sin(lambda)+cos(zeta)*cos(lambda), sin(zeta)*cos(phi);
    -cos(zeta)*sin(phi)*cos(lambda)+sin(zeta)*sin(lambda), -cos(zeta)*sin(phi)*sin(lambda)-sin(zeta)*cos(lambda), cos(zeta)*cos(phi)]
