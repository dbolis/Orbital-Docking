function out=ECI2LVLH313(theta, i, omega)
% 313 rotation matrix
out=[cos(theta)*cos(omega)-sin(theta)*cos(i)*sin(omega), cos(theta)*sin(omega)+sin(theta)*cos(i)*cos(omega), sin(theta)*sin(i);
    -sin(theta)*cos(omega)-cos(theta)*cos(i)*sin(omega), -sin(theta)*sin(omega)+cos(theta)*cos(i)*cos(omega), cos(theta)*sin(i);
     sin(i)*sin(omega), -sin(i)*cos(omega), cos(i)]