function der = integrationTrackGeoSat(t,x,A,B,Jc,w_c,q0e_0)
w=x(1:3)
q0=x(4)
q=x(5:7)
q0c=x(8)
qc=x(9:11)



wc_tilde = [0, -w_c(3), w_c(2);
             w_c(3), 0, -w_c(1);
              -w_c(2), w_c(1), 0]
          
w_tilde = [0, -w(3), w(2);
             w(3), 0, -w(1);
              -w(2), w(1), 0]
we=w-w_c

qc_tilde = [0, -qc(3), qc(2);
            qc(3), 0, -qc(1);
            -qc(2), qc(1), 0]

qe=[q0c, transpose(qc); 
    -qc, q0c*eye(3)-qc_tilde]*[q0;q]

Tc=w_tilde*Jc*w-Jc*(A\B)*we-sign(q0e_0)*Jc*(A\qe(2:4))



der(1:3,1)=Jc\(-w_tilde*Jc*w+Tc)
der(4,1)=-0.5*transpose(w)*q
der(5:7,1)=-0.5*w_tilde*q+0.5*q0*w
der(8,1)=-0.5*transpose(w_c)*qc
der(9:11,1)=-0.5*wc_tilde*qc+0.5*q0c*w_c

