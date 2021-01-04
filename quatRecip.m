function qr=quatRecip(q)

qcon=[q(1);
      -q(2);
      -q(3);
      -q(4)]
  
qnorm=norm(q)

qr=qcon/qnorm^2
