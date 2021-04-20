sysm alpha beta delta
g=[0 delta^(-1) -beta*(delta^(-2))];
syms a11 a12 a13 a21 a22 a23 a31 a32 a33
O=[a11 a12 a13; a21 a22 a23; a31 a32 a33];
gt=[0; delta^(-1); -beta*(delta^(-2))];
g*O*gt