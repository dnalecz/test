clear // Brilouin-magnet V6 04-02-15 refactor
        // E:\PROGRAMY-Dawid\ 3/06.16 dodaje external field Bex
function M=f(m,T,Bex)
J=7/2;
N=1;      g=2;
kB=1;      uB=1;
t=(2*J+1)/(2*J);
Lam=0.2619;// K/UB2 = 49 dla SMO daje 5*49=245K(TN)
x=(g*uB*J*(Lam*m+Bex*0.67171))/(kB*T)
BJ=t*coth(t*x)-(1)/(2*J)*coth(x/(2*J))
M=N*g*J*BJ-m
endfunction

iter=850
Lam=0.2619;
Bex=0.1 //[Tesla]

for x=1:iter
clc;
T=x/100.0
x0=7;
xsol=fsolve(x0,f);
//plot(T,xsol,"+g")
temp(x)=T;
mom(x)=xsol;
end
w=[temp,mom];

for x=1:iter
w(x,3)=w(x,2)*w(x,2);
end

for x=1:iter-1
w(x,4)=-(w(x+1,3)-w(x,3))/(w(x+1,1)-w(x,1))*(Lam/2 *8.31); //pochodna * Lam/2 dlaczego nie wiem                           // nie dziala w poluzew
w(x,5)=w(x,4)/w(x,1); // c/T ??                     // niezalezne od g zawsze 2
end

for x=1:iter-2
w(x,6)=(w(x+1,1)-w(x,1))*w(x,5);
w(x,7)=sum(w(:,6)); 
end

// w1 temp      // w2 moment
// w3 m^2       // w4 dm2/dT * Lam
// w5 c/T       // w7 entropy

//clf() // czysci ekran
subplot(2,2,1); 
title('m(T) [uB]'); 
plot(w(:,1),w(:,2),"+b")
subplot(2,2,2); 
title('dm2/dT * Lam [R]'); 
plot(w(:,1),w(:,4),"+b")
subplot(2,2,3); 
title('dm2/dT * Lam/T [R/K]'); 
plot(w(:,1),w(:,5),"+b")
subplot(2,2,4); 
title('Entropy [R]'); 
plot(w(:,1),w(:,7),"+b")
w(:,7)
log(8)






