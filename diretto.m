function[a,e,i,ohm,omega,teta]=diretto(vett_r,vett_v,mu)

% Dati i vettori r,v e il parametro costante mu individua
% (a,e,i,asc,omega,teta), con gli angoli espressi in radianti

I=[1 0 0];
J=[0 1 0];
K=[0 0 1];

r=norm(vett_r);
v=norm(vett_v);

a=1/((2/r)-(v^2/mu));

vett_h=cross(vett_r,vett_v);
h=norm(vett_h);

vett_e=(1/mu)*(cross(vett_v,vett_h)-mu*vett_r/r);
e=norm(vett_e);

i=acos(dot(vett_h,K)/h);

vett_n=cross(K,vett_h)/norm(cross(K,vett_h));

if dot(vett_n,J)>0
    ohm=acos(dot(I,vett_n));
else
    ohm=2*pi-acos(dot(I,vett_n));
end

if dot(vett_e,K)>0
    omega=acos(dot(vett_n,vett_e)/e); 
else
    omega=2*pi-acos(dot(vett_n,vett_e)/e);
end

v_rad=dot(vett_v,vett_r);
if v_rad>0
    fprintf('Maggiore')
    teta=acos(dot(vett_r,vett_e)/(r*e));
else
    fprintf('Minore')
    teta=2*pi-acos(dot(vett_r,vett_e)/(r*e));
end

end