function [dvA,dtA,etA,atA,dvB,dtB,etB,atB] = Btg (a1,e1,a2,e2,mhu)

p1=a1*(1-e1^2);
p2=a2*(1-e2^2);
rp1=a1*(1-e1);
rp2=a2*(1-e2);
ra1=a1*(1+e1);
ra2=a2*(1+e2);

%% Caso A
% parametri dell'orbita di trasferimento
atA=(rp1+ra2)/2;
etA=(ra2-rp1)/(ra2+rp1);
ptA=atA*(1-etA^2);

% velocit� di partenza e arrivo: prima sono al pericentro poi all'apocentro
v1A=sqrt(mhu/p1)*(1+e1);
v2A=sqrt(mhu/p2)*(1-e2);

% velocit� sull'orbita di trasferimento
vt1A=sqrt(mhu/ptA)*(1+etA);
vt2A=sqrt(mhu/ptA)*(1-etA);

% dv e dt caso A
dvA=abs(vt1A-v1A)+abs(v2A-vt2A);
dtA=(pi*atA^1.5)/sqrt(mhu);

%% Caso B
% parametri dell'orbita di trasferimento
atB=(ra1+rp2)/2;
etB=abs(ra1-rp2)/(ra1+rp2);
ptB=atB*(1-etB^2);

% velocit� di partenza e arrivo: prima sono al pericentro poi all'apocentro
v1B=sqrt(mhu/p1)*(1-e1);
v2B=sqrt(mhu/p2)*(1+e2);

% velocit� sull'orbita di trasferimento con controllo
if ra1>rp2
    disp('l''apocentro della prima � pi� grande del pericentro della seconda')
    vt1B=sqrt(mhu/ptB)*(1-etB);
    vt2B=sqrt(mhu/ptB)*(1+etB);
else
    disp('il pericentro della seconda � pi� grande dell''apocentro della prima, caso invertito')
    vt1B=sqrt(mhu/ptB)*(1+etB);
    vt2B=sqrt(mhu/ptB)*(1-etB);
end


% dv e dt caso A
dvB=abs(vt1B-v1B)+abs(v2B-vt2B);
dtB=(pi*atB^1.5)/sqrt(mhu);