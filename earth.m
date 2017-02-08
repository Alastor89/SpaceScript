function [eh] = earth ()
% Questa function devo ancora finire di aggiustarla, plotta la terra o la
% terra rotante, a seconda del valore immesso nello switch. Il problema è
% che la rotazione non funziona al di fuori di un ciclo. Il mio consiglio è
% quello di chiamare earth(0) facendo ciclare l'animazione all'interno del
% main-script !

RT = 6378;
L = linspace(-RT-2500,RT+2500,1000);
N = length(L);
a = zeros(N,1); % x-component
b = zeros(N,1); % y-component
c = L'; % z-component, l'asse
load('topo.mat','topo','topomap1');
[x,y,z] = sphere(100);
x = x*RT;
y = y*RT;
z = z*RT;
props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;
hold on
axis equal
eh=surface(x,y,z,props);
light('position',[-1 0 1]);
light('position',[-1.5 0.5 -0.5], 'color', [.6 .2 .2]);
view(3);
plot3(a,b,c,'k--','LineWidth',2);
%hold off

end

     

