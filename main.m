clear all
close all
clc

% File in cui salverò i dati
filename = 'results.txt';

% Asse dei poli
l = 1000;
Z = zeros(l);
Z(:,3) = linspace(-l/2,l/2,l);

% Importo i dati
load data.mat

% Costruisco i vettori colonna per posizione e velocità
r_i = [r_x_i; r_y_i; r_z_i];
v_i = [v_x_i ; v_y_i ; v_z_i];

% Calcolola norma di r_i e v_i
r_iniziale = norm(r_i);
v_iniziale = norm(v_i);
fileID = fopen(filename,'w+');
fprintf(fileID,'Norma di r_iniziale: %f \n Norma di v_iniziale: %f\n',r_iniziale,v_iniziale);
fclose(fileID);


% Parametro mu terra
mu_t = 398600;

% Calcolo i parametri orbiati per il punto iniziale
[a_i,e_i,i_i,omega_i,OMEGA_i,theta_i] = pos2par(r_i,v_i,mu_t);
e_i = norm(e_i);

% Calcolo i vettori r,v per il punto finale
[r_f,v_f] = par2pos(a_f,e_f,i_f,OMEGA_f,omega_f,theta_f,mu_t);

% Calcolola norma di r_f e v_f
r_finale = norm(r_f);
v_finale = norm(v_f);
fileID = fopen(filename,'a+');
fprintf(fileID,'Norma di r_finale: %f \n Norma di v_finale: %f\n',r_finale,v_finale);
fclose(fileID);

% Calcolo del periodo delle due orbite
T_iniziale = (2*pi*a_i^(3/2))/sqrt(mu_t);
T_finale = (2*pi*a_f^(3/2))/sqrt(mu_t)
fileID = fopen(filename,'a+');
fprintf(fileID,'Periodo orbita iniziale: %f  ore \n Periodo orbita finale: %f ore \n',T_iniziale/3600,T_finale/3600);
fclose(fileID);


% Calcolo dei parametri orbitali in gradi per le nostre orbite
i_i_deg = rad2deg(i_i);
OMEGA_i_deg = rad2deg(OMEGA_i);
omega_i_deg = rad2deg(omega_i);
i_f_deg = rad2deg(i_f);
OMEGA_f_deg = rad2deg(OMEGA_f);
omega_f_deg = rad2deg(omega_f);
fileID = fopen(filename,'a+');
fprintf(fileID,'i_i: %f\n OMEGA_i: %f\n omega_i: %f\n i_f: %f\n OMEGA_f: %f\n omega_f: %f\n  \n',i_i_deg,OMEGA_i_deg,omega_i_deg,i_f_deg,OMEGA_f_deg,omega_f_deg);
fclose(fileID);

% Calcolo momento della quantita di monto
h_i = cross(r_i,v_i);
h_f = cross(r_f,v_f);

% Calcolo posizioni iniziali e finali al variare di theta, da 0 a 360
theta_deg = [1:1:360];
theta_rad = deg2rad(theta_deg);

% Preallocazione memoria
r_iniz = zeros(3,length(theta_deg));
r_fin = zeros(3,length(theta_deg));
v_iniz = zeros(3,length(theta_deg));
v_fin = zeros(3,length(theta_deg));

for k = 1:length(theta_deg)
    [r_iniz(:,k),v_iniz(:,k)] = par2pos(a_i,e_i,i_i,OMEGA_i,omega_i,theta_rad(k),mu_t);
    [r_fin(:,k),v_fin(:,k)] = par2pos(a_f,e_f,i_f,OMEGA_f,omega_f,theta_rad(k),mu_t);    
end
    
% Plot orbita iniziale
f(1) = figure(1); 
hold on
view(3)
plot3(r_iniz(1,:),r_iniz(2,:),r_iniz(3,:),'r','linewidth',3.5)
% plot3(r_iniz(1,180),r_iniz(2,180),r_iniz(3,180),'k*','linewidth',5)
% plot3(r_iniz(1,end),r_iniz(2,end),r_iniz(3,end),'k*','linewidth',5)
grid on
zoom on
axis vis3d
earth
hleg1 = legend('Initial Orbit');
xlabel('x [Km]')
ylabel('y [Km]')
zlabel('z [Km]')

% Plot orbita finale
f(2) = figure(2);
hold on
plot3(r_iniz(1,:),r_iniz(2,:),r_iniz(3,:),'r','linewidth',3.5)
plot3(r_fin(1,:),r_fin(2,:),r_fin(3,:),'linewidth',3.5)
% plot3(r_iniz(1,180),r_iniz(2,180),r_iniz(3,180),'k*','linewidth',5)
% plot3(r_iniz(1,end),r_iniz(2,end),r_iniz(3,end),'k*','linewidth',5)
% plot3(r_fin(1,180),r_fin(2,180),r_fin(3,180),'k*','linewidth',5)
% plot3(r_fin(1,end),r_fin(2,end),r_fin(3,end),'k*','linewidth',5)
grid on
zoom on
axis vis3d
hleg1 = legend('Initial Orbit','Final Orbit');
earth;
xlabel('x [Km]')
ylabel('y [Km]')
zlabel('z [Km]')

[delta_v_inc_change_A,theta_inc_change_A,omega_inc_change_A]=inclinationchange(a_i,e_i,i_i,OMEGA_i,i_f,OMEGA_f,omega_i,0,mu_t);
[r_inc_change_A,v_inc_change_A] = par2pos(a_i,e_i,i_f,OMEGA_f,omega_inc_change_A,theta_inc_change_A,mu_t)

fileID = fopen(filename,'a+');
fprintf(fileID,'Dati cambio di piano punto A \n Delta_V = %f \n Anomalia Vera = %f \n Anomalia di pericentro = %f \n',delta_v_inc_change_A,theta_inc_change_A,omega_inc_change_A);
fclose(fileID);

[delta_v_inc_change_B,theta_inc_change_B,omega_inc_change_B]=inclinationchange(a_i,e_i,i_i,OMEGA_i,i_f,OMEGA_f,omega_i,1,mu_t);
[r_inc_change_B,v_inc_change_B] = par2pos(a_i,e_i,i_f,OMEGA_f,omega_inc_change_B,theta_inc_change_B,mu_t)

fileID = fopen(filename,'a+');
fprintf(fileID,'Dati cambio di piano punto B \n Delta_V = %f \n Anomalia Vera = %f \n Anomalia di pericentro = %f \n',delta_v_inc_change_B,theta_inc_change_B,omega_inc_change_B);
fclose(fileID);


% Preallocazione memoria
r_inc_change = zeros(3,length(theta_deg));
v_inc_change = zeros(3,length(theta_deg));

for k = 1:length(theta_deg)
     [r_inc_change(:,k),v_inc_change(:,k)] = par2pos(a_i,e_i,i_f,OMEGA_f,omega_inc_change_A,theta_rad(k),mu_t);    
end 

f(3) = figure(3);
% Plot orbita di trasfermineto e punto in cui effettuo la manovra
hold on
plot3(r_iniz(1,:),r_iniz(2,:),r_iniz(3,:),'r','linewidth',3.5)
plot3(r_inc_change(1,:),r_inc_change(2,:),r_inc_change(3,:),'g','linewidth',2);
plot3(r_iniz(1,234),r_iniz(2,234),r_iniz(3,234),'g*','linewidth',12)
plot3(r_iniz(1,54),r_iniz(2,54),r_iniz(3,54),'g*','linewidth',12)
%plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12)
grid on
zoom on
axis vis3d
earth
hleg2 = legend('Initial Orbit','Inclination Change Orbit','Operation point 1', 'Operation point 2');
hold on

% Calcolo dei tempi
[delta_t_inc_change_A] = time(a_i,e_i,theta_i,theta_inc_change_A,mu_t)
[delta_t_inc_change_B] = time(a_i,e_i,theta_i,theta_inc_change_B,mu_t)

fileID = fopen(filename, 'a+');
fprintf(fileID,'Tempi \nTempo impiegato ad arrivare al punto A : %f\nTempo impiegato ad arrivare al punto B : %f\n',delta_t_inc_change_A/3600,delta_t_inc_change_B/3600);
fclose(fileID);

% Cambio anomalia di pericentro

%[delta_v, theta_man, theta_after_man] = perianomaly_change(a_i, e_i, omega_inc_change_A, omega_f, theta_inc_change_A, mu_t)
[delta_v,theta_man,theta_after_man]= anoperichange(a_i, e_i, omega_inc_change_A, omega_f, theta_inc_change_A, mu_t)
[r_per_change_man,v_per_change_man] = par2pos(a_i,e_i,i_f,OMEGA_f,omega_f,theta_man,mu_t); 
[r_per_change_after,v_per_change_after] = par2pos(a_i,e_i,i_f,OMEGA_f,omega_f,theta_after_man,mu_t); 
[delta_t_inc_change_A_AFTER] = time(a_i,e_i,theta_inc_change_A,theta_man,mu_t)

delta_t_inc_change_A_AFTER/3600

for k = 1:length(theta_deg)
    [r_per_change(:,k),v_per_change(:,k)] = par2pos(a_i,e_i,i_f,OMEGA_f,omega_f,theta_rad(k),mu_t);    
end 

f(4) = figure(4);
% Plot orbita di trasferumentoo e punto in cui effettuo la manovra
hold on
plot3(r_inc_change(1,:),r_inc_change(2,:),r_inc_change(3,:),'g','linewidth',2);
plot3(r_per_change(1,:),r_per_change(2,:),r_per_change(3,:),'m','linewidth',2);
plot3(r_per_change_man(1),r_per_change_man(2),r_per_change_man(3),'m*','linewidth',12)
plot3(r_per_change(1,77),r_per_change(2,77),r_per_change(3,77),'m*','linewidth',12);
grid on
zoom on
axis vis3d
earth
view(2)
hleg3 = legend('Inclination Change Orbit','Argument of Pericenter Change Orbit', 'Operation Point 1', 'Operation point 2' );
hold on

% Trasferimento BITANGENTE 

[dvA,dtA,e_bit_A,a_bit_A,dvB,dtB,e_bit_B,a_bit_B] = bitangente (a_i,e_i,a_f,e_f,mu_t)

for k = 1:length(theta_deg)
     [r_bitan_A(:,k),v_bitan_B(:,k)] = par2pos(a_bit_A,e_bit_A,i_f,OMEGA_f,omega_f,theta_rad(k),mu_t); 
     [r_bitan_B(:,k),v_bitan_B(:,k)] = par2pos(a_bit_B,e_bit_B,i_f,OMEGA_f,omega_f,theta_rad(k),mu_t); 
end 


f(5) = figure(5);
% Plot orbita di trasferumentoo e punto in cui effettuo la manovra
hold on
plot3(r_per_change(1,:),r_per_change(2,:),r_per_change(3,:),'m','linewidth',2);
plot3(r_bitan_A(1,1:180),r_bitan_A(2,1:180),r_bitan_A(3,1:180),'c','linewidth',2);
plot3(r_bitan_B(1,180:end),-r_bitan_B(2,180:end),r_bitan_B(3,180:end),'c','linewidth',2);
plot3(r_fin(1,:),r_fin(2,:),r_fin(3,:),'linewidth',2)
grid on
zoom on
axis vis3d
earth
view(2)
hleg4 = legend('Argument of Pericenter Change Orbit', 'Bitangent Transfer 1', 'Bitangent Trasfer 2', 'Final Orbit');
hold on

% Orgia finale 

% f(6) = figure(6);
% hold on
% plot3(r_iniz(1,:),r_iniz(2,:),r_iniz(3,:),'r','linewidth',3.5)
% plot3(r_inc_change(1,:),r_inc_change(2,:),r_inc_change(3,:),'g-.','linewidth',2);
% plot3(r_per_change(1,:),r_per_change(2,:),r_per_change(3,:),'m-.','linewidth',2);
% plot3(r_bitan_A(1,1:180),r_bitan_A(2,1:180),r_bitan_A(3,1:180),'c','linewidth',2);
% plot3(r_fin(1,:),r_fin(2,:),r_fin(3,:),'linewidth',3.5)
% % Plot del percorso
% plot3(r_per_change(1,258:end),r_per_change(2,258:end),r_per_change(3,258:end),'m','linewidth',2);
% plot3(r_inc_change(1,236:end),r_inc_change(2,236:end),r_inc_change(3,236:end),'g','linewidth',2);
% plot3(r_inc_change(1,236:-1:360),r_inc_change(2,236:-1:360),r_inc_change(3,236:-1:360),'g','linewidth',2);
% plot3(r_inc_change(1,1:103),r_inc_change(2,1:103),r_inc_change(3,1:103),'g','linewidth',2);
% % Fine plot percorso
% plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12)
% plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12);
% plot3(r_per_change_man(1),r_per_change_man(2),r_per_change_man(3),'m*','linewidth',12)
% plot3(r_bitan_A(1,1),r_bitan_A(2,1),r_bitan_A(3,1),'c*','linewidth',12);
% plot3(r_bitan_A(1,180),r_bitan_A(2,180),r_bitan_A(3,180),'b*','linewidth',12);
% grid on
% zoom on
% axis vis3d
% earth
% hleg4 = legend('Initial Orbit','Inclination Change Orbit','Argument of Pericenter Change Orbit', 'Bitangent Transfer', 'Final Orbit');
% hold on

% Orgia finale 

f(6) = figure(6);
hold on
plot3(r_iniz(1,:),r_iniz(2,:),r_iniz(3,:),'r','linewidth',3.5)
plot3(r_inc_change(1,:),r_inc_change(2,:),r_inc_change(3,:),'g-.','linewidth',2);
plot3(r_per_change(1,:),r_per_change(2,:),r_per_change(3,:),'m-.','linewidth',2);
plot3(r_bitan_A(1,180:end),r_bitan_A(2,180:end),r_bitan_A(3,180:end),'c-.','linewidth',2);
plot3(r_fin(1,:),r_fin(2,:),r_fin(3,:),'linewidth',3.5)
% Plot del percorso
plot3(r_bitan_A(1,1:180),r_bitan_A(2,1:180),r_bitan_A(3,1:180),'c','linewidth',2);

%plot3(r_per_change(1,258:end),r_per_change(2,258:end),r_per_change(3,258:end),'m','linewidth',2);

plot3(r_per_change(1,85:end),r_per_change(2,85:end),r_per_change(3,85:end),'m','linewidth',2);


%plot3(r_inc_change(1,236:end),r_inc_change(2,236:end),r_inc_change(3,236:end),'g','linewidth',2);

plot3(r_inc_change(1,236:285),r_inc_change(2,236:285),r_inc_change(3,236:285),'g','linewidth',2);

%plot3(r_inc_change(285),r_inc_change(2,285),r_inc_change(3,285),'g*','linewidth',12);

%plot3(r_inc_change(1,236:-1:360),r_inc_change(2,236:-1:360),r_inc_change(3,236:-1:360),'g','linewidth',2);
%plot3(r_inc_change(1,1:103),r_inc_change(2,1:103),r_inc_change(3,1:103),'g','linewidth',2);
% Fine plot percorso
plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12)
%plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12);
%plot3(r_per_change_man(1,180),r_per_change_man(2,180),r_per_change_man(3,180),'m*','linewidth',12)

plot3(r_per_change(1,80),r_per_change(2,80),r_per_change(3,80),'m*','linewidth',12)


plot3(r_bitan_A(1,1),r_bitan_A(2,1),r_bitan_A(3,1),'c*','linewidth',12);
plot3(r_bitan_A(1,180),r_bitan_A(2,180),r_bitan_A(3,180),'c*','linewidth',12);
grid on
zoom on
axis vis3d
earth



