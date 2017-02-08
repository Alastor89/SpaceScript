clear all
close all
clc

% File in cui salverò i dati
filename = 'results.txt';

% Asse dei poli
l = 25000;
Z = zeros(l);
Z(:,3) = linspace(-l/2,l/2,l);

% Importo i dati
load data.mat

% Costruisco i vettori colonna per posizione e velocità
r_i = [r_x_i; r_y_i; r_z_i];
v_i = [v_x_i ; v_y_i ; v_z_i];

% Parametro mu terra
mu_t = 398600;

% Calcolo i parametri orbiati per il punto iniziale
[a_i,e_i,i_i,omega_i,OMEGA_i,theta_i] = pos2par(r_i,v_i,mu_t);
e_i = norm(e_i);

% Calcolo i vettori r,v per il punto finale
[r_f,v_f] = par2pos(a_f,e_f,i_f,OMEGA_f,omega_f,theta_f,mu_t);


% Calcolo posizioni iniziali e finali al variare di theta, da 0 a 360
theta_deg = [1:1:360];
theta_rad = deg2rad(theta_deg);

% Calcolo posizioni iniziali e finali al variare di theta, da 0 a 360
theta_deg_vid = [1:1:720];
theta_rad_vid = deg2rad(theta_deg_vid);



for k = 1:length(theta_deg)
    [r_iniz(:,k),v_iniz(:,k)] = par2pos(a_i,e_i,i_i,OMEGA_i,omega_i,theta_rad(k),mu_t);
    [r_fin(:,k),v_fin(:,k)] = par2pos(a_f,e_f,i_f,OMEGA_f,omega_f,theta_rad(k),mu_t);
end

theta_i_vid = floor(rad2deg(theta_i));

% Cambio di inclinazione

[delta_v_inc_change_A,theta_inc_change_A,omega_inc_change_A]=inclinationchange(a_i,e_i,i_i,OMEGA_i,i_f,OMEGA_f,omega_i,0,mu_t);
[r_inc_change_A,v_inc_change_A] = par2pos(a_i,e_i,i_f,OMEGA_f,omega_inc_change_A,theta_inc_change_A,mu_t);

theta_inc_change_vid = floor(rad2deg(theta_inc_change_A));

% Preallocazione memoria
r_inc_change = zeros(3,length(theta_deg));

% Vettore contenente l'orbita di cambio piano

for k = 1:length(theta_deg)
     [r_inc_change(:,k),v_inc_change(:,k)] = par2pos(a_i,e_i,i_f,OMEGA_f,omega_inc_change_A,theta_rad(k),mu_t);    
end 


% Cambio anomalia di pericentro con calcolo del punto di arrivo e del punto
% di manovra

[delta_v,theta_man,theta_after_man]= anoperichange(a_i, e_i, omega_inc_change_A, omega_f, theta_inc_change_A, mu_t)
[r_per_change_man,v_per_change_man] = par2pos(a_i,e_i,i_f,OMEGA_f,omega_f,theta_man,mu_t);
[r_per_change_after,v_per_change_after] = par2pos(a_i,e_i,i_f,OMEGA_f,omega_f,theta_after_man,mu_t);

theta_man_vid = floor(rad2deg(theta_man));
theta_after_man_vid = (floor(rad2deg(theta_after_man))+360);
%theta_per_anoper = 

% Vettore contenente l'orbita di cambio anomalia di pericentro

for k = 1:length(theta_deg)
    [r_per_change(:,k),v_per_change(:,k)] = par2pos(a_i,e_i,i_f,OMEGA_f,omega_f,theta_rad(k),mu_t);
end


% Trasferimento all'orbita finale

[a_bit , e_bit , delta_v_bit, delta_t ] = bitang(a_i,e_i,a_f,e_f,mu_t);


for k = 1:length(theta_deg)
    [r_bitan(:,k),v_bitan(:,k)] = par2pos(a_bit,e_bit,i_f,OMEGA_f,omega_f,theta_rad(k),mu_t);
end

p_v = 3;

f(1) = figure(1);
% Setto i nomi degli assi
xlabel('x [Km]')
ylabel('y [Km]')
zlabel('z [Km]')
% Apre il file video e lo comincia a catturare

vidObj = VideoWriter('orbit_1.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);

% Plot orbita iniziale

w=1;

for k = theta_i_vid:1:theta_inc_change_vid
    clf
    cla reset
    hold on
    grid on
    axis vis3d
    view(3)
    % Plot complessivo delle orbite
    plot3(r_iniz(1,:),r_iniz(2,:),r_iniz(3,:),'r','linewidth',2)
    plot3(r_inc_change(1,:),r_inc_change(2,:),r_inc_change(3,:),'g-.','linewidth',2);
    plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12)
    plot3(r_per_change(1,:),r_per_change(2,:),r_per_change(3,:),'m-.','linewidth',2);
    plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12);
    plot3(r_per_change_man(1),r_per_change_man(2),r_per_change_man(3),'m*','linewidth',12)
    plot3(r_bitan(1,:),r_bitan(2,:),r_bitan(3,:),'c-.','linewidth',2);
    plot3(r_fin(1,:),r_fin(2,:),r_fin(3,:),'linewidth',2)
    plot3(r_fin(1,end),r_fin(2,end),r_fin(3,end),'*b','linewidth',12)
    % Plot animazione che mi interessa
    [r_iniz_vid(w,:),v_iniz_vid(w,:)]=par2pos(a_i,e_i,i_i,OMEGA_i,omega_i,theta_rad_vid(k),mu_t);
    plot3(r_iniz_vid(w,1),r_iniz_vid(w,2),r_iniz_vid(w,3),'ok','LineWidth',3);             % Pallino
    plot3(r_iniz_vid(1:w,1),r_iniz_vid(1:w,2),r_iniz_vid(1:w,3),'k','LineWidth',3);        % Scia
    eh = earth;
    view(3)
    rotate(eh,[0 0 1],3*w,[0 0 0]);
    writeVideo(vidObj, getframe(gca));
    w = w+1;
end

clear w
close(vidObj);
close(f(1));

% Plot plot cambio di piano
f(2) = figure(2);
% Setto i nomi degli assi
xlabel('x [Km]')
ylabel('y [Km]')
zlabel('z [Km]')
vidObj = VideoWriter('orbit_2.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
w=1;

for k = 236:1:theta_after_man_vid   
    clf
    cla reset
    hold on
    axis vis3d
    grid on
    view(3)
    % Plot complessivo delle orbite
    plot3(r_iniz(1,:),r_iniz(2,:),r_iniz(3,:),'r','linewidth',2)
    plot3(r_iniz(1,151:234),r_iniz(2,151:234),r_iniz(3,151:234),'k','linewidth',2)
    plot3(r_inc_change(1,:),r_inc_change(2,:),r_inc_change(3,:),'g-.','linewidth',2);
    plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12)
    plot3(r_per_change(1,:),r_per_change(2,:),r_per_change(3,:),'m-.','linewidth',2);
    plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12);
    plot3(r_per_change_man(1),r_per_change_man(2),r_per_change_man(3),'m*','linewidth',12)
    plot3(r_bitan(1,:),r_bitan(2,:),r_bitan(3,:),'c-.','linewidth',2);
    plot3(r_fin(1,:),r_fin(2,:),r_fin(3,:),'linewidth',2)
    plot3(r_fin(1,end),r_fin(2,end),r_fin(3,end),'*b','linewidth',12)
    % Plot Animazione che mi interessa
    [r_inc_change_vid(w,:),v_inc_change_vid(w,:)]=par2pos(a_i,e_i,i_f,OMEGA_f,omega_inc_change_A,theta_rad_vid(k),mu_t);
    plot3(r_inc_change_vid(w,1),r_inc_change_vid(w,2),r_inc_change_vid(w,3),'ok','LineWidth',3);                   %Pallino
    plot3(r_inc_change_vid(1:w,1),r_inc_change_vid(1:w,2),r_inc_change_vid(1:w,3),'k','LineWidth',3);  % Scia
    eh = earth;
    rotate(eh,[0 0 1],3*w,[0 0 0]);
    view(3)
    writeVideo(vidObj, getframe(gca));
    w = w+1;
end

clear w
close(vidObj);
close(f(2));

% Plot arrivo all'apocentro orbita magenta

f(3) = figure(3);
% Setto i nomi degli assi
xlabel('x [Km]')
ylabel('y [Km]')
zlabel('z [Km]')
vidObj = VideoWriter('orbit_3.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
w=1;

for k = 262:1:360
    clf
    cla reset
    hold on
    axis vis3d
    grid on
    view(3)
    % Plot complessivo delle orbite
    plot3(r_iniz(1,:),r_iniz(2,:),r_iniz(3,:),'r','linewidth',2)
    plot3(r_iniz(1,151:234),r_iniz(2,151:234),r_iniz(3,151:234),'k','linewidth',2)
    plot3(r_inc_change(1,:),r_inc_change(2,:),r_inc_change(3,:),'g-.','linewidth',2);
    plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12)
    plot3(r_per_change(1,:),r_per_change(2,:),r_per_change(3,:),'m-.','linewidth',2);
    plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12);
    
    plot3(r_inc_change(1,236:end),r_inc_change(2,236:end),r_inc_change(3,236:end),'k','linewidth',2)
    plot3(r_inc_change(1,1:103),r_inc_change(2,1:103),r_inc_change(3,1:103),'k','linewidth',2);
    
    plot3(r_per_change_man(1),r_per_change_man(2),r_per_change_man(3),'m*','linewidth',12)
    plot3(r_bitan(1,:),r_bitan(2,:),r_bitan(3,:),'c-.','linewidth',2);
    plot3(r_fin(1,:),r_fin(2,:),r_fin(3,:),'linewidth',2)
    plot3(r_fin(1,end),r_fin(2,end),r_fin(3,end),'*b','linewidth',12)
    % Plot Animazione che mi interessa
    [r_per_change_vid(w,:),v_per_change_vid(w,:)]=par2pos(a_i,e_i,i_f,OMEGA_f,omega_f,theta_rad_vid(k),mu_t);
    plot3(r_per_change_vid(w,1),r_per_change_vid(w,2),r_per_change_vid(w,3),'ok','LineWidth',3);                   %Pallino
    plot3(r_per_change_vid(1:w,1),r_per_change_vid(1:w,2),r_per_change_vid(1:w,3),'k','LineWidth',3);  % Scia
    eh = earth;
    rotate(eh,[0 0 1],3*w,[0 0 0]);
    view(3)
    writeVideo(vidObj, getframe(gca));
    w = w+1;
end

clear w
close(vidObj);
close(f(3));

% Quarto plot
f(4) = figure(4);
% Setto i nomi degli assi
xlabel('x [Km]')
ylabel('y [Km]')
zlabel('z [Km]')
vidObj = VideoWriter('orbit_4.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 15;
open(vidObj);
w=1;

for k = 1:1:180
    clf
    cla reset
    hold on
    axis vis3d
    grid on
    view(3)
    % Plot complessivo delle orbite
    plot3(r_iniz(1,:),r_iniz(2,:),r_iniz(3,:),'r','linewidth',2)
    plot3(r_iniz(1,151:234),r_iniz(2,151:234),r_iniz(3,151:234),'k','linewidth',2)
    plot3(r_inc_change(1,:),r_inc_change(2,:),r_inc_change(3,:),'g-.','linewidth',2);
    plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12)
    plot3(r_per_change(1,:),r_per_change(2,:),r_per_change(3,:),'m-.','linewidth',2);
    plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12);
    

    plot3(r_inc_change(1,236:end),r_inc_change(2,236:end),r_inc_change(3,236:end),'k','linewidth',2)
    plot3(r_inc_change(1,1:103),r_inc_change(2,1:103),r_inc_change(3,1:103),'k','linewidth',2);
    plot3(r_inc_change(1,236:end),r_inc_change(2,236:end),r_inc_change(3,236:end),'k','linewidth',2);
    plot3(r_inc_change(1,1:103),r_inc_change(2,1:103),r_inc_change(3,1:103),'k','linewidth',2);
    
    plot3(r_per_change_man(1),r_per_change_man(2),r_per_change_man(3),'m*','linewidth',12)
    
    plot3(r_per_change(1,262:360),r_per_change(2,262:360),r_per_change(3,262:360),'k','linewidth',2);
    
    
    plot3(r_bitan(1,:),r_bitan(2,:),r_bitan(3,:),'c-.','linewidth',2);
    plot3(r_fin(1,:),r_fin(2,:),r_fin(3,:),'linewidth',2)
    plot3(r_fin(1,end),r_fin(2,end),r_fin(3,end),'*b','linewidth',12)
    % Plot Animazione che mi interessa
    [r_bitan_vid(w,:),v_bitan_vid(w,:)]=par2pos(a_bit,e_bit,i_f,OMEGA_f,omega_f,theta_rad_vid(k),mu_t);
    plot3(r_bitan_vid(w,1),r_bitan_vid(w,2),r_bitan_vid(w,3),'ok','LineWidth',3);                   %Pallino
    plot3(r_bitan_vid(1:w,1),r_bitan_vid(1:w,2),r_bitan_vid(1:w,3),'k','LineWidth',3);  % Scia
    eh = earth;
    rotate(eh,[0 0 1],3*w,[0 0 0]);
    view(3)
    writeVideo(vidObj, getframe(gca));
    w = w+1;
end

clear w
close(vidObj);
close(f(4));

% Plot finale

f(5) = figure(5);
% Setto i nomi degli assi
xlabel('x [Km]')
ylabel('y [Km]')
zlabel('z [Km]')
vidObj = VideoWriter('orbit_5.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 15;
open(vidObj);
w=1;

for k = 180:1:360
    clf
    cla reset
    hold on
    axis vis3d
    grid on
    view(3)
    % Plot complessivo delle orbite
    plot3(r_iniz(1,:),r_iniz(2,:),r_iniz(3,:),'r','linewidth',2)
    plot3(r_iniz(1,151:234),r_iniz(2,151:234),r_iniz(3,151:234),'k','linewidth',2)
    plot3(r_inc_change(1,:),r_inc_change(2,:),r_inc_change(3,:),'g-.','linewidth',2);
    plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12)
    plot3(r_per_change(1,:),r_per_change(2,:),r_per_change(3,:),'m-.','linewidth',2);
    plot3(r_inc_change_A(1),r_inc_change_A(2),r_inc_change(3),'g*','linewidth',12);
    plot3(r_per_change_man(1),r_per_change_man(2),r_per_change_man(3),'m*','linewidth',12)
    plot3(r_bitan(1,:),r_bitan(2,:),r_bitan(3,:),'c-.','linewidth',2);
    
    
    
    
    
    plot3(r_per_change(1,262:360),r_per_change(2,262:360),r_per_change(3,262:360),'k','linewidth',2);
    
    
    
    plot3(r_bitan(1,1:180),r_bitan(2,1:180),r_bitan(3,1:180),'k','linewidth',2);
    
    plot3(r_fin(1,:),r_fin(2,:),r_fin(3,:),'linewidth',2)
    plot3(r_fin(1,end),r_fin(2,end),r_fin(3,end),'*b','linewidth',12)
    % Plot Animazione che mi interessa
    [r_fin_vid(w,:),v_fin_vid(w,:)]=par2pos(a_f,e_f,i_f,OMEGA_f,omega_f,theta_rad_vid(k),mu_t);
    plot3(r_fin_vid(w,1),r_fin_vid(w,2),r_fin_vid(w,3),'ok','LineWidth',3);                   %Pallino
    plot3(r_fin_vid(1:w,1),r_fin_vid(1:w,2),r_fin_vid(1:w,3),'k','LineWidth',3);  % Scia
    eh = earth;
    rotate(eh,[0 0 1],3*w,[0 0 0]);
    view(3)
    writeVideo(vidObj, getframe(gca));
    w = w+1;
end

close(vidObj);
close(f(5));

