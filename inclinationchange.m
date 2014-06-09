function [delta_v,theta_1,omega_2]=inclinationchange(a,e,i_1,OMEGA_1,i_2,OMEGA_2,omega_1,mu)
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
% 
% Questa function implementa il cambio di piano, restituisce la variazione di i e OMEGA
%
% PARAMETRI DI USCITA :
% 
% delta_v                         Costo della manovra in termini di delta_v
% theta_1			  Anomalia vera all'inizio della manovra
% omega_2			  Anomalia di pericentro alla fine della manovra
%
% % Author : Francescodario Cuzzocrea
% emal : francescodario.cuzzocrea@mail.polimi.it
% (C) 2014


% Calcolo i delta_i e delta_OMEGA 

delta_i = i_2 - i_1;
delta_OMEGA = OMEGA_2 - OMEGA_1;
delta_check = delta_i*delta_omega

% Calcolo alpha sfruttando il teorma del coseno visto a lezione


