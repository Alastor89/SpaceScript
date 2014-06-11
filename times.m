function [delta_t] = times(a,e,theta_1,theta_2,mu)
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
% Questa function implementa il cambio di anomalia del pericentro tra due
% orbite aventi lo stesso semiasse maggiore e la stessa eccentricità. 
% Fornire il parametro gravitazionale di riferimento
% 
% Questa function serve a calcolare il tempo trascorso per spostarsi da un
% anomalia vera theta_1 ad un anomalia vera theta_2, dati semiasse maggiore
% ed eccentricità 



% Calcolo l'anomalia eccentrica    http://en.wikipedia.org/wiki/Eccentric_anomaly

E_1 = 2*atan(sqrt((1-e)/(1+e)))*tan(theta_1/2);
E_2 = 2*atan(sqrt((1-e)/(1+e)))*tan(theta_2/2);

% Calcolo i tempi dal pericentro, la formula l'abbiamo vista all'inizio del
% corso !!!

t_1 = sqrt(a^3/mu)*(E_1 - e*sin(E_1));
t_2 = sqrt(a^3/mu)*(E_2 - e*sin(E_2));

% Posso quindi calcolarmi il periodo che mi servirà più avanti 

T =  2*pi/sqrt(mu) * a^(3/2);

% Risolvo la questione del tempo negativo aggiungendoci il periodo al tempo

if t_1 < 0
    t_1 = t_1 + T;
end

if t_2 < 0
    t_2 = t_2 + T;
end

% A questo punto posso finalmente calcolare i tempi di percorrenza

if theta_2 > theta_1
    delta_t = t_2 - t_1;
else
    delta_t = t_2 - t_1 + T;
end



