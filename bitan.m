function [delta_V,delta_T] = bitan(a1,e1,a2,e2,mu)
% Questa funzione esegue un trasferimento da un orbita con semiasse maggiore
% a1 ed eccentricità e1 ad un'orbita con semiasse maggiore a2 ed
% eccentrità e2. Viene utilizzato un trasferimento BITANGENTE

h = waitbar(0,'Computating times...');

if nargin == 4
    w = msgbox('Hai dimenticato mu, lo sto automaticamente settando a 398600');
    mu = 398600;
end

