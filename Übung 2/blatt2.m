
kmax=7;  % Hier sollten Sie bei einem ersten Test zun�chst mit einem k< 7 starten

% Vektoren zeit und fehler f�r die Teilaufgabe 3
zeit=zeros(kmax,1); fehler=zeros(kmax,1);  

    k=2;  %  hier nur ein kleines k als Beispiel f�r die erste Matrix A
    m=2^k;  n=m^2;
    
    % Blockmatrizen bauen mit Hilfe von Kroneckerprodukten:
    % sparse Matrizen  (Tipp: mit full(B) k�nnen Sie wieder
    % die volle Matrix erhalten um Sie zu pr�fen/visualisieren)
    e = ones(m,1);
    B = spdiags([-e 4*e -e], -1:1, m, m);   % Bandmatrix B
    E=speye(m);                             % Identit�t                     
    D=spdiags([-e -e], [-1,1],m,m);
    A=kron(E,B) + kron(D,E);                % Matrix A als Kroneckerprodukt

    
    % Nun k�nnen Sie starten mit Ihrer LR-Zerlegung
    % Hilfreich k�nnen folgende Befehler sein:
    % tic  und toc f�r die Zeitmessung 
    % tril und triu  f�r Dreiecksmatrizen
    
    



%% Graphische Ausgabe: Hier k�nnen Sie Ihr Ergebnis zeichen und speichern
% Speichen  des Fenters 1 k�nnen Sie zum Beispiel mit Hilfe von
% saveas(1, 'ergebnis_blatt2', 'png')


    
%% Funktion LR zur Berechnung der LR-Zerleung der Matrix A
function [A] = LR(A)



end