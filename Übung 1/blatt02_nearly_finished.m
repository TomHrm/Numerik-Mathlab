kmax=7;  % Hier sollten Sie bei einem ersten Test zunächst mit einem k< 7 starten

% Vektoren zeit und fehler für die Teilaufgabe 3
zeit=zeros(kmax,1); fehler=zeros(kmax,1);  

    k=2;  %  hier nur ein kleines k als Beispiel für die erste Matrix A
    m=2^k;  n=m^2;
    
    % Blockmatrizen bauen mit Hilfe von Kroneckerprodukten:
    % sparse Matrizen  (Tipp: mit full(B) können Sie wieder
    % die volle Matrix erhalten um Sie zu prüfen/visualisieren)
    e = ones(m,1);
    B = spdiags([-e 4*e -e], -1:1, m, m);   % Bandmatrix B
    E = speye(m);                             % Identität                     
    D = spdiags([-e -e], [-1,1],m,m);
    A = kron(E,B) + kron(D,E);                % Matrix A als Kroneckerprodukt

    
    % Nun können Sie starten mit Ihrer LR-Zerlegung
    % Hilfreich können folgende Befehler sein:
    % tic  und toc für die Zeitmessung 
    % tril und triu  für Dreiecksmatrizen
    
    



%% Graphische Ausgabe: Hier können Sie Ihr Ergebnis zeichen und speichern
% Speichen  des Fenters 1 können Sie zum Beispiel mit Hilfe von
% saveas(1, 'ergebnis_blatt2', 'png')

m = LR(A)
    
%% Funktion LR zur Berechnung der LR-Zerleung der Matrix A
function [m] = LR(A)
    n = size(A, 2);
    m = 
    for(k=1:n-1)
        if(A(k,k) == 0)
            error('Error: all Akk must != 0')
        end
        R(k,k+1:n) = A(k,k+1:n);
        A(k,k+1:n) = R(k,k+1:n);
        L(k+1:n,k) = A(k+1:n,k)/A(k,k);
        A(k+1:n,k) = L(k+1:n,k);
        A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - L(k+1:n,k) * R(k,k+1:n);
    end
    R(n,n) = A(n,n);
    A(n,n) = R(n,n);
end