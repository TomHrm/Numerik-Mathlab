%Blatt 04

kmax=6;                     % bestimmt die maximale Groesse von A 
Nmax=20000;                 % maximale Anzahl von Iterationen fuer Ihre Schleife
eps=1e-8;


    k=2;
    
    
    for k=1:6
        m=2^k;  n=m^2;
        % Matrix A als Blockmatrix (vom Blatt 2)
        e = ones(m,1);
        B = spdiags([-e 4*e -e], -1:1, m, m);   % Bandmatrix B
        E=speye(m);                             % Identitaet                     
        D=spdiags([-e -e], [-1,1],m,m);
        A=kron(E,B) + kron(D,E);                % Matrix A als Kroneckerprodukt
   
        xb=ones(n,1);
        b=A*xb;                                 % Rechte Seite des LGS
    
        x = zeros(n,1);
    
        x = Jacobi(A, x, b, eps)
    end
    

    % Hinweise
    % statt einer for-schleife  kann auch eine while-schleife verwendet werden
    % Setzten Sie sicherheitshalber eine maximale Anzahl Nmax maximaler
    % Iterationen fest
    % Zum Schaetzer fuer q -  setzen Sie ein sinnvolles q_0 fest
    
    
   
%% Implementation des Jacobi-Verfahrens

% Anmerkung zur Implementation des Jacobi-Verfahrens:
% Laut Algorithmus 6 des Skripts, wird eine Matrix B von der Größe n x n,
% also von der gleichen Größe wie Matrix A übergeben.
% in der gegebenen Matlab-Datei wurde ein x angegeben, welcher bei
% Nachfrage beim Gruppentutor der Startvektor x sein soll. Warum???

% Da nicht Matrix B nicht bekannt war bzw. ich nicht wusste, ob ich die
% irgendwie errechnen kann o.Ä. wurde für Kompilationszwecke stattdessen
% Matrix A an den (angegebenen) Stellen benutzt.
% Diese Umstände haben das Bearbeiten diese Woche leider nicht so spaßig
% gemacht :/
% Da dadurch keine vernünftigen Ergebnisse zustande kamen, werden diese
% auch nicht visualisiert.
function [x] = Jacobi(A,x,b, eps)
    abbr_krit = false;
    r_0 = b - A * x;
    %in der folgenden Zuweisung müsste statt der Matrix A, Matrix B benutzt
    %werden. Dies war leider nicht möglich (s.o.)
    delta_x = r_0' / A;
    x_old = x;
    x = x + delta_x;
    while (abbr_krit == false)
        r = b - A * x;
        %in der folgenden Zuweisung müsste statt der Matrix A, Matrix B benutzt
        %werden. Dies war leider nicht möglich (s.o.)
        delta_x = r' / A;
        x_old_2 = x_old;
        x_old = x;
        x = x + delta_x;
        if ((max(r) <= eps * r_0) | (max(delta_x) <= eps * max(x)) | (max(x - x_old) <= eps * (1 - ((x - x_old) / (x_old - x_old_2))) * max(x_old)))
            abbr_krit = true;
        end
    end
end
