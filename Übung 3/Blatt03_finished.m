% Blatt 03

% Hilbertmatrizen koennen mit der Funktion hilb erstellt werden
n=8;      H=hilb(n);


%Hinweise: 
% 1) y=noisyVector(b,p)   kann zum Generieren gestoerter rechter Seiten
% verwendet werden
% 2) die Funktionen vorwaerts und rueckwaerts vom Blatt 1 koennen ebenfalls
% verwendet werden
% 3) semilogy  koennen Sie zum plotten verwenden
% 4) Hilfreich koennen auch cond und norm sein
% Beachten Sie bei allen Matlabfunktionen die optionalen Parameter


Nmax=10; 
kond=zeros(Nmax,1);
fehlerH=zeros(Nmax,1);
fehlerX=zeros(Nmax,1); fehlerB=fehlerH;



% Aufgabe 4.2

L = cholesky(H)

%kond_H ist die Konditionszahl von H
kond_H = s_max(H) * s_max(inv(H))

%F_rel = relativer Fehler der Zerlegung
F_rel = s_max(H - L * L.') / s_max(H)


% Aufgabe 4.3
for i=2:Nmax
    H = hilb(i);
    L = cholesky(H);
    b = s_vec(H) + noisyVector(s_vec(H),1);

    x_s = rueckwaerts(L.',vorwaerts(L, b));

    fehlerH(i) = s_max(H - L * L.') / s_max(H);
    fehlerX(i) = 1 / max(x_s);
    
end

vec_x = [1 2 3 4 5 6 7 8 9 10];
semilogy(vec_x,fehlerH,vec_x,fehlerX)

% Funktion, welche das Maximum der Zeilensummen von L als return Wert hat
% (also ||.||_1)
function [Ze]=z_max(L)
    n = size(L,2);
    Z = zeros(n,1);
    for i = 1:n
        Z(i) = sum(L(i,1:n));
    end
    Ze = max(Z);
end

% Funktion, welche das Maximum der Spaltensummen von L als return Wert hat
% (also ||.||_infty)
function [Ze]=s_max(L)
    n = size(L,2);
    Z = zeros(n,1);
    for i = 1:n
        Z(i) = sum(L(i:n,i));
    end
    Ze  = max(Z);
end

% Funktion, welche die Spaltensummen von L als return Vektor hat
% (also ||.||_1)
function [Z]=s_vec(L)
    n = size(L,2);
    Z = zeros(n,1);
    for i = 1:n
        Z(i) = sum(L(i:n,i));
    end
end

%Implementieren Sie hier den Cholesky-Algorithmus
function [L]=cholesky(A)
    n = size(A,2);
    L = zeros(n,2);
    for i = 1:n
        L(i,i) = sqrt(A(i,i) - sum(L(i,1:i-1).^2));
        if (L(i,i) <= 0 | isreal(L) == 0)
            error('Error: A is not positive definite')
        end
        for j = i+1:n
            L(j,i) = (A(j,i) - sum(L(j,1:i-1) .* L(i,1:i-1))) / L(i,i);
        end
    end
end





% Funktion zum Addidieren von p% Rauschen 
function y=noisyVector(b,p)
% y=b+noise, wobei noise normierter Fehlervektor (normalverteilt) ist
% der relative Fehler ist ||y-b||/||b|| = p/100

    delta=randn(size(b));                            % Zufallsvektor delta
    noise=delta/norm(delta,'inf')*norm(b,'inf');     % normieren des Fehlervektors   
    y = b+ noise* p/100;                             % additiver Fehler (p%)
end



% Implementation der Vorwaerts- und Rueckwaertselimination 
function [x]=rueckwaerts(R,b)
n=size(R,2);
x=zeros(n,1);
for i=n:-1:1
    x(i)=(b(i)-R(i,i+1:n)*x(i+1:n))/R(i,i);
end
end

% Im Vergleich zu Blatt 1 ist die vorwaerts Funktion leicht modifiziert
% worden fuer allgemeine untere Dreiecksmatrizen L 
function [x]=vorwaerts(L,b)
n=size(L,2);
x=zeros(n,1);
for i=1:n
    x(i)=(b(i)- L(i,1:i-1)*x(1:i-1))/L(i,i);
end
end

