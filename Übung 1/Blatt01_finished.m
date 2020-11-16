% Blatt 01

%Beispielmatrizen aus Beispiel 2.3, A=LR
L= [1, 0, 0;
    2/3, 1, 0;
    1/3, 2, 1];
R= [3, 1, 6;
    0, 1/3, -1;
    0, 0, 1];
b= [2; 7; 4];

%Loesen Sie hier die Gleichung Ax=b mit Hilfe der beiden Funktionen
%vorwaerts und rueckwaerts


x = rueckwaerts(R,vorwaerts(L,b))


%Implementieren Sie hier die Funktion vorwaerts(L,b), die eine Gleichung
%der Form Lx=b loest. L normierte untere linke Dreiecksmatrix
function [y]=vorwaerts(L,b)
    n=size(L,2);
    y=zeros(n,1);
    y(1) = b(1)/L(1,1);
    for j=2:n
        y(j) = (b(j) - (L(j,1:j-1)*y(1:j-1))) / L(j,j);
    end
end

%Implementieren Sie hier die Funktion rueckwaerts(L,b), die eine Gleichung
%der Form Rx=b loest. R rechte obere Dreiecksmatrix
function [x]=rueckwaerts(R,y)
    n=size(R,2);
    x=zeros(n,1);
    x(n) = y(n)/R(n,n);
    for j=1:n-1
        x(n-j) = (y(n-j) - (R(n-j,n-j+1:n)*x(n-j+1:n))) / R(n-j,n-j);
    end
end