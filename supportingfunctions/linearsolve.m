function x=linearsolve(L,U,b)
y=ltri(L,b);
x=utri(U,y);

function y = ltri(L,b)
n=size(b,1);
y=zeros(n,1);
for j = 1:n-1
    y(j) = b(j)/L(j,j);
    b(j+1:n) = b(j+1:n) - y(j) * L(j+1:n,j);
end
y(n) = b(n)/L(n,n);
end

function x = utri(U,y)
n=size(y,1);
x=zeros(n,1);
for j = n:-1:2
    x(j) = y(j)/U(j,j);
    y(1:j-1) = y(1:j-1) - x(j) * U(1:j-1,j);
end
x(1) = y(1)/U(1,1);
end
end