function [x,c] = householderreduction(A,b)

c = 0;
[m,n]=size(A);
if m~=n | n~=size(b,1), error('not a square matrix problem'); end;
B=[A b];
N=size(B,2);

bdash = b;

for k=1:n-1
    if A(k,k) == 0
        sig= norm(A((k:n),k));
        c = c + ( n-(k -1) ) + 1;
    else
        sig = sign(A(k,k)) * norm(A((k:n),k));
        c = c + ( n-(k -1) ) + 2;
        
    end
    
    u = A(k:m, k);
    u(1,1) = u(1,1) + sig;
    
    H = sig * ( sig + A(k,k) );
    c = c + 1;
    
    p = eye(n-(k-1)) - u* (u'./H);
    c = c + ( 2 * ((n-(k-1))) ); 
    A (k,k) =-sig; 
    A(k+1:n,k) = 0;
    
    A(k:n,k+1:n) = p * A(k:n,k+1:n);
    c = c + ( (n-(k-1)) *(n-(k-1))*(n-k) );
    
    bdash(k:n,[1:size(b,2)]) = p * bdash(k:n,[1:size(b,2)]);
    c = c +  ( (n-(k-1))^2 * size(b,2) );
    
end
    
B = [A bdash];

% Solution by backsubstitution :
x=zeros(size(b)); % predefinition of x
for k=n:-1:1
  x(k,:)=B(k,n+1:N);
  for j=k+1:n
    x(k,:)=x(k,:)-B(k,j)*x(j,:);
  end
  x(k,:)=x(k,:)/B(k,k);
end