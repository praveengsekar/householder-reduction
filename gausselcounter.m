%
%  file gausselcounter.m
%
%  [x,mulel,mulsub] = gausselcounter(A,b) returns the solution of the set of
%  linear algebraic equations A*x=b using Gauss elimination
%  without pivoting and the count of mul/division. 
%  b and x can have several columns
%
function [x,mulel] = gausselcounter(A,b)

mulel = 0;

[m,n]=size(A);
if m~=n || n~=size(b,1), error('not a square matrix problem'); end;

B=[A b];
N=size(B,2);

% bring the matrix into triangular form (Gauss elimination):

for k=1:n-1    % loop over columns where the zeros will appear
  fac=1/B(k,k);
  
  for i=k+1:n   % loop over rows where subtractions take place
    fac1=fac*B(i,k); % factor B(i,k)/B(k,k);
    mulel = mulel+1;
    B(i,k)=0; % new zero by construction
    B(i,k+1:N)=B(i,k+1:N)-B(k,k+1:N)*fac1; % subtraction
    mulel = mulel+ N - k;
  end
end

% Solution by backsubstitution :
x=zeros(size(b)); % predefinition of x
for k=n:-1:1
  x(k,:)=B(k,n+1:N);
  for j=k+1:n
    x(k,:)=x(k,:)-B(k,j)*x(j,:);
  end
  x(k,:)=x(k,:)/B(k,k);
end