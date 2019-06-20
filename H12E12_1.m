n = 3;
A = rand(n);
b = rand(n,1);
[x1,c] = householderreduction(A,b);
fprintf ('The solution for a random matrix\n');
fprintf ('%d\n', x1);

A1 = [0 4 5; 4 5 1; 1 6 7];
b1 = [1 3 2]';
[x2,c] = householderreduction(A1,b1);
[y2,c] = gausselcounter(A1,b1);
fprintf ('\nhouseholders\n');
fprintf ('%d\n', x2);
fprintf ('\ngauss-elimination\n');
fprintf ('%d\n', y2);

for n = 2:100
    A = rand(n,n); 
    b = rand(n,1); 
    
    [x,c] = gausselcounter(A,b);
    [y,c] = householderreduction(A,b);
    z = A\b;
    
    dx(n-1) = max ( abs (x-z) );
    dy(n-1) = max ( abs (y-z) );
end

p1 = semilogy([2:100],dx,'b*');
hold on
p2 = semilogy([2:100],dy,'ro');
xlabel('n') 
ylabel('Error')
legend([p1 p2],'Gaussian Elimination','Householder Reduction')

n = 1000;
A = rand(n);
b = rand(n,1);
[xh1,ch] = householderreduction(A,b);
[xg1,cg] = gausselcounter(A,b);
fprintf('\nThe ratio of counters for large n = %d' , ch/cg);