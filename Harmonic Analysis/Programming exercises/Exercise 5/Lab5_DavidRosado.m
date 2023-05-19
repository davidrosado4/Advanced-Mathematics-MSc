% David Rosado Rodrígez


%--------------------------------EXERCISE 1-------------------------------------

% Generate 35 random points between 0 and 400
points = rand(35, 1) * 400;

% Define the polynomial
function output = p(x)
    output = 2*sin(x) + sin(2*x) + 21*sin(300*x);
endfunction

% Create the matrix A to solve Ax = y
A = zeros(35,400);
for i=1:35
  for j=2:400
    A(i,1) = 1;
    A(i,j) = sin((j-1)*points(i));
  endfor
endfor

% Create the vector y to solve Ax = y
y = zeros(1,35);
for i=1:35
  y(i) = p(points(i));
endfor

% Solve the minimization problem
c = ones(1,400);
x = glpk(c,A,y);

% Print the values different from 0
disp("The coefficients different from 0 are");
for i=1:length(x)
  if abs(x(i))>1e-10
    disp(['x(', num2str(i-1), ') = ', num2str(x(i))]);
  endif
endfor
