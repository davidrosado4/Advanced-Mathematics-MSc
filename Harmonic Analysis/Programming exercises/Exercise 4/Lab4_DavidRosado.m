% David Rosado Rodríguez


%---------------------------------EXERCISE 1------------------------------------

% Read the image and visualize it
N = 256;
F = phantom(N);
figure(1);
imshow(F);

% Using pre-built functions...
[G,xp] = radon(F,0:179);
FR = iradon(G,0:179);
figure(2);
imshow(FR);

% Building our functions...

% Start the timer
tic;

% Settings
omega = 0.5;
length_grid = length(xp);
h = zeros(length_grid,180);
thetas = zeros(1,180);
q = (length_grid-1)/2;
MAT_recover = zeros(N,N);

% Function that compute Ram_Lak filter
function output = Ram_Lak(s,omega)
  if s == 0
    output = 0.5*omega*omega;
  else
    value = 2*pi*omega*s;
    output = sin(value)/value - 0.5*((sin(value/2)/(value/2)) *(sin(value/2)/(value/2)));
    output = output*omega*omega;
  endif
endfunction

% Step 1
disp("Computing the first step of the algorithm...");
for j = 1:180
  for k = 1:length_grid
    for l = 1:length_grid
      h(k,j) = h(k,j) + Ram_Lak(xp(k)-xp(l),omega)*G(l,j);
    endfor
  endfor    
endfor

% Step2
disp("Computing the second step of the algorithm...");
% Create the theta vector
for j = 1:180
  thetas(j) = exp(1i*(j*pi/180));
endfor
% Recovering the image...
for row = 1:N
  for col = 1:N
    x = -((N+1)/2 - col);
    y = - (row-(N+1)/2);
    for j = 1:180
      eta = x*real(thetas(j)) + y*imag(thetas(j));
      k = floor(eta);
      eta = eta - k;
      k = k + q + 1; % k goes from -q to q --> move it to 1 to length_grid
      MAT_recover(row,col) = MAT_recover(row,col) + (1-eta)*h(k,j) + eta*h(k+1,j);
    endfor
    MAT_recover(row,col) = MAT_recover(row,col)*(2*pi/180);
  endfor
endfor

% Show the generated image with our functions
figure(3);
imshow(MAT_recover);

% Stop the timer and display the elapsed time
elapsed_time = toc;
disp(['Time to compute our iradon function: ' num2str(elapsed_time) ' seconds']);
  
  
  
  