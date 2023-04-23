% First lab. David Rosado Rodríguez

% -------------------------------FIRST EXERCISE---------------------------------



% Read the file

x = load('newcases.txt');

% lengthen the vector to apply the filter
x = [(x(1)*ones(1,2))';x;(x(length(x))*ones(1,2))'];

% Create and perform the filter
b = ones(1,5)/5;
a = zeros(1,5);
a(1) = 1;
y = filter(b,a,x);

% Plot the data before filter and after filtering
figure(1);
subplot(1, 2, 1);
x_axis_vector = 1:length(x);
plot(x_axis_vector,x);
title('Evolution of new cases before filtering');
subplot(1, 2, 2);
x_axis_vector = 1:length(y);
plot(x_axis_vector,y);
title('Evolution of new cases after filtering');

% Compute the absolute value of the Fourier transform of the filter h
function output = H(w)
  output = abs(2*cos(2*w) + 2*cos(w) + 1)/5;
endfunction

% Plot it
figure(2);
x_axis_vector = linspace(-pi, pi, 100);
y = zeros(1, length(x_axis_vector)); 

for i = 1:length(x_axis_vector)
  y(i) = H(x_axis_vector(i));
end

plot(x_axis_vector,y);
title('Abslute value of the Fourier transform of the filter','fontsize',10);

% Analyze why the noise is reduced
figure(3);
freqz(b);

% Propose different filters...

x = load('newcases.txt');

% lengthen the vector to apply the filter
x = [(x(1)*ones(1,3))';x;(x(length(x))*ones(1,3))'];

% Seven days average filter

% Create the filter
b = ones(1,7)/7;
a = zeros(1,7);
a(1) = 1;
y = filter(b,a,x);

% Plot the data before filter and after filtering
figure(4);
subplot(1, 2, 1);
x_axis_vector = 1:length(x);
plot(x_axis_vector,x);
title('Evolution of new cases before filtering');
subplot(1, 2, 2);
x_axis_vector = 1:length(y);
plot(x_axis_vector,y);
title('Evolution of new cases after filtering');

% Compute the absolute value of the Fourier transform of the new filter
function output = H_1(w)
  output = abs(2*cos(3*w)+2*cos(2*w) + 2*cos(w) + 1)/7;
endfunction

% Plot it
figure(5);
x_axis_vector = linspace(-pi, pi, 100);
y = zeros(1, length(x_axis_vector)); 

for i = 1:length(x_axis_vector)
  y(i) = H_1(x_axis_vector(i));
end

plot(x_axis_vector,y);
title('Abslute value of the Fourier transform of the new filter','fontsize',10);




% -------------------------------SECOND EXERCISE--------------------------------

% Read the file
x = load('heights.txt');

% Compute the total height climbed and descended
climbed = 0;
descended = 0;

for n = 1:length(x)-1
    dif = x(n+1) - x(n);
    if dif > 0
        climbed += dif;
    endif
    if dif < 0
        descended -= dif;
    endif
end

% Print the climbed and descended distance

fprintf('The total climbed distance before filtering is: %.2f\n',climbed);
fprintf('The total descended distance before filtering is: %.2f\n', descended);

% Apply the previous filter 

% lengthen the vector to apply the filter
x = [(x(1)*ones(1,2))';x;(x(length(x))*ones(1,2))'];

% Create and perform the filter
b = ones(1,5)/5;
a = zeros(1,5);
a(1) = 1;
y = filter(b,a,x);

% Plot the data before filter and after filtering
figure(6);
subplot(1, 2, 1);
x_axis_vector = 1:length(x);
plot(x_axis_vector,x);
title('Evolution of heights before filtering');
subplot(1, 2, 2);
x_axis_vector = 1:length(y(5:end));
plot(x_axis_vector,y(5:end));
title('Evolution of heights after filtering');

% Compute the total height climbed and descended
climbed = 0;
descended = 0;

for n = 5:length(y)-1
    dif = y(n+1) - y(n);
    if dif > 0
        climbed += dif;
    endif
    if dif < 0
        descended -= dif;
    endif
end

% Print the climbed and descended distance

fprintf('The total climbed distance after filtering is: %.2f\n',climbed);
fprintf('The total descended distance after filtering is: %.2f\n', descended);


