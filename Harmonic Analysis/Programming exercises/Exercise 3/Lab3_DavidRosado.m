% David Rosado Rodríguez


% -------------------------------EXERCISE 1-------------------------------------
printf("\n\nComputing exercise 1...\n\n");
% Compare our fft with the already implemented fft

% Create  a random variable
a = 1:32;

y_myfft = recursive_fft(a);
y_fft = fft(a);
printf("The error between my fft and the already implemented is : %e\n", norm(y_myfft - y_fft));

% Compare our ifft with the already implemented ifft

a_myfft = recursive_ifft(y_fft)/length(y_fft);
a_fft = ifft(y_fft);
printf("The error between my ifft and the already implemented is : %e\n", norm(a_myfft - a_fft));

% ------------------------------EXERCISE 2--------------------------------------

% Start the timer
tic

printf("\n\nComputing exercise 2...\n\n");
% Find the first 100.000 primes using Eratosthenes function

n = 1299709;      % This number is the prime number 100.000
is_prime = Eratosthenes(n);

% Find the indices of all the primes
primes = find(is_prime);

% COMPUTATION OF P^2

% Coefficients of the polynomial p
a = is_prime;

% Pad the polynomial with zeros to the nearest power of 2
deg = 2*n + 2;
m = 2^nextpow2(deg);
a = [0 ,a, zeros(1, m-length(is_prime)-1)];

% Compute the FFT of the padded polynomial
A = fft(a);

% Element-wise multiplication of the FFT
B = A .* A;

% Compute the inverse FFT of the result
b = ifft(B);

% Extract the coefficients of the resulting polynomial
b = real(b(1:deg));

% Substract 1 and divide by 2 the coefficients
b = (b-1)/2;

% Count of primes
twotimes_primes = 2*primes;
count = 0;
for i = 1:length(twotimes_primes)
  triplets = b(twotimes_primes(i) + 1);
  count = count + triplets;
end

printf("The number of triplets found in the first 100.000 primes is: %d\n", count);

% Stop the timer and display the elapsed time

elapsed_time = toc;
disp(['Time required to perform the second exercise: ' num2str(elapsed_time) ' seconds']);
