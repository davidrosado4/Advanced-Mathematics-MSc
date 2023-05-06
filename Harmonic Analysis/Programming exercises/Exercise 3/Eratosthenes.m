% SIEVE OF ERATOSTHENES ALGORITHM 

function primes = Eratosthenes(n)

    % Create an array of boolean values to indicate whether each number is prime or not
    is_prime = true(1, n);

    % Mark 1 as not prime (it's the only even prime)
    is_prime(1) = false;

    % Mark all multiples of 2 as not prime (except 2)
    for i = 4:2:n
        is_prime(i) = false;
    end

    % Loop over odd numbers up to sqrt(n), marking their multiples as not prime
    for i = 3:sqrt(n)
        if is_prime(i)
            for j = i^2:i*2:n
                is_prime(j) = false;
            end
        end
    end
    primes = is_prime;
end
