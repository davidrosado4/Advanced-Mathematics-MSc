% IMPLEMENTATION OF THE INVERSE FAST FOURIER TRANSFORM ALGORITHM

function y = recursive_ifft(a)
  n = length(a);
  if n == 1
    y = a;
  else
    w_n = exp((2*pi*1i)/n);
    w = 1;
    a_even = a(1:2:end);
    a_odd = a(2:2:end);
    y_even = recursive_ifft(a_even);
    y_odd = recursive_ifft(a_odd);
    y = zeros(1, n);
    for k = 0:(n/2)-1
      y(k+1) = y_even(k+1) + w*y_odd(k+1);
      y(k+(n/2)+1) = y_even(k+1) - w*y_odd(k+1);
      w = w*w_n;
    endfor
  endif
endfunction
