function [x] = MnormVector(M, v)

x = sqrt(v' * M *v);

if imag(x) > 10^-12
    warning('Imaginary part might be to big')
    imag(x)
end
x = real(x);
end