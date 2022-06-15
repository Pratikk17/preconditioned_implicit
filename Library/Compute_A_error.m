function [err]=Compute_A_error(x,x_exact,A)

err=zeros(size(x,2),1);
for i=1:size(x,2)
    z(i) = sqrt( (x(:,i)-x_exact)'*A*(x(:,i)-x_exact) );
    err(i) = real(z(i));
end
if max(abs(imag(z))) > 10^-14
        warning('Imaginray part to big')
        max(abs(imag(z)))
end
end
