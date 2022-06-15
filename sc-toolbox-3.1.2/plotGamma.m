gamma =  linspace(1.01,2.1,10)
l = length(gamma);
R = 1:l;

for j = 1:1
    gammaN = 1.001

    p = polygon([1 gammaN gammaN+sqrt(3)*gammaN*1i 1+sqrt(3)*gammaN*1i]);
    plot(p)
    f = extermap(p);
    g=inv(f);

    %R(j) = abs(g(0));
    R(j) = abs(evalinv(f,0,10^-30));
end

plot(gamma, R)
