function u1 = TVD(u,lam)
f = inline('u.^2/2');
N = length(u); 
du1 = zeros(N,1); du2 = zeros(N,1); Du1 = zeros(N,1); Du2 = zeros(N,1);
r = zeros(N,1); u1 = zeros(N,1);

um1 = circshift(u,1); um2 = circshift(um1,1); up1 = circshift(u,-1);
ubar = 1/2*(u+um1); Lam = lam*ubar;
for i = 1:N
    if u(i) == um1(i)
        r(i) = 1;
    elseif ubar(i)>0
        du1(i) = um1(i)-um2(i); du2(i) = u(i)-um1(i);
        r(i) = du1(i)/du2(i);
    elseif ubar<0
        Du1(i) = up1(i)-u(i); Du2(i) = u(i)-um1(i);
        r(i) = Du1(i)/Du2(i);
    end
end


phim = vanLeer(r); phip = circshift(phim,-1);

Id1 = ubar>0; Id2 = ubar<0;
if ~isempty(Id1)
    u1(Id1) = u(Id1)-lam*(f(u(Id1))-f(um1(Id1)))-...
        1/2*Lam(Id1).*(1-Lam(Id1)).*(phip(Id1).*(up1(Id1)-u(Id1))-phim(Id1).*(u(Id1)-um1(Id1)));
end
if ~isempty(Id2)
    u1(Id2) = u(Id2)-lam*(f(up1(Id2))-f(u(Id2)))-...
        1/2*Lam(Id2).*(1+Lam(Id2)).*(phip(Id2).*(up1(Id2)-u(Id2))-phim(Id2).*(u(Id2)-um1(Id2)));
end

