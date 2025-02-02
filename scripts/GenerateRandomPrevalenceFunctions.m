
function [pl, pu, p] = generateRandomPrevalenceFunctions(T)

prd = rand()+1;
a = 0.5*rand();
b = a + (1-2*a)*rand();

w = 2*prd*pi/T;
phase = T*rand();
phi = w*phase;

k = floor((w*T+phi)/2/pi);
t = (2*pi*k-phi)/w;

cmin = min([1/T,1/t*(1/(a+b)-1)]);

c = cmin*rand();

dmin = min([1-c*T, 1/(a+b)-1-c*t]);

d = dmin*rand();

x = 1:T;
p = a*cos(w*x+phi)+b;
r = c*x + d;
pl = p(x).*(1-r(x));
pu = p(x).*(1+r(x));

pmax = max(pu);

if pmax>1
    p = p/pmax;
    pl = pl/pmax;
    pu = pu/pmax;
end
% close all
% figure
% hold on
% plot(x,p);
% plot(x,pl);
% plot(x,pu);
% ylim([0,1])

end

% p = @(x) a*cos(w*x+phi)+b;
% r = @(x) c*x + d;
% pl = @(x) p(x).*(1-r(x));
% pu = @(x) p(x).*(1+r(x));

% close all
% figure
% hold on
% ezplot(p,[0,T]);
% ezplot(pl,[0,T]);
% ezplot(pu,[0,T]);
% ylim([0,1])
