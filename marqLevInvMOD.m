%% Function Forward Solution for Marq-Levn Algorithm

function[g] = marqLevInvMOD(r,t,s)
q = 13;
f = 10;
m = 4.438;
x = 0;
e = exp(0.5*log(10)/m);
h = 2*q-2;
u = s*exp(-f*log(10)/m-x);
l = length(r);
n = 1;
for i = 1:n+h
w = l;
v = r(l);
while w>1
w = w-1;
aa = tanh(t(w)/u);
v = (v+r(w)*aa)/(1+v*aa/r(w));
end
a(i) = v;
u = u*e;
end
i = 1;
g = 105*a(i)-262*a(i+2)+416*a(i+4)-746*a(i+6)+1605*a(i+8);
g = g-4390*a(i+10)+13396*a(i+12)-27841*a(i+14);
g = g+16448*a(i+16)+8183*a(i+18)+2525*a(i+20);
g = (g+336*a(i+22)+225*a(i+24))/10000;
end
