clear all

syms a b c
data.a = a;
data.b = b;
data.c = c;

l = a+b+c;
lf = matlabFunction(l,'Vars',data)