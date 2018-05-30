function [i] = LPP(x)
n=x;
x(1)=n(2);
x(2)=n(1);
i=x;
end

