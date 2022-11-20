function [arg,f] = stageFirstOrderLP(fc,n)
syms f
arg = 1./(1i.*f+fc).^n ;
end