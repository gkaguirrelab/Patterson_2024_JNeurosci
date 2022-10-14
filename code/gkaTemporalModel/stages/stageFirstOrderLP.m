function [arg,f] = stageFirstOrderLP(f,fc,n)
arg = 1./(1i.*f+fc).^n ;
end