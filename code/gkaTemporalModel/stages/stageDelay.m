function arg = stageDelay(d)
syms f
arg = exp(-1i.*(d*2*pi)*f);
end