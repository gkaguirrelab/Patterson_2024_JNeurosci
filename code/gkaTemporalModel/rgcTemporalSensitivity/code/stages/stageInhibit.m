function arg = stageInhibit(fc,k)
syms f
arg = (1i.*f+(1-k)*fc).^2 ./ (1i.*f+fc).^2 ;
end