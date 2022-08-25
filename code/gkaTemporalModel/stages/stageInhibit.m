function arg = stageInhibit(f,fc,k)
arg = (1i.*f+(1-k)*fc).^2 ./ (1i.*f+fc).^2 ;
end