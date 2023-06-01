function arg = stageSecondOrderLP(fc,Q)
syms f
% Sallen-Key filter
    arg = 1 ./ ( -(f./fc).^2 + ((1i.*f)./(Q*fc)) +1 );
end