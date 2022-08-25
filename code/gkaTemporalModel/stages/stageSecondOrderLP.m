function arg = stageSecondOrderLP(f,fc,Q)
% Sallen-Key filter
    arg = 1 ./ ( -(f./fc).^2 + ((1i.*f)./(Q*fc)) +1 );
end