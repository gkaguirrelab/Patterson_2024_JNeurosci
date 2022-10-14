
n=22;

[phat, ~, mu] = polyfit(myTime', irfLumLGN', n);
flipPhat = flip(phat); % Flip order to [p0, ..., pn]
p = zeros(size(flipPhat));
for ii = 0:n
    for k = ii:n
        p(ii+1) = p(ii+1) + nchoosek(k, k-ii) * flipPhat(k+1)/mu(2)^k * (-mu(1))^(k-ii);
    end
end
p = flip(p); % Back to original order [pn, ..., p0]

yFit = polyval(phat,myTime',[],mu);
fitIRFSym = poly2sym(p,x);
ySymFit = eval(subs(fitIRFSym,x,myTime));

figure
plot(myTime,irfLumLGN,'.k')
hold on
plot(myTime,yFit,'-r')
plot(myTime,ySymFit,'-g')

fitTTFSym = subs(fourier(fitIRFSym,x,w),w,f*2*pi);
plotRF(fitTTFSym)
