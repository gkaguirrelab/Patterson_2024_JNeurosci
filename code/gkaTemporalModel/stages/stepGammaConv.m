function irfOut = stepGammaConv(irfIn,deltaT,t1,w)

gammaTime = 0:deltaT:1;
gammaKernel = gammaTime.*exp(-gammaTime/t1)-w.*gammaTime.*exp(-gammaTime/(1.5*t1));
gammaKernel = gammaKernel / sum(gammaKernel);
temp = conv(irfIn,gammaKernel);
irfOut = temp(1:end-length(gammaKernel)+1);

end