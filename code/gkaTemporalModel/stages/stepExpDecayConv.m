function irfOut = stepExpDecayConv(irfIn,deltaT,t2)

expTime = 0:deltaT:5;
expKernel = exp(-expTime/t2);
expStart = find(irfIn>max(irfIn)/1e3,1);
expKernel = circshift(expKernel,expStart);
expKernel(1:expStart)=0;
expKernel = expKernel./sum(expKernel);
temp = conv(irfIn,expKernel);
irfOut = temp(1:end-length(expKernel)+1);

end