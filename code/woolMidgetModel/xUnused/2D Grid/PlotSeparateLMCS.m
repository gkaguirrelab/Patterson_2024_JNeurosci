figure;

subplot(2,2,1);
imagesc(X,X,AllLc);
title('L Center Response Profile');
xlabel('x (deg)');
ylabel('y (deg)'); 
axis square;
axis([-SurroundRadiusDeg*2 SurroundRadiusDeg*2 -SurroundRadiusDeg*2 SurroundRadiusDeg*2]);

subplot(2,2,2);
imagesc(X,X,AllLs);
title('L Surround Response Profile');
xlabel('x (deg)');
ylabel('y (deg)'); 
axis square;
axis([-SurroundRadiusDeg*2 SurroundRadiusDeg*2 -SurroundRadiusDeg*2 SurroundRadiusDeg*2]);

subplot(2,2,3);
imagesc(Freq,Freq,Lcpower);
title('L Center Frequency Domain');
xlabel('Frequency (cpd)');
ylabel('Frequency (cpd)'); 
axis square;

subplot(2,2,4);
imagesc(Freq,Freq,Lspower);
title('L Surround Frequency Domain');
xlabel('Frequency (cpd)');
ylabel('Frequency (cpd)'); 
axis square;

figure;

subplot(2,2,1);
imagesc(X,X,AllMc);
title('M Center Response Profile');
xlabel('x (deg)');
ylabel('y (deg)'); 
axis square;
axis([-SurroundRadiusDeg*2 SurroundRadiusDeg*2 -SurroundRadiusDeg*2 SurroundRadiusDeg*2]);

subplot(2,2,2);
imagesc(X,X,AllMs);
title('M Surround Response Profile');
xlabel('x (deg)');
ylabel('y (deg)'); 
axis square;
axis([-SurroundRadiusDeg*2 SurroundRadiusDeg*2 -SurroundRadiusDeg*2 SurroundRadiusDeg*2]);

subplot(2,2,3);
imagesc(Freq,Freq,Mcpower);
title('M Center Frequency Domain');
xlabel('Frequency (cpd)');
ylabel('Frequency (cpd)'); 
axis square;

subplot(2,2,4);
imagesc(Freq,Freq,Mspower);
title('M Surround Frequency Domain');
xlabel('Frequency (cpd)');
ylabel('Frequency (cpd)'); 
axis square;