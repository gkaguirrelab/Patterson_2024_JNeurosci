%Round the cone locations to some number of decimals places after zero
RoundedLocations=ceil(SortedTotalMosaic.*1000)/1000;
RoundedLCoords=ceil(LCoords.*1000)/1000;

Fx=zeros(10001);
%lin=linspace(-.5,.5,1000);
lin=(-.5:.0001:.5);
for f=1:length(Fx(:,1))
Fx(f,:)=lin;
end
%dx=ceil(dx.*100)/100;

Fy=zeros(10001);
%lin=linspace(-.5,.5,1000);
lin=(-.5:.0001:.5);
for f=1:length(Fy(1,:))
Fy(:,f)=lin';
end
Fy=flipud(Fy(:,1:end));
%dy=ceil(dy.*100)/100;

% dd=zeros(10001,10001,3);
% dd(:,:,1)=dx;
% dd(:,:,2)=dy;

FXY=zeros(10001,10001);

for p=1:length(RoundedLCoords)
    if isnan(RoundedLCoords(p,1))==0
        testx(p,:)=find(round(Fx(1,:)-RoundedLCoords(p,1),4)==0);
        testy(p,:)=find(round(Fy(:,1)-RoundedLCoords(p,2),4)==0);
        FXY(testy(p,:),testx(p,:))=1;
    else
        testx(p,:)=0;
        testy(p,:)=0;
    end
end


figure;imagesc(FXY);

cross=grating.*FXY;
crossgabor=gabor.*FXY;
crossgauss=gauss.*FXY;
normcrossgauss=crossgauss/(sum(sum(crossgauss)));

for c=1:length(testx)
    values(c,:)=crossgauss(testy(c),testx(c)); %LDoG, etc. from GenerateCell
end

crossgauss=zeros(10001,10001);
for d=1:length(testx)
    if testx(d,:)>0
        crossgauss(testy(d),testx(d))=LDoG(d);
    end
end

figure;
plot3(RoundedLocations(:,1),RoundedLocations(:,2),values,'r.')
axis square
figure;
scatter(RoundedLocations(1:ConesToSurround,1),RoundedLocations(1:ConesToSurround,2),30,values(1:ConesToSurround),'fill')
axis square