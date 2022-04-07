distances=pdist2(inc_pix_coordsADJ,[0,0]);
NearestNeighbors=[];
Sortedinc_pix_coordsADJ(:,1:2)=inc_pix_coordsADJ;
Sortedinc_pix_coordsADJ(:,3)=distances;
Sortedinc_pix_coordsADJ=sortrows(Sortedinc_pix_coordsADJ,3);
NearestNeighbors=zeros(length(Sortedinc_pix_coordsADJ));
for q=1:length(Sortedinc_pix_coordsADJ)
dist=pdist2(Sortedinc_pix_coordsADJ(q,:),Sortedinc_pix_coordsADJ);
NearestNeighbors(q,:)=dist;
end

done=0;
qualify=[];
qualify(1)=randi(10);
vector=[];
vector=find(NearestNeighbors(qualify(1):end,qualify(1))>0.06);
qualify(2)=vector(1)+qualify(1)-1;

for x=2:length(Sortedinc_pix_coordsADJ)
    if x>length(qualify)
        break
    end
    vector=find(NearestNeighbors(qualify(x):end,qualify(x))>0.06);
    for y=1:length(vector)
        if NearestNeighbors(vector(y)+qualify(x)-1,qualify(1:x-1))>0.06
            qualify(x+1)=vector(y)+qualify(x)-1;
            vector=[];
            break
        end
    end
end

figure;
blues=Sortedinc_pix_coordsADJ(qualify',1:2);
plot(inc_pix_coordsADJ(:,1),inc_pix_coordsADJ(:,2),'k.');
hold on;plot(blues(:,1),blues(:,2),'ro');
plot(0,0,'go');