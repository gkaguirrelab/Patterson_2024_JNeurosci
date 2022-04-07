%Generate a hex grid of cell positions in our universe, starting from the
%origin and working outward.
coords=[0 0 0];
radius=100;
for i=1:radius
    for j=-i:i
        for k=-i:i
            for l=-i:i
                if abs(j)+abs(k)+abs(l)==i*2 && j+k+l==0
                    %disp(strcat(num2str(j),num2str(k),num2str(l)));
                    coords(end+1,:)=[j k l];
                end
            end
        end
    end
end

for m=1:length(coords)
x = (sqrt(3)*(coords(m,1) + (coords(m,2)/2)));
y = ((3/2)*coords(m,2));
pixel_coords(m,:)=[x y];
end

%pixel_coords(:,1)=jitter(pixel_coords(:,1));
%pixel_coords(:,2)=jitter(pixel_coords(:,2));

plot(pixel_coords(:,1),pixel_coords(:,2),'bo')
axis equal

Em=10;
ConesPerMM=ceil((52170*exp(-1.057.*Em))+(15550*exp(-0.1605.*Em))); %From 'cones per midget gc data' Excel file
DFRadius=0.002738*Em.^1.327;

r=1;
ConeCount=0;
while ConeCount<ConesPerMM
    r=r+.1;
    counter=0;
    inc_pix_coords=[];
    for p=1:length(pixel_coords)
        if pixel_coords(p,1)>-r && pixel_coords(p,1)<r
            if pixel_coords(p,2)>-r && pixel_coords(p,2)<r
                inc_pix_coords(end+1,1)=pixel_coords(p,1);
                inc_pix_coords(end,2)=pixel_coords(p,2);
                counter=counter+1;
            end
        end
    end
    ConeCount=counter;
end
rectangle('Position',[-r,-r,2*r,2*r],'EdgeColor','r');

%Adjust the aperture window to a 1-mm2 portion 
RFactor=1/(2*r);
pixel_coordsADJ=pixel_coords*RFactor;
inc_pix_coordsADJ=inc_pix_coords*RFactor;
figure;hold on;
%plot(pixel_coordsADJd(:,1),pixel_coordsADJd(:,2),'bo');
plot(inc_pix_coordsADJ(:,1),inc_pix_coordsADJ(:,2),'ro');
%rectangle('Position',[-r*RFactor,-r*RFactor,2*r*RFactor,2*r*RFactor],'EdgeColor','r');
axis([-.5 .5 -.5 .5]);
axis square
