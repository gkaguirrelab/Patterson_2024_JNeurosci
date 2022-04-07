for c=1:length(Cells)
    ImportCell=Data.(matlab.lang.makeValidName(Cells{c}));

    Em=ImportCell.Eccentricity;
    ChromTag=ImportCell.ChromTag;
    TotalLData=ImportCell.ResponseFunctions.LIsoResponse;
    TotalMData=ImportCell.ResponseFunctions.MIsoResponse;
    LMSumData=ImportCell.ResponseFunctions.LMSumResponse;
    LMDiffData=ImportCell.ResponseFunctions.LMDiffResponse;

    if strcmp(ChromTag,'L-dominated')==1
        WeakCone=TotalMData;         
    elseif strcmp(ChromTag,'M-dominated')==1
        WeakCone=TotalLData;
    elseif strcmp(ChromTag,'Achromatic')==1
        WeakCone=LMSumData;
    end
    
    if strcmp(ChromTag,'Achromatic')==1
        OSI(c,:)=NaN;
    elseif strcmp(ChromTag,'Achromatic')==0
        for t=1:length(thetas)
            OrientationLabel=strcat('Theta',num2str(thetas(t)),'deg');
            FirstAmp=WeakCone.(matlab.lang.makeValidName(OrientationLabel)).Amplitude(1);
            MaxAmp=max(WeakCone.(matlab.lang.makeValidName(OrientationLabel)).Amplitude);
            MaxIndex=find(WeakCone.(matlab.lang.makeValidName(OrientationLabel)).Amplitude==MaxAmp);
            NotchAmp=min(WeakCone.(matlab.lang.makeValidName(OrientationLabel)).Amplitude(1:MaxIndex));
            NotchIndex=find(WeakCone.(matlab.lang.makeValidName(OrientationLabel)).Amplitude(1:MaxIndex)==NotchAmp);
            DeltaAmp(t,1)=(FirstAmp-NotchAmp)/FirstAmp;
            DeltaAmp(t,2)=MaxIndex;
            DeltaAmp(t,3)=NotchIndex;
        end

        prefMod=find(DeltaAmp(:,1)==max(DeltaAmp(:,1)));
        if length(prefMod)>1
            OSI(c,:)=NaN;
        else
            if length(thetas)==6
                if prefMod==1
                    orthMod=4;
                elseif prefMod==2
                    orthMod=5;
                elseif prefMod==3
                    orthMod=6;
                elseif prefMod==4
                    orthMod=1;
                elseif prefMod==5
                    orthMod=2; 
                elseif prefMod==6
                    orthMod=3;
                end
            elseif length(thetas)==4
                if prefMod==1
                    orthMod=3;
                elseif prefMod==2
                    orthMod=4;
                elseif prefMod==3
                    orthMod=1;
                elseif prefMod==4
                    orthMod=2;
                end
            end

            prefModIndex1=DeltaAmp(prefMod,1);
            prefModIndex2=DeltaAmp(prefMod,2);
            prefModIndex3=DeltaAmp(prefMod,3);
            orthModIndex1=DeltaAmp(orthMod,1);
            orthModIndex2=DeltaAmp(orthMod,2);
            orthModIndex3=DeltaAmp(orthMod,3);

            OSI(c,:)=(prefModIndex1-orthModIndex1)/prefModIndex1;
        end

%         for f=1:length(thetas)
%             OrientationLabel=strcat('Theta',num2str(thetas(f)),'deg');
%             MaxAmp(f)=WeakCone.(matlab.lang.makeValidName(OrientationLabel)).Amplitude(prefModIndex1(1));
%             NotchAmp(f)=WeakCone.(matlab.lang.makeValidName(OrientationLabel)).Amplitude(prefModIndex2(1));
%         end

%         NotchMod=max(NotchAmp)-min(NotchAmp);
%         OSI(c,:)=NotchMod/(mean(NotchAmp));
    end
    %clearvars -except Data Cells thetas OSI;
end