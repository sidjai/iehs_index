function [msg, allclear,u] = makeUserInput(u,chem, cfg)
%makeUserInput Partialy fill in the Safety info and unit list to make the
%filling them in easier
%   flag is a vector that says if the safetyfile is filled in or not and if
%   the unit list is filled in or not.

findLine = @(c , string) (find(~cellfun('isempty', strfind(c,string)),1));

chmLen= length(chem);
untLen= size(u,1);
cfg.parametersarLen= length(cfg.parameters);
name{1}= cfg.InfoWkstFile;
name{2}= cfg.UnitAlterations;
name{3}= cfg.weightFile;

des{1}='Safety info wkst';
des{2}='Unit Alterations wkst';
des{3}='User weight wkst';
%%%%%%%%%%%%%%
%Check if the files are there and filled 
%%%%%%%%%%%%%%

for i=1:3
    info=dir(name{i});
    flag(i)=isempty(info);
    
    
    
    if flag(i)
        filled(i)=0;
    else
        %look into the file and see if it hasn't been filled
        [nums,junk,raw]= xlsread(name{i},1);
        switch i
          case 1 
            %Criteria for filled: Safety info
                %whether the boiling point is filled
            [junk indVal junk indField]=findColNums(raw);
            inTb=findField(raw,'Tb',indField);
            temp = raw{inTb,indVal}; 
            if (isempty(temp) ||  isnan(temp))
                filled(1)=0;
            else
                filled(1)=1;
            end

                
        
           
            
          case 2
            %Criteria for filled: Unit Alterations
                %If anything is turned on or off
            boo=nums(:,end);
            
            if cfg.altFlag
                filled(2)=max(boo);
            else
                filled(2)=1;
            end
            
            %Also parse the data
            if filled(2) && cfg.altFlag
                for n=1:size(u,1)
                    u(n).altFlag=boo(n);
                end
            else
                for n=1:size(u,1)
                    u(n).altFlag=0;
                end
            end
            
           
            
          case 3 
            %Unit weights
            %criteria: if anything is not =1
            for ws=1:3
                if ws ==1
                    rawWeight{1}=nums;
                else
                    rawWeight{ws}=xlsread(name{i},ws);
                end
                fillTemp(ws,1)=~isempty(find(rawWeight{ws}~=1,1));
            end
            if cfg.weightFlag
                
                filled(3)=max(fillTemp);
            else
                filled(3)=1;
            end
            
            %also parse the weights
            if filled(3) && cfg.weightFlag
                %ordered by combination steps
                %1)SecEff add, EDP; 2)Agro chems, EDPAchems; 3)Mix add, MDP
                %4)Vessel add, CDP; 6)Agro units, CDPAunits; 7)UUs add, IDP
                weight = cell(6,1);
                
                if fillTemp(1)
                    weight{2}=numNonz(rawWeight{1});
%                     unitweight{2}='Chem';
                    
                end
                
                if fillTemp(2)
                    for type=1:4
                        if ~isempty(find(rawWeight{2}(:,type)~=1,1))
                            avgfied = numNonz(rawWeight{2});
                            switch type
                                case 1
                                    
                                    weight{1}=avgfied(:,1);
%                                     unitweight{1}='Add SecEffs';
                                case 2
                                    weight{3}=avgfied(:,2);
                                case 3
                                    weight{4}=avgfied(:,3);
                                case 4
                                    weight{6}=avgfied(:,4);
                            end
                        end
                        
                    end
                    
                    
                end
                    
                if fillTemp(3)
                    weight{5}=numNonz(rawWeight{3});
%                     unitweight{5}='Unit';
                end
                
                
                u(1).userWeights=weight;
                
            else
                u(1).userWeights=0;
            end
                
            
           
        end
        
    end
    
    if i==2 && cfg.altFlag==0
        flag(2)=0;
        filled(2)=1;
    end
    
    if flag(i)
        msg{i,1} = ['A new ' des{i} ' has been created and partially filled' char(10) 'located at ' name{i} char(10) 'go fill in the rest before proceeding'];
    else
        if ~filled(i)
            msg{i,1}=['You have not filled the ' des{i} char(10) 'located at ' name{i} char(10) 'Go fill it before proceeding'];
        else
            msg{i,1}=[des{i} ' is filled and ready to go'];
        end
    end
    
                
           
end

% Find out if all files are ready to go.
allclear= 0;
if ~max(flag) && min(filled)
    allclear=1;
%     msg='All clear: Safety file and Alt wkst are ready to go';
end
    
%check to see if there is a database list

check = dir(cfg.MasterChemListName);
prevFlag=~isempty(check);


if prevFlag
    fid = fopen(cfg.MasterChemListName);
    master=textscan(fid,'%s %s %d','Delimiter', ',');
    fclose(fid);
    %Check the database list to see if the chem is there
    for c=1:chmLen
        temp = findLine(master{1},chem(c).Formula);
        if isempty(temp)
            ind(c)=0;
        else
            ind(c)=temp;
        end
    end
else
    ind=zeros(1,chmLen);
end
    

    
%Make the sheets that need to be made

%%%%%%%%%%%%%%%%%%%%%
%Make the safety info sheet, from templete or database
%%%%%%%%%%%%%%%%%%%%%
if flag(1)
%     [junk junk raw]=xlsread('safetyInfoTemplate.xltx');
    
    [indVar indVal indUnit indField]=findColNums(cfg.IOnotes{4});

    for c=1:chmLen
        
       

        %Check the database list to see if it is already done
        
        if ind(c)>0
            Csheet{c}=xlsread(master{2}{ind(c)},master{3}{ind(c)});
        else
            %Use the templete to make the sheet with some filled in info
            Csheet{c}=cfg.IOnotes{4};
            labels{c}=chem(c).ID;
            ind=findField(cfg.IOnotes{4},'Name',indField);
            Csheet{c}{ind,indVal}=chem(c).ID;

            ind=findField(cfg.IOnotes{4},'formula',indField);
            Csheet{c}{ind,indVal}=chem(c).Formula;
        end
        


    end

    for c=1:chmLen
        if c<4
            xlswrite(cfg.InfoWkstFile,Csheet{c},c)
        else
            xlswrite(cfg.InfoWkstFile,Csheet{c},labels{c})
        end
    end


elseif filled(1)
    
    %add the safety sheets to the database list
    if prevFlag
        fid = fopen(cfg.MasterChemListName,'at');
        fseek(fid,0,'eof');
    else
        fid = fopen(cfg.MasterChemListName,'wt');
    end
    
    for c=1:length(chem)
        if ind(c)==0
%             safe=strrep(safetyfile, '\','\\');
            fprintf(fid,'%s \r\n',[chem(c).Formula ', ' cfg.InfoWkstFile ', ' num2str(c)]);
        end
    end
    fclose(fid);
    
end


%%%%%%%%%%%%%%%%%%%
%Make the Unit Alterations sheet
%%%%%%%%%%%%%%%%%%
if flag(2)
    

    %header
    
    Usheet={'Number','Name','Type','Input','','Output','','','Any Alterations?'};
%     Usheet{2,1} ='Whole design';
%     Usheet{2,7} =1;
    add=1;
    for n=1:untLen
        Usheet{n+add,1}=n;
        Usheet{n+add,2}=u(n).name;
        Usheet{n+add,3}=u(n).type;
        Usheet{n+add,4}=strrep(u(n).Indes,'Input: ','');
        Usheet{n+add,6}=strrep(u(n).Outdes,'Output: ','');
        
        Usheet{n+add,9}=0;
        
        



    end


    xlswrite(cfg.unitAlterations,Usheet,1)
end


%%%%%%%%%%%%%%%%%%%
%Make the weight wkst
%%%%%%%%%%%%%%%%%%
if flag(3)
    %4 sheets
    
    %1)chem weights
    Wsheet{1}=[{''};cfg.parameters'];
    
    for c=1:chmLen
        Wsheet{1}{1,1+c}=chem(c).ID;
        %Fill it with the initial weights to be used (1)
        for p=2:parLen+1
            
            Wsheet{1}{p,1+c}=1;
        end
    end
    
    %2)agro step for all the methods
    Wsheet{2}=[{''},p',{ones(parLen,3)}];
    Meth ={'SecEff','Mix','Vessel'};
    for ag=1:3
        Wsheet{2}{1,1+ag} = Meth{ag};
        for p=2:parLen+1
            
            Wsheet{2}{p,1+ag}=1;
        end
    end
    
    %3)unit weights
    Wsheet{3}=[{''};p'];
    
    for n=1:untLen
        Wsheet{3}{1,1+n}=u(n).name;
        %Fill it with the initial weights to be used (1)
        for p=2:parLen+1
            
            Wsheet{3}{p,1+n}=1;
        end
    end
    
    %4)overall control??? (mass weight etc)
    
        
    for ws =1:length(Wsheet)
        xlswrite(cfg.weightFile,Wsheet,ws)
    end
    
    
end


end

function [val] = numNonz(vec)

for c = 1:size(vec,2)
    tm = find(vec(:,c));
    val(:,c)=vec/length(tm);
end

end

function [indVar indVal indUnit indField] = findColNums(raw)
for h=1:size(raw,2)
    if (ischar(raw{1,h}) || ~isnan(raw{1,h}))
        if ~isempty(strfind(raw{1,h},'ble'))
            indVar=h;
        elseif ~isempty(strfind(raw{1,h},'alue'))
            indVal=h;
        elseif ~isempty(strfind(raw{1,h},'nit'))
            indUnit=h;
        elseif ~isempty(strfind(raw{1,h},'ield'))
            indField=h;
        end
    end
end


end

function [val] = findField(raw, field, indField )

ind=2;
while ind< size(raw,1)
    temp=strcmp(field,raw{ind,indField});
    if temp
        val=ind;
        ind=999;
    else
        ind=ind+1;
    end
end
    


end