function [ chem, u,s] = loadAddProp( chem, u, s,par )
%loadAddProp Parse information that can only be done when all the reports
%are loaded
%   Some variables require post processing that can only be done when
%   everything is loaded. This includes dealing with units, and
%   designtaions.
%   The main addition is the change in conditions across a unit.
%   dim(u): n->(n,# of different possible changes)
%   Also contains sizing calculations for the units in a function below

findLine = @(c , string) (find(~cellfun('isempty', strfind(c,string))));
dH=@(T,ref,CP)(ref+CP*(T-25)); 
untLen=size(u,1);
parLen=length(par);
chmLen=length(chem);

u(1).chemNames={chem.Name};

% First make sure the streams are nice or estimate stuff

for w = 1:length(s)
    chemNums=s(w).compList;
    
    if estVar(s(w),'MW')
        s(w).MW= s(w).x(1,chemNums) * [chem(chemNums).MW]';
    end
    
    if estVar(s(w),'CP')
        if s(w).Phase(1)>0
            vr='CPG';
        else
            vr='CPL';
        end
        
        temp=0;
        for ch=chemNums
            temp = s(w).x(1,ch)*chem(ch).(vr)(s(w).T+273.15)+temp;
        end
        s(w).CP=temp;
        

    end
    
    if estVar(s(w),'viscous')
        s(w).viscous= s(w).x(1,chemNums) * [chem(chemNums).muv]';
    end
    s(w).flagPsu=0;
end


% create the changes within a unit using unit alternatives

divMean = {{'T','rho','CP'},'P',{'x','Activity','rho'},'','','',''};
possChg = {'T','rho','CP','P','x','Activity'};

%     if min(u(n).class{2}>3)

for n = 1:untLen
    
    %Make a combined stream for all units with multiple inlets
    
    if u(n).multIn
        in=length(s)+1;
        u(n).inPsu= in;
        s(in)=crossStream(s,u(n).in);
        
    else
        in=u(n).in;
    end
    
    u(n).Pin=s(in).P;
    
    %%%%%%%%%
    %Change the failure rates
    %%%%%%%%%
    
    %max fail for base case
    
    if u(n).P>2
        highP=5*10^-6;
    else
        highP=100*10^-6;
    end

    if u(n).classFunc==4
        highT=50*10^-6;
        highP=5*10^-6;
    else
        highT=100*10^-6;
    end
    
    %adjust compressors base fail rate
    if u(n).classFunc==3 && u(n).flagCompr
        if u(n).power<1000
            u(n).fail=550.66*10^-6;
        elseif u(n).power<3000
            u(n).fail=880.21*10^-6;
        else
            u(n).fail=2433.02*10^-6;
        end
    end
    
    u(n).failbaseTM=highT;
    u(n).failbasePM=highP;

    % Adjust for the model alterations
%     u(n).failadj=u(n).failbase*exp(a*u(n).failbase+b); 
    detectBase=[.2;.4;.3;.1;.5;.2;.3;.2;.3;.3;.7;.3;.3]; %composition or amount gets a .3
    
    
    if u(n).altFlag>0
        
        %Make the designation string
        dispDes{1,1}=['Unit ' num2str(n) '(' u(n).name ')'];
        dispDes{2,1}=[u(n).type 'T=' num2str(u(n).T) 'C , P=' num2str(u(n).P) 'atm'];
        dispDes{3,1}=u(n).Indes;
        dispDes{4,1}=u(n).Outdes;
        
        %make the GUI and get the user input
        ui = loadAlterations(u(n).altFlag,dispDes,par);
        
        %intialize the base armor
        armor{1}=zeros(2,1); %[max T; maxP]
        armor{2}=0; %base fail rate
        armor{3}=zeros(parLen,1); %[detection(paras)]
        armor{4}=zeros(parLen,1); %[consequences(paras)]
        
        %go through all the alterations and add up effects
        for a = 1: u(n).altFlag
            for d = 1:length(ui{1,a})
                if ui{1,a}==1
                    temp=length(ui{3,a});
                else
                    temp=1;
                end
                
                for k = 1:temp
                    
                    armor{ui{1,a}(d)}(ui{3,a}(k)) = armor{ui{1,a}(d)}(ui{3,a}(k))...
                        +(ui{4,a}*ui{5,a}*ui{6,a});
                end
            end
            
            

        end
        
        


        u(n).failadj = (1-indize(armor{2}))*u(n).failbase;
        
            
        u(n).failadjTM = (1-indize(armor{1}(1)))*highT;
        u(n).failadjPM = (1-indize(armor{1}(2)))*highP;
        
        for p=1:parLen
            
            u(n).detect(p,1) = detectBase(p)-indize(armor{3}(p))/2;
            
            u(n).armorConq(p,1) = (indize(armor{4}(p))/2);
        end
        
        
        
%         [u(n).failadj,u(n).failadjTM, u(n).failadjPM] = ...
%             loadAlterations(u(n).flagAlt,dispDes,u(n).failbase,u(n).failbaseTM, u(n).failbasePM);
    else
        %If there is no alterations just use the base variables
        u(n).failadj=u(n).failbase;
        u(n).failadjTM=highT;
        u(n).failadjPM=highP;
        
        u(n).detect = detectBase;
        u(n).armorConq =zeros(parLen,1);
    end

    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % do the different conditions within the unit
    %%%%%%%%%%%%%%%%%%%%%%%
    
    base=find(u(n).chg(:,1),1);
    %don't want to average every field just those that change
    if base < 0
        base = u(n).out(1);
    end
        
        
    for ch=1:u(n).numCon
        
        if ch~=1
            u(n,ch)=u(n);
        end
        
        nc = sub2ind(size(u),n,ch);
        u(nc).des=u(n).maDes{ch};
    
        % Assign the variables that potentially can change
        
        [u] = assignVar(possChg,u(n).chg(:,ch) ,u,nc,s);
        
        %Do the special cases that require some attention
        
%         u(nc).fail =
%         augfailTP(u(nc).T,u(nc).P,u(nc).failadj,1*10^-6,.1*10^-6,u(n).classFunc) ;
        chBase = find(u(n).chg(:,ch)~=base,1);
        if isempty(chBase) || chBase < 0 
            chBase = base;
        end
        u(nc).dH=s(chBase).H-s(in).H;
       
        
        u(nc).Pout = u(nc).P;
        if  ch~=1 && u(n).classFunc==4 && estVar(u(n), 'dHreac')
            
            temp=(u(nc).dH+u(nc).duty)/(u(nc).F(1)*u(nc).conv);
            temp=4184*temp*3600;
            u(nc).dHreac=temp;
        end
        
        % Assign the rest of the variables
        
        potFields = fieldnames(s(base));
%         extFields = fieldnames(u(nc));
        exclude = {'labels','from','to'};
        lenPot=length(potFields);
        Bpot = ones(lenPot,1).*base;
        
        for p = 1:lenPot
            var = potFields{p};
%             temp=find(strcmp(extFields, var),1);
            
            
            if (~estVar(u(nc), var)) ...
                    || ~isempty(find(strcmp(exclude,var),1))
                Bpot(p)=0;
            end
        end
        
        [u] = assignVar(potFields,Bpot ,u,nc,s);
        
    %Miscellaneous changes to the units as well as sizing.
    
           
    %change the unit chemical list to represent the new conditions.
        %if its below 10-10 its not a chemical in the unit
    u(n).compList = find(u(n).x(1,:) > 10^-8);
    u(n).numComp = length(u(n).compList);
        
    if ch==1
        %do sizing calculations the first time through


        u(nc).Fv=u(nc).F(1)/(u(nc).rho(1)*60);
%         u(n).Fv=pmass{2}*u(n).rho(1)*(1000/60);
        u(n).F(3)=u(n).Fv;
        u(n).unitF{3}='L/MIN';
        
        %Now size the first unit then assign it to all the alternatives

        u = sizing(u,n,s);
    end
    end
end



   

end

function [out] = indize(val)
if val ==0
    out=0;
else
    
    temp=sort([-1,1,(1/log10(25))*log10(val)]);
    out=temp(2);
end

end

function [u] = assignVar(vars, streams,u,n,s)

if length(streams)==1
    temp=ones(length(vars),1);
    streams=temp*streams;
end

for p=1:length(vars)
    w=streams(p);
    if w~=0
        if w > 0
            u(n).(vars{p})=s(w).(vars{p});
        else
            if w == -888
                %average two outlets
                fir = u(n).out(1);
                sec = u(n).out(2);
            else
                %average outlet and inlet
                fir = u(n).out(1);
                sec = u(n).in(1);
            end
            u(n).(vars{p}) = (s(fir).(vars{p})+s(sec).(vars{p}))/2;
        end 
    end
end

end

function mix = crossStream(s,comb)

mix = s(comb(1));
mix.flagPsu=1;

% mix.CP=([s(comb).F(2)]')*[s(comb).CP(1)];

mix.T=mix.H(2)/mix.CP;


tP=0;
tCP=0;
tMW=0;
tx1=zeros(1,length(s(comb(1)).x));
tx2=tx1;
ta=tx1;
tr=[0,0];
tF=0;
tH=0;
for si=comb
    tF=tF+s(si).F;
    tH= tH+s(si).H;
end
mix.F=tF;
mix.H= tH;

xF= @(w,type)(s(w).F(type)/mix.F(type));
for si=comb
    tP= tP+s(si).P*xF(si,3);
    tx1= tx1+s(si).x(1,:).*xF(si,1);
    tx2=tx2+s(si).x(2,:).*xF(si,2);
    ta(1,:)= ta+s(si).Activity.*xF(si,1);
    tr(1)= tr(1)+s(si).rho(1).*xF(si,1);
    tr(2)= tr(2)+s(si).rho(2).*xF(si,2);
    tMW = tMW+s(si).MW*xF(si,1);
    tCP = tCP + s(si).CP*xF(si,2);
    
%     
%     tP= tP+s(si).P*s(si).F(3);
%     tH= tH+s(si).H.*s(si).F;
%     tx(1)= tx+s(si).x.*s(si).F(1);
%     tx(2)=tx+s(si).x.*s(si).F(2);
    
%     for type=1:3
%         tH(type)= tH(type)+s(si).H(type)*s(si).F(type);
%     end

    
end


mix.P= tP;

mix.CP= tCP;
mix.MW= tMW;
mix.x= [tx1;tx2];
mix.Activity= ta;
mix.rho= tr;
mix.T=mix.H(1)/mix.CP;


% for so=sot
%     tPo=tPo+s(so).P*s(so).F(3);
% end


end
            

function [u ] = sizing(u,n,s)

if u(n).multIn
    sin=u(n).inPsu;
else
    sin=u(n).in;
end
sot=u(n).out;

switch u(n).classFunc
case 1
    %mixers , tanks
    u(n).tau = 0.25; %15 min residence time

    tempA = (u(n).F*u(n).tau);



    u(n).Vvapor=(tempA(1)/u(n).rho(1))*.25; % L
    u(n).V = (tempA(1)/u(n).rho(1)) + u(n).Vvapor; %add vapor space in sizing

    u(n).Amount = tempA;

case 2
    %separtions
    vap=s(u(n).outTop);
    liq=s(u(n).outBot);
    area=(1/3.28^2)*0.21*vap.F(1)*(u(n).MW(1) /(u(n).rho(1)*28316.85))^.5;

    if ~estVar(u,'stages')
        %2 feet per stage then 15% of the stages for the space at
        %top and bottom

        H=2.3*u(n).stages; 
        u(n).V=area*H;
    else
        %flash so use 5 feet (3+1+diameter of feed) + the liquid
        %pool which has 2 min of residence time.

        H=5/3.28+(vap.F(3)*2*1000/(area));

        u(n).V=area*H*1000;%L
        u(n).Vvapor=u(n).V-(vap.F(3)*2*1000);
    end

    u(n).Amount=u(n).rho.*(u(n).V*1000)*1000;

%             u(nc).Amount =[50;500]; %DO THIS
%             u(nc).V = [1000]; %DO THIS
%             u(nc).Vvapor= 1000; %Do this
case 3
    %Heaters / compressors

    if ~u(n).flagCompr
        diff(1)=(s(sin).T+273.15)-(298.15);
        if u(n).duty>0
            diff(2)=(s(sot).T+273.15)-(30+273.15);
        else
            diff(2)=(s(sot).T+273.15)-(15+273.15);
        end
        diff=abs(diff);
        u(n).dTLM=(diff(1)-diff(2))/(log(diff(1)/diff(2)));
        %Q=UAT with U coming from the average for organic and steam
        %/water

        u(n).SurArea=4.184*abs(u(n).duty)/((100*5.673)*u(n).dTLM); 
        %pipe W is around 7.5*L
        u(n).V = 1000*(pi*(7.5/4)*(u(n).SurArea/(pi*7.5))^(1.5));
        u(n).Vvapor= 0;
        u(n).tau=(u(n).V/u(n).F(3))/60;


        u(n).Amount(2)= u(n).V*u(n).rho(1);
        u(n).Amount(1)= u(n).Amount(2)/u(n).MW(1);
    else
        %Compressor
        
        u(n).tau = 5/36000; %5 seconds

        tempA = (u(n).F*u(n).tau);
        u(n).Vvapor=(tempA(1)/u(n).rho(1))*.25; % L
        u(n).V = (tempA(1)/u(n).rho(1)) + u(n).Vvapor; %add vapor space in sizing

        u(n).Amount = tempA;

        exp=u(n).gamma/(u(n).gamma-1);
        u(n).surge = (u(n).Pout/u(n).Pin)*(((u(n).effIsen*u(n).dH(1))/(u(n).CP(1)*(u(n).T+273.15))) +1)^(exp);

    end

case 4 
    %Reactors
    

    tin=0;
    tout=0;
    tx=[];
    for in=sin
        tin=tin+s(in).F(1);


        tx(in,:)=(s(in).F(1).*s(in).x(1,:));
    end
    Ftin=tin;

    Fin=sum(tx,1);
    xin=Fin./Ftin;
    for out=sot
        tout=tout+s(out).F(1);
        temp(out,:)=s(out).F(1).*s(out).x(1,:);
    end
    Ftout=tout;
    Fot=sum(temp,1);
    xot=Fot./Ftout;


%         diff=(s(sot).F(1).*s(sot).x(1,:))-(s(sin).F(1).*s(sin).x(1,:));
    diff=Fot-Fin;
    if ~estVar(u,'stoCel')
        stoic=zeros(chmLen,1);
        for r=1:length(u(n).stoCel)
            temp(:,r)=getChem(chem,u(n).stoCel{r});

            lim(1,r)=findLine({chem.ID},u(nc).unitconv);
            temp(:,r)=temp./abs(temp(lim(r)));
        end
        stoic(:,1)=sum(temp,2);

        
    else
        %get the stoic matrix through other means
        %glhf

%             if ~estVar(u(n), 'inert')
%                 temp=getChem(chem,u(n).inert);
%             end
        
        [junk,lim]=min(diff);
        diff =abs(diff); %moles reacted
        stoic(:,1)=diff./abs(diff(lim));
        u(n).conv =diff(lim)/Fin(lim);
        phi = Fot./Fot(lim);

        


    end
    
    
    if estVar(u, 'dHreac')
        reacEng=(u(n).dH+u(n).duty)/(u(n).F(1)*u(n).conv);
        reacEng=4184*reacEng*3600;
        u(n).dHreac=reacEng;
    end

    SA=abs(u(n).duty)/265000; %min heat out through a hexchanger
    
    u(n).V=SA*0.05*1000; %ID of .1m volume in L
    u(n).Vvapor=u(n).Phase(1) * u(n).V;
    u(n).Amount(2)=u(n).V*u(n).rho(1);
    u(n).Amount(1)=u(n).Amount(2)*u(n).MW(1);
    
    

case 5
    %Transport
    u(n).tau = 1/3600; %1 second

    tempA = (u(n).F*u(n).tau);



    u(n).Vvapor=(tempA(1)/u(n).rho(1))*.25; % L
    u(n).V = (tempA(1)/u(n).rho(1)) + u(n).Vvapor; %add vapor space in sizing

    u(n).Amount = tempA;


end
       


end


function [u ] = assignUnit(val,var,u,s,w)


n= s(w).from;
numChg= size(u,2);
if n~=0
if u(n).chgFlag
    %Appropiately assign stream val to the right units that change
    ch=1;
    
    while (ch<=numChg && ~estVar(u(n,ch),'name'))
        %Test to see if the stream is right
        if u(n).multOut
            if u(n).outTop==w
                str='Top';
            else
                str='Bot';
            end
        else
            str='out';
        end

        req{1,1}=str;
        req{1,2}= 1;

        if strcmp(var,'CP')
            req{2,1} = 'T';
            if strcmp(str,'Top')
                req{2,2}= 1;
            else
                req{2,2}= 0;
            end
        end
        
        %Take the requirements and use it on the designation of the unit
        flag=1;
        r=1;
        while (flag==1 && r<=size(req,1))
            flag=~isempty(strfind(u(n,ch).des,req{r,1}));
            if ~req{r,2}
                flag=~flag;
            end
            r=r+1;
        end
        
        %ifthe unit does use the stream then set it equal to the var
        if flag
            u(n,ch).(var)=val;
        end
        
        ch=ch+1;
    end
    
    
else
    %just assign the stream val to the units 1:1
    if estVar(u(n),var)
        u(n).(var)=val;
    end
    
    
end
end
end


function  [mat] = getChem(chem,list)
m=size(list,1);
mat=zeros(m,1);
for fl= 1:m
    c= ~cellfun('isempty', strfind({chem.ID},list{fl,1}));
    mat(c,1)=list{fl,2};
end

end