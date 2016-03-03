function [ u, chem ] = calcIntVal( u, chem,n, ch, parameters, justify)
%Does the calculations to get the EHS indices for one unit
%   At this point plant is loaded up with the aspen and the safety
%   information which will be used to calculate the index
%   n = unit number
%   ch = the number indicating which internal change in variables should be
%   used. 1 is the default steady state / outlet set of variables and it
%   changes as you go higher or across the columns in the large matrix.


 %if the variable changes and you want to figure out the top vals
 %, change the linear index, n, to the new index with the second u

     n = sub2ind(size(u),n,ch);

chemNums=u(n).compList;
chmLen=u(n).numComp;
sep = u(n).multOut;

T=u(n).T;
Tk=T+273.15;
P=u(n).P;
% chem=plant.chem;
temp=0;
tempw=0;
temptb=0;

%set up the CRW values for reactive groups
noncombust = [2;44;100;37;98;38;39];
oxidizers = [2;44];
organics = [33;47;17;27;32;42;20;30];



simple = @(c)(find(chemNums==c,1));

%%%
%Calculate additional physical parameters
%%%
k=1;
u(n).vapAmount =[];
for c=chemNums
    u(n).psat(c,1)=chem(c).Psat(Tk);
    u(n).ppar(c,1)=(P)*u(n).x(1,c)*u(n).Activity(1,c);
    u(n).Kval(c,1)=u(n).psat(c)*u(n).Activity(1,c);
    u(n).vapAmount(k,1)=u(n).ppar(c,1)*u(n).Vvapor/(.08206* Tk);
    temp=temp+u(n).ppar(c,1);
    temptb=temptb+chem(c).Tb;


    %Add up all the combustion energies for the tank
    tempw= tempw+u(n).x(1,c)*u(n).Amount(1)*chem(c).dHburn;

    k=k+1;
end
avgTb=temptb/chmLen;
u(n).combEner=tempw;
u(n).ytheo=u(n).ppar/temp; %theoretical composition of the vapor in eq at this temp
yt= [u(n).ytheo(chemNums)];
u(n).relExp{1}=u(n).vapAmount/100;

u(n).relExp{2}=u(n).vapAmount/100*size(u,1);
temp=0;



%mixture bubble point
% tol=.05;
% result=0;
% m=1;
% while (abs(1-result)>tol) && m<50
%     if result==0
%         Tguess=avgTb;
%     else
%         Tguess=((1-result)*Tguess)+Tguess;
%     end
%
%     result=0;
%     for cc=chemNums
%
%         tek=chem(c).Psat(Tguess)/P;
%         result=result+tek*u(n).x(1,c);
%     end
%     m=m+1;
% end
% u(n).Tbmix=Tguess;


%take out noncombustible and calc from the ytheo
%also see if there is any organics or oxidizers

k=0;
ycomb=zeros(length(chem),1);
for c=chemNums;

    if chem(c).flagFuel
        ycomb(c,1)=u(n).ytheo(c);
    else
        ycomb(c,1)=0;
        k=k+1;
    end
end
    ycomb=ycomb/(chmLen-k);





%%%%%%%%%%%%%%%%%%%%
% Individual Values
%%%%%%%%%%%%%%%%%%%%

k=1;
for c=chemNums

%mobility
pa=1;
    if u(n).Phase(1)>0
        %gas dispersion

        potVal{pa,1}(1,k)=chVal(.5+((18/chem(c).MW(1))^(.5))...
            *(1.18/chem(c).rho(1))*(abs(chem(c).bouy))^...
            (.5*sign(chem(c).bouy)));

    end

    if u(n).Phase(1)<1
        %liquid or two phase flow release
        fv=real(u(n).CP(1)*(T-chem(c).Tb)/chem(c).dHvap(Tk));
        if (P>1 && u(n).T>chem(c).Tb)
            %flashing liquid
%             u(n).value(1,k)=.4*(1+fv);
            potVal{pa,1}(2,k)=.4*(1+fv);
        else
            %vaporizing liquid (mass transfer - get flux)
            %average of Index(K(psat-p))
%             potVal{pa,1}(2,k)=min(.4+.004*((18/chem(c).MW)^(1/3)*(u(n).psat(c)-u(n).ppar(c))/P),0.8);
            potVal{pa,1}(3,k)=chVal(min(0.6,0.2 + (.4*u(n).Kval(c)/50)...
                *(1/u(n).viscous(end))));

        end
    else
        %gas release so use chocked flow
%         if ~isempty(strfind(u(n).P{2},'psi'))
%
%         end

%         potVal{pa,1}(2,k)=chVal(.8+(u(n).ppar(c)*.02));
    end
    if chem(c).hflag==1
        potVal{pa,1}(3,k)=chem(c).hval(pa);
    end
%     value(1,k)=chkLowLimit(potVal{1});
%fire and explosion
pa=2;
    %Fuel
    if ~estVar(chem(c),'Flash')
        potVal{pa,1}(1,k)=chVal(1-(chem(c).Flash-T)*.5/100);
    end
        %Predict flash point from boiling point using Rao eq 6-1
        %need to seperate by chemical clas
        %then use the table to get the constants.
    if ~estVar(chem(c),'FlashPred')
%         isfield(u(n),'FlashPred') && c<length(u(n).FlashPred) && u(n).FlashPred(c)~=0
        potVal{pa,1}(3,k)=chVal(1-(chem(c).FlashPred-T)*.5/100);
    else
        chem(c).FlashPred=((chem(c).Tb+273.15)*(.726-.0000715*...
            (chem(c).dHburn*1000/chem(c).dHvap(300))))-273.15;

        potVal{pa,1}(3,k)=chVal(1-(chem(c).FlashPred-T)*.5/100);

    end

    Cst=100/(1+(chem(c).oxyStoic/.21));
    if estVar(chem(c),'LFL')
        chem(c).LFL=.55*Cst;
    end

    if estVar(chem(c),'UFL')
        chem(c).UFL=3.5*Cst;
    end

    LFLtemp=chem(c).LFL-(0.75/(chem(c).dHburn*(10^-6)/4.184))*(T-25);


        potVal{pa,1}(2,k)= chVal(1-LFLtemp/15);
    %Mix-



    %Oxygen
    if ~estVar(chem(c),'LOC')
        potVal{pa,1}(3,k)=chVal(1-(chem(c).LOC-5)/25);
    else
        if ~estVar(chem(c),'LFL')
        chem(c).LOC=chem(c).LFL*(chem(c).oxyStoic);
        potVal{pa,1}(3,k)=chVal(1-(chem(c).LOC-5)/25);
        %or use 6-16 if you know UFL
        %chem(u).LOC=((LFL+1.11*UFL)/2.11)*(.21*(100-UFL))/UFL
        end

    end

    if chem(c).hflag==1
        potVal{pa,1}(4,k)=chem(c).hval(pa);
    end

%     potVal{pa,1}(:,k)=[mean(valFuel);mean(valOxy)];
%     value(2,k)=chkLowLimit([valFuel,valOxy]);

%Acute Toxicity

pa=3;
    % Use TWA-C instead
    if ~estVar(chem(c),'TLVc')
%         potVal{3}(1,k)= 1-.25*(log10(chem(c).ERPG)-1);
        potVal{pa,1}(1,k)=chVal(1-log10(chem(c).TLVc/10)/4);
    end
    if ~estVar(chem(c),'IDLH')
        potVal{pa,1}(2,k)= chVal(1-(log10(chem(c).IDLH/10))/4);
    end
    if ~estVar(chem(c),'LCInhal')
        potVal{pa,1}(3,k)=chVal(1-.5*log10(chem(c).LCInhal/100));
%         potVal{pa,1}(2,k)=.6560-.3212*log10(chem(c).LD);
    end

    if chem(c).hflag==1
        potVal{pa,1}(4,k)=chem(c).hval(pa);
    end

%Prob Runaway Decomp

pa=4;

    if ~estVar(chem(c),'Tdecomp')
      potVal{pa,1}(2,k)=chVal(1-(chem(c).Tdecomp-T)/200);
    end

    if ~estVar(chem(c),'AIT')
        potVal{pa,1}(4,k)=chVal(1-(chem(c).AIT-T)/200);
    elseif chem(c).hflag==1
        potVal{pa,1}(3,k)=chem(c).hval(pa);
    end

%     varChoice='';
    dHInd=[0,0,0,0];

    energy={'dHreac','dHdecomp','dHburn','dGform'};
    for er=1:length(energy)
        if ~estVar(chem(c),energy{er})
            dHInd(er)=1;
        end
    end

    if estVar(u(n),'dTad') || c>length(u(n).dTad) || u(n).dTad(c)==0

        vInd= find(dHInd,1);
        if ~isempty(vInd)
%             u(n).dTad(1,c)=-chem(c).(energy{vInd})*u(n).x(1,c)*u(n).Amount(1)...
%                 /((u(n).CP(1)/100)*u(n).Amount(2));
            u(n).dTad(1,c)=TNTeq(chem,u,c,n,energy{vInd})...
                /((u(n).CP(1)/1000)*u(n).Amount(1));
        else
            u(n).dTad(1,c)=10^-10;
        end
    end

    asdf(1)=chVal((1/.6)*log10(abs(u(n).dTad(1,c))/25));

    %Self reactivity CRW

    asdf(2)=chem(1).ReacMat(c,c)+chem(1).GasMat(c,c);


    [potVal{pa,1}(1,k),jhk]=max(asdf);

    flagMax(pa,1)=jhk-1;

        % Estimate AIT maybe?
%         potVal{pa,1}(1,k)=valAdTime;


%Expected damage
pa=5;


    for hi=find(dHInd)
        potVal{pa,1}(hi,k)=chVal((1/50)*...
            ProbTNT(TNTeq(chem,u,c,n,energy{hi})/300,15));
    end

    if chem(c).hflag==1
        potVal{pa,1}(5,k)=chem(c).hval(pa);
    end


%Explosion
pa=6;
    %If contains all C,H,N,O, Use Oxygen content

     if ~estVar(chem(c),'velsu')

         potVal{pa,1}(1,k)=chVal((chem(c).velsu/50)-.9);
     end


     if ~estVar(chem(c),'dHexp')
         potVal{pa,1}(2,k)=chVal((1/50)*...
            ProbTNT(TNTeq(chem,u,c,n,'dHexp')/300,15));
     end


    if ~estVar(chem(c),'oxyCont') && chem(c).contCHNO==1
        if chem(c).oxyCont(1)>0
            potVal{pa,1}(4,k)=2-chem(c).oxyCont/80;
        else
            potVal{pa,1}(4,k)=2+chem(c).oxyCont/120;
        end
    end


    if chem(c).hflag==1
        potVal{pa,1}(5,k)=chem(c).hval(pa);
    end

    if ~estVar(chem(c),'MIE')
        potVal{pa,1}(6,k)=chVal((100 - chem(c).MIE)/80);
    end

    if ~estVar(chem(c),'hdet50')
        potVal{pa,1}(6,k)=chVal(1-((1/log10(50)) * log10(chem(c).hdet50)));
    end
%HEALTH
%Irritation
pa=7;

    %DOT classification
    if ~estVar(chem(c),'LFL') && ~estVar(chem(c),'DOT')
        if ~isempty(strfind(chem(c).DOT,'rrit'))
            potVal{pa,1}(1,k)=(1);
        end
    end

    if chem(c).hflag==1
        potVal{pa,1}(2,k)=chem(c).hval(pa);
    end

    if ~estVar(chem(c), 'pH')
        potVal{pa,1}(3,k)=chVal((abs(chem(c).pH-7)-2)/5);
    end

    if ~estVar(chem(c),'Ldderm')
        if chem(c).Ldderm==0
            potVal{pa,1}(4,k)=chVal(0);
        else
            potVal{pa,1}(4,k)=chVal(1-.5*log10(chem(c).Ldderm/100));
        end
    end
    %EC classification
    %Test Output
    %pH
    %LD50 dermal


%Chronic Toxicity
pa=8;

    if ~estVar(chem(c),'TLV')
        potVal{pa,1}(1,k)=1-log10(chem(c).TLV)/3;
    end



    %TLA-TWA MSDS or ACGIH
    %REL-NIOSH

    %TEST output

    if chem(c).hflag==1
        potVal{pa,1}(3,k)=chem(c).hval(pa);
    end

%ENVIROMENT
%Water effects
pa=9;

    if ~estVar(chem(c),'FMLC50')
        potVal{pa,1}(1,k)=chVal(chem(c).FMLC50/5);
    end
    if ~estVar(chem(c),'MGLC50')
        potVal{pa,1}(2,k)=chVal(chem(c).MGLC50/5);
    end
    if ~estVar(chem(c),'IGC50')
        potVal{pa,1}(3,k)=chVal(chem(c).IGC50/5);
    end


    if chem(c).hflag
        potVal{pa,1}(4,k)=chem(c).hval(pa);
    end

    %LC50 mg/l
    %TEST Output
    %DOT dangerous to marine life symbol
    %R- Codes

%Air effects
pa=10;
    %ERPG-2
    if ~estVar(chem(c),'ERPG2')
        potVal{pa,1}(1,k)=1-log10(chem(c).ERPG2/5)/3;
    end
    %RfD
    if ~estVar(chem(c),'Rfd')
        potVal{pa,1}(2,k)=chVal(1-log10(chem(c).Rfd)/3);
    end
    %chronic tox index
    if chem(c).hflag
        potVal{pa,1}(3,k)=chem(c).hval(pa);
    end

    if length(potVal)>7
        if ~isempty(potVal{8}) && size(potVal{8},2)>=k
            potVal{pa ,1}(4,k)=max(potVal{8,1}(:,k));
        else
            potVal{pa ,1}(4,k)=chVal(0);
        end
    end
%Solid effects
pa=11;
    %Yes or no
    potVal{pa,1}(1,k)=chVal(0);

%Bioaccumalation
pa=12;

    if ~estVar(chem(c),'BCFreal')
        potVal{pa,1}(1,k)=chVal(0.5*(chem(c).BCFreal)-1);
    elseif ~estVar(chem(c),'BCFpred')

        potVal{pa,1}(2,k)=chVal(0.5*(chem(c).BCFpred)-1);
    end

    if ~estVar(chem(c),'Kow')
        potVal{pa,1}(3,k)=chVal(0.5*(chem(c).Kow)-1.5);
    end

    if chem(c).hflag
        potVal{pa,1}(4,k)=chem(c).hval(pa);
    end
    %log BCF TEST output
    %log Kow octanol / water partition coefficent
    %R-Codes

%Degredation
pa=13;

    if ~estVar(chem(c),'Percistency')
        potVal{pa,1}(1,k)=chVal(.5*log10(chem(c).Percistency));
    end

    %persistency Ecotox

    if ~estVar(chem(c),'Halflifereal')
        potVal{pa,1}(2,k)=.5*log10(chem(c).Halflifereal);
    end
    %Halflife -obs
    %halflife -literature

    if ~estVar(chem(c),'OCED')
        %learn what OCED is
%         potVal{pa,1}(3,k)=1-1.1*log10(chem(c).OCED/20);
    end
    %OCED biodegradability after 28 days
    %BOD/COD ratio of biological oxygen need to chemical
    %
    %
    if chem(c).hflag==1
        potVal{pa,1}(4,k)=chem(c).hval(pa);
    end

    %just as a place holder
    maxFlag(pa,1)=0;
    potVal{pa,1}=0;
    k=k+1;
end



%Combine inputs to a single index.
paraLen=13;
IntVal = zeros(paraLen, chmLen);
for pw = 1:paraLen
    diff=size(potVal{pw},2);
    if diff~=chmLen;

        for err=diff+1:chmLen
            IntValJust{pw,err}='no data';
        end

    end
    for dd=1:diff
        tr=find(potVal{pw}(:,dd),1);
        if isempty(tr)
            IntVal(pw,dd)=0;
            IntValJust{pw,dd}='no data';
        else
            IntVal(pw,dd)= potVal{pw}(tr,dd);
            if maxFlag(pa)
                IntValJust{pw,dd}=justify{2}{pw,justify{1}(pw,2)};
            else
                IntValJust{pw,dd}=justify{2}{pw,tr};
            end
%             temp=temp+IntVal(pw,dd);
        end
    end
    IntValRaw(pw,1)=mean(nonzeros(IntVal(pw,:)));
    IntValRaw(isnan(IntValRaw))=0;
%     IntVal(pw,1:diff)= max(potVal{pw}, [],1);Max
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fate = zeros(paraLen,chmLen);
k=1;
for c = chemNums
    fl=0;
    for r=1:length(chem(c).RG)
        if ~isempty(find(organics==chem(c).RG(r),1))
            fl=1;
        end
    end
    %Enviro effects are affected by degraddation and accumalation
    if fl==1
        fate(9,k)=IntVal(13,k)*IntVal(12,k)*IntVal(9,k);
        fate(10,k)=fate(9,k);
    else
        fate(9,k)=IntVal(12,k)*IntVal(9,k);
    end

    %fire, explosion and irritation are affected by mobility
    %toxcity is very effected by mobility.
    for pe = [2,3,6,7,8]
        if pe == 3 || pe == 8
            fate(pe,k) = 1.25 *IntVal(1,k)*IntVal(pe,k);
        else
            fate(pe,k) = IntVal(1,k)*IntVal(pe,k);
        end

    end

    %Reaction potential and reaction probability are affected if both exist
    fate(4,k)= IntVal(4,k)*IntVal(5,k);
    fate(5,k)= IntVal(4,k)*IntVal(5,k);

    k=k+1;
end

% Enforce maximum of 0.5
lo = fate>0.5;
fate(lo) = 0.5;
edp=real(IntVal+fate);

for pw = 1:paraLen
    fateRaw(pw,1)=mean(nonzeros(fate(pw,:)));
    fateRaw(isnan(fateRaw))=0;
    edpRaw(pw,1)=mean(nonzeros(edp(pw,:)));
    edpRaw(isnan(edpRaw))=0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MIX parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%could use accentric factor to lower bounds due to uncertainty in ideality

ydiff= u(n).ytheo-(1/chmLen);

pa=1;
    %if electrolytes weight using (z1+z2)/((z1/D)+(z2/D))
    %otherwise get the inverse molar masses (MM-1+MM2-1...)^.5
    %what if modeled as a puff? here finds an effective mobility
    temp=0;

    for ch=chemNums
        temp=temp+1/chem(ch).MW(1);

    end

    potmix(pa,1)=chVal(((temp/(chmLen/18))^.5)-1);

    %delf
    %Anything that requires filling (aerosol)- check var=flowrate in
    %Phase change or interaction between flows - checkvar= delH vap
    %Change in temperature or pressure (greater than 40 deg breaks heat)
    %Corrosion

%fire
pa=2;

    %Use componnent with lowest flash
    %1)get psat at flash point
    %2)psatmix=psatflash/x
    %3)Antoine (psatmix)->Tflash mixture
    %4)a*(Tflashmix/Tflash)+b

    %Fuel + oxidizer
    fuel=0;
    oxid=0;
    for ch=chemNums
        fuel=fuel+ u(n).ytheo(ch)*chem(ch).flagFuel;
        oxid=oxid+ u(n).ytheo(ch)*chem(ch).flagOx;
    end
    potmix(pa,1)=chVal(2*fuel*oxid);
    %Use Le chatliers principle to get LFLmix (combustible basis mole
    %fraction in air 1/(ycom/LFL)





    %effective combustion energy
%     [chem.dHburn]*ydiff
    u(n).LFLmixburn=(11200*100*4.184/-([chem(chemNums).dHburn]*yt));
    potmix(pa,2)=chVal((u(n).LFLmixburn-mean([chem(chemNums).LFL]))/5,999);

    if ~max(estVar(chem, 'LFL'))

        temp=0;
        for ch=chemNums;
            temp=temp+ycomb(ch)/chem(ch).LFL;
        end
        u(n).LFLmix=1/temp;
        potmix(pa,3)=chVal(1-(1/(15*temp)));
    end

%     mix(pa,1)=max(potmix(pa,:));
    %continue to get the index as before then reduce it so it can be used
    %as a penalty

    %average both values to get a final penalty

    %delf
    %
    %Corrosion

%Acute Toxicity
pa=3;

    %p87 Use a flash calculation with roults law to get a mixture value of

    la=estVar(chem, 'TLVc');
    la=[la(chemNums)];
    if max(la)
        for li=find(la)'
            chem(chemNums(li)).TLVc=10*(10^(4*(1-IntVal(pa,li))));
        end
    end

    %Exposure

    u(n).acExposure=(1./([chem(chemNums).TLVc]))*u(n).relExp{1};
%     u(n).acExposure=((10^6)/5000)*yWeight(chem, u ,n, 3 , 'TLVc', 'all');
    potmix(pa,1)=chVal(.5*(u(n).acExposure-1),999);

    %Le chatliers LFLmix
    temp=0;
    for ch=chemNums;
        temp=temp+(u(n).ytheo(ch)/chem(ch).TLVc);

    end
    u(n).TLVcmix=1/temp;
    potmix(pa,2)=chVal(-.25*log10(u(n).TLVcmix/(10)),999);
%         potmix(pa,2)=chVal(-(u(n).TLVcmix-mean([chem(chemNums).TLVc]))/50,999);
%     else
%         temp = 0;
%         kl=1;
%         for ch=chemNums;
%             temp=temp+(u(n).ytheo(ch)/IntVal(pa,kl));
%             kl=kl+1;
%         end
%         mix(pa,1)=1/temp;
%     end

    %delf
    %
    %Corrosion


%Prob Runaway Decomp
pa=4;

    %use boxes in CRW
    lo=tril(ones(chmLen),-1);
    exclude=[];
    for eq=1:length(chem)
        if isempty(simple(eq))
            exclude=[exclude simple(eq)];
        end
    end

    if ~isempty(exclude)
        for hu=exclude
            lo(hu,:)=zeros(1,chmLen);
            lo(:,hu)=zeros(chmLen,1);
        end

%         for re=1:length(lo)
%             for ce=1:length(lo)
%                 if mat(re,ce)==1
%                     if ~isempty(find(exclude,re)) || ~isempty(find(exclude,ce))
%                         lo(re,ce)=0;
%                     end
%                 end
%             end
%         end


    end

    numBox= length(find(lo==1));
    lo=logical(lo);
    u(n).CRWMat=chem(1).ReacMat(lo)+chem(1).GasMat(lo);
    potmix(pa,1)=chVal(((sum(u(n).CRWMat))/numBox)-0.2,999);

    %UU
    %Combine all the heat of combustion of the surrounding units to
    %simulate an external fire leading to a BLEVE

%Expected damage
pa=5;

%Explosion
pa=6;

    %mixture fundamental burning velocity
    la=estVar(chem, 'velsu');
    la=[la(chemNums)];
    if max(la)
        for li=find(la)'
%             chem(c).velsu/50)-.9);
            chem(chemNums(li)).velsu=(IntVal(pa,li)+.9)*50;
        end

    end

    temp=[];
    res=yWeight(chem,u,n,2,'MW','Fuel');
    nums=[find(res{1})]';
    if ~isempty(nums)
        fueloxmix=res{2}/(32*mean([chem(res{1}).oxyStoic]));
        k=1;
        for fl=nums
            eqv=fueloxmix/chem(fl).fueloxopt;
            temp(k)=chem(fl).velsu*...
                (Tk/400)^(1.783-.375*(eqv))*...
                (P)^(-.17*eqv^(.5*sign(eqv)));
            k=k+1;
        end

        u(n).velsumix=(temp*u(n).ytheo(res{1}));
        potmix(pa,1)=chVal((u(n).velsumix-mean([chem(res{1}).velsu]))/50,999);
    end




    %Use delGibbsmixture (gg)

    %delf: use mechanical explosions / confined space /missile damage
    %potential


    %delf
    %Turbulence of mixture CSTRs, Extractors etc
    %Pressure
    %Corrosion

    %BLEVE pot -if Tb<T gamma heat capacity ratio * vapor space
    %flash mixture

    %POIUHGSPIG NEED: Vvapor, totmoles, Ewithstand (how much energy can the
    %vessel withstand...Alternative just use the second pressure with hte
    %pressure rating)...or use P/Psec ratio as the indicator


%     tol=.05;
%     flag=0;
% result=0;
% while abs(1-result)>tol
%     if result==0
%
%     else
%         Vguess=((1-result)*Vguess)/10+Vguess;
%         if Vguess < 10^-10
%             Vguess =.9;
%         end
%         result=0;
%     end
%

%     for cc=chemNums
%
%         tek=u(n).psat(cc)/P;
%         result=result+(tek*u(n).x(1,cc))/(1+Vguess*(tek-1));
%     end
%     (1-result)
% end

%BLEVE
    if u(n).Phase(1)<1;
        if fv <1 && fv>0
            Vguess=fv;
        else
            Vguess=.3;
        end

        fun = @(x)(mixFlash(u(n),x));
        if isreal(mixFlash(u(n),Vguess))
            u(n).fvmix= fzero(fun,Vguess);
        else
            u(n).fvmix =-999;
        end
        if u(n).fvmix <1 && u(n).fvmix>0
            moladd=u(n).Amount(1)*u(n).fvmix;
            Psec=moladd*.0821*Tk/u(n).Vvapor;
            E=(Psec*u(n).Vvapor*log(Psec/P))/101.325; %Joules
        %     u(n).mix=log(E/u(n).Ewithstand);
            potmix(pa,2)=chVal(0.9*log(.05*E/300000)-1,999);%2*the range of exps with a range of -1 to 1
        else
            potmix(pa,2)=0;
        end
    end
    %UU
    %Combine all the heat of combustion of the surrounding units to
    %simulate an external fire leading to a BLEVE

%HEALTH
%Irritation
pa=7;

    %DOT classification
    %EC classification
    %Test Output
    %pH
    %LD50 dermal
    %Mix- combined pH? Combined effect?

%Chronic Toxicity
pa=8;

    %Mix the TLV to get a TLV for the whole mixture
    la=estVar(chem, 'TLV');
    la=[la(chemNums)];
    if max(la)
        for li=find(la)'
            chem(chemNums(li)).TLV=10^((1-IntVal(pa,li))*3);
        end
    end

    %Exposure
    u(n).chExposure=(1./([chem(chemNums).TLV]))*u(n).relExp{2};
%     u(n).chExposure=((10^6)/500)*yWeight(chem, u ,n, 3 , 'TLV', 'all');
    potmix(pa,1)=chVal(.5*(u(n).chExposure-1),999);

    %Le chatliers LFLmix
    temp=0;
    for ch=chemNums;
        temp=temp+(u(n).ytheo(ch)/chem(ch).TLV);

    end
    if temp>0 && ~isnan(temp)
        u(n).TLVmix=1/temp;
        potmix(pa,2)=chVal(1-(log10(u(n).TLVmix)/3),999);
    end


%         potmix(pa,1)=chVal(-(u(n).TLVmix-mean([chem(chemNums).TLV]))/50,999);
%     else
%         temp = 0;
%         kl=1;
%         for ch=chemNums;
%             temp=temp+(u(n).ytheo(ch)/IntVal(pa,kl));
%             kl=kl+1;
%         end
%         if temp>0 && ~isnan(temp)
%             mix(pa,1)=1/temp;
%         end
%     end

%ENVIROMENT
%Water effects
pa=9;

    %LC50 mg/l
    %TEST Output
    %DOT dangerous to marine life symbol
    %R- Codes
    %Mix- check EcoTox for datam

%Air effects
pa=10;

    %ERPG-2
    %chronic tox index



%Solid effects
pa=11;
    %Yes or no

%Bioaccumalation
pa=12;

    %log BCF TEST output
    %log Kow octanol / water partition coefficent
    %R-Codes

%Degredation
pa=13;

    %persistency Ecotox
    %Halflife -obs
    %halflife -literature
    %OCED biodegradability after 28 days
    %BOD/COD ratio of biological oxygen need to chemical
    %
    %
    potmix(pa,1)=0; %place holder
%
%     mix=max(potmix,[],2);
%     [junk,junk,lo]=find(potmix);

    for pi=1:pa
%         [temp in]=max(nonzeros(potmix(pi,:)));
        [junk in temp]=find(potmix(pi,:),1);


        if isempty(temp)
            mix(pi,1)=0;
            mixJust{pi,1}='No data';
        else
            mix(pi,1)=temp;
            mixJust{pi,1}=justify{3}{pi,in};
        end
    end

%%%%%%%%
%VESSEL
%%%%%%%%.

%Fail rate, f
a=1;% get these from user input showing how different their technologies are compared to the ones in the literature.
b=.1; %nominal value
c=1; %from the average fail rate of the industry
failBar= a*log10(u(n).failbase/c)+b;

%delf, change in failure rate of the system

%temperature
pf=1;

    if u(n).classFunc ==3
        %reactor so temperature is a big indicator of bad things
        mult=1.2;
    else
        mult=1;
    end

%     delf(pf,1)=chVal(mult*log10(max(25-T, T-100)/25),1);


%Pressure
pf=2;
    % use thermodynamic availability from Crowl "calc the engery of
        %its most conservative and allows for heat loss to the enviro
    MechExp = (P*u(n).V*(log(P)-(1-1/P)))/101.325; %Joules
    expansion(1)=.001304*P*101.325*u(n).V;
    expansion(2)=1/(Tk*1000) * (P-mean(u(n).ppar))^2 *u(n).V;
%     delf(pf,1)=chVal(0.9*log10(.05*MechExp/300000)-1,1); %2000 of mechExp is about 1 ton of TNT




%Alternative
if u(n).multIn
    in=u(n).inPsu;
else
    in=u(n).in;
end



% adjust for the process extremity
if u(n).classFunc~=3
    %mixers, seps, transport and reactors
    %Atmospheric vs refrig vs pressurized

    [u(n).fail u(n).df]= augfailTP(u(n).T,u(n).P,u(n).V,u(n).failadj,u(n).failadjTM,u(n).failadjPM) ;
else
    if u(n).classVar(1)==1
        %heater
        %temp range


        vin=sort([(1+log10(u(n).dTLM/40)),.8,.4]);
        u(n).fail=u(n).failadj*vin(2);
        u(n).df(1)=vin(2);
        u(n).df(2)=0;


    else
        %compressor
        %power req.

        u(n).df(1)=0;
        u(n).df(2)=0;
    end
end




%delfDanger, turning delf into dangerous properties
% lossCont = delf(1)+delf(2);

% lossCont = chVal((1/90)*(u(n).fail*10^6)-(1/9),999);
% lossCont=zeros(13,1);

u(n).lossCont = chVal((1/3)*log10(u(n).fail*10^5),999);

% delfDanger([2,3,5,6,12,13],1) = lossCont;


% delfDanger(1,1)=0;
% delfDanger(2,1)=lossCont;
% delfDanger(3,1)=lossCont;
% delfDanger(4,1)=0;
% delfDanger(5,1)=lossCont;
% delfDanger(6,1)=lossCont;
% delfDanger(7,1)=0;
% delfDanger(8,1)=0;
% delfDanger(9,1)=0;%cause sig figs
% delfDanger(10,1)=0;%cause sig figs
% delfDanger(11,1)=0;%cause sig figs
% delfDanger(12,1)=lossCont;
% delfDanger(13,1)=lossCont;

% Vessel= delfDanger*failBar;
conq=IntValRaw>0.05;

u(n).Vessel = u(n).lossCont.*u(n).detect.*(conq-u(n).armorConq);

%END CALC
%No more indexes to calculate


    u(n).IntVal=IntVal;
    u(n).IntValJust=IntValJust;
    u(n).potVal=potVal;
    u(n).fate= fate;
    u(n).edp=edp;
    u(n).mix=mix;
    u(n).mixJust=mixJust;
    u(n).potmix=potmix;
    u(n).IntValRaw=IntValRaw;
    u(n).fateRaw=fateRaw;
    u(n).edpRaw=edpRaw;

    clear('IntVal', 'IntValJust', 'potVal', 'fate', 'mix', 'mixJust', 'potmix', 'failBar',...
        'IntValRaw', 'fateRaw', 'edp', 'edpRaw', 'conq');



%
%     plant.u(num)=new;
%     plant.chem=chem;
end

function [val] = yWeight(chem, u ,n, flag , var, subset)
chmNums=u(n).compList;
lo=zeros(length(chem),1);

if strcmp(subset,'all')
    for ie=chmNums
        lo(ie)=1;
    end
else
    for ie=chmNums

        if chem(ie).(['flag' subset])
            lo(ie)=1;
        end
    end
%     lo=chmNums.*logical([chem(chmNums).(['flag' subset])]);
end
lo=logical(lo);
if flag==3
    exp=-1;
%     res=((1./[chem(lo).(var)])*[u(n).ytheo(lo)]);
    flag=1;
else
    exp=1;
%     res=([chem(lo).(var)]*[u(n).ytheo(lo)]);
end
temp=[];
k=1;
for li=[find(lo)]'
    temp(k)=chem(li).(var)^exp;
    k=k+1;
end


    res= temp*[u(n).ytheo(lo)];
%     temp=temp+((chem(li).(var))^exp)*u(n).ytheo(li);


if flag==2
    val{1}=lo;
    val{2}=res;
elseif flag
     val=res;
else
    val=lo;
end

end

function val = mixFlash (u,Vguess)
result =0;
for cc=u.compList

        tek=u.psat(cc)/u.P;
        result=result+(tek*u.x(1,cc))/(1+Vguess*(tek-1));
end
val = 1-result;
end
%TNTeq(chem,u,c,n,'dHburn')
function val =TNTeq(chem,u,c,n,name)
if strfind(name,'reac')
    in=u(n).(name);
else
    in= chem(c).(name);
end
if sign(in)>0
    val=0;
else
    val=-.01*in*...
        u(n).x(1,c)*u(n).Amount(1)*u(n).Activity(1,c);
end

end

function [per ] = ProbTNT(TNTamt,dist)

z=dist/TNTamt^(1/3);

mag=@(x,c) (1+(x/c)^2)^(.5);

overp = (1616*mag(z,4.5)^2)/(mag(z,.048)*mag(z,.32)*mag(z,1.35));
% overp*101.325
% probt= -23.8+2.92*log(overp*101325)
probt= -17.79+2.18*log(overp*101325);
per=50*(1+((probt-5)/abs(probt-5))*erf(abs(probt-5)/2^(.5)));
end

function val =chVal (varargin)
if length(varargin)>1
    flag=1;
else
    flag=0;
end

in=real(varargin{1});
if isnan(in)
    in=0;
end

if flag
    temp = sort([-1,in,1]);
    val= temp(2);
else
    temp = sort([10^-10,in,1]);
    val= temp(2);

end
end



%
%
%  function [val]=antoine(input,constants,type)
%  %antoine Get saturated pres. or temp. according to NIST
%  %  log10(P) = A ? (B / (T + C))
% %   val is the output in bar or T
% %   T is in Kelvin
% %   P is in bar
% %
% % A=[13.7819,13.9320,14.0579,4.19927];
% % B=[2726.81,3056.96,3331.45,1569.622];
% % C=[27.572,217.625,214.627,-63.572];
%
% A=constants(1);
% B=constants(2);
% C=constants(3);
%
% %figure out if user wants P or T out
% if ~isempty(strfind(type,'T')) || ~isempty(strfind(type,'t'))
%     flag=0;
%     if ~isempty(strfind(type,'C'))
%         T=input+273.15;
%     else
%         T=input;
%     end
% else
%     P=input;
%     flag=1;
% end
%
% %use antoine eq to get the unknown variable
% if flag==0
%     val=10^(A-(B/(T+C)));
% else
%     val= (B/(A-log10(P)))-C-273.15;
% end
% %  end

%
% function [val]= chkLowLimit (val)
% val=max(val);
% if val <0
%     val=0;
% end
%
% end
