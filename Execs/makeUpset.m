function [ newB, val ] = makeUpset( u,oldS, var, g,time)
%makeUpset from plant data simulates a 15 min upset in the given stream
%with the given variable under the guideword g.
%   Uses mass balances and heat balances for mixing of the old unit and the
%   amount added with the upset in a steady state calculation. The dynamic
%   case of addition as semi-batch is beyond teh scope especially with this
%   time frame
%   Things to do:
%   Use the upstream's class to determine which variables to deviate

findLine = @(c , string) (find(~cellfun('isempty', strfind(c,string))));

timehr = time/60;

% oldS=s(stream);

oldB=u(oldS.to);

    
devVars ={'T','P','F'};

ind= findLine(devVars,var);
dev=zeros(length(devVars),1);
% g
% dev=[0;1;1]
if oldS.from ~= 0
    
    val=mean(nonZo(u(oldS.from).Vessel));
    
   
    %dev(ind)= (g*.5*(mean(u(oldS.from).Vessel)+1));
    
    vex = sort([-.9,.9,(g*val)]);
    dev(ind)=vex(2);
   
else
%     g
    dev(ind)= g*0.05; %feed fluctuation
end

val=dev(ind);
% g={1;-1};

% dev
addMole = (oldS.F(1))*dev(3)*timehr; %kmol
addMass = (oldS.F(2))*dev(3)*timehr; %kg
addVol = (oldS.F(3))*dev(3)*timehr*60; %L

%tsMole = (oldS.F(1))*(1+dev(3))*timehr; %kmol
%tsMass = (oldS.F(2))*(1+dev(3))*timehr; %kg
%tsVol = (oldS.F(3))*(1+dev(3))*timehr*60; %L
%  = oldS.F{1}*dev(1)*timehr;
addStuff{1}=[addMole; oldB.Amount(1)];
addStuff{2}=[addMass; oldB.Amount(2)];
addXmole = oldS.x(1,:);
addXmass = oldS.x(2,:); % later add the composition deviation var
comp{1} = [addXmole; oldS.x(1,:)];
comp{2} = [addXmass; oldS.x(2,:)];
Q=(oldS.CP(1) *(oldS.T(1)*dev(1))); %J/kmol
Qi=Q*(addMole + oldS.F(1)*timehr); %J
if ~isempty(strfind(oldB.unitH{1},'cal'))
    Q=Q/4184; %j-cals then kg- gm
end

Qs=[Q*oldS.MW, Q, Q* (oldS.F(2)/3.6) ];

delHi = (Qi/4.184)+ addMass*oldS.H(2)*1000; %cal
% addEnth(1) = oldS.H(1)*addMass*1000 + Qi;
% addEnth(2) = oldS.H(2)*addMass*1000 + Qi;
% addEnth(3) = oldS.H(3)*addMass*1000 + Qi;
%start mixing things to get a new unit with the upset. 
newB=oldB;
    newB.IntVal=[];
    newB.potVal={};
    newB.fate= [];
    newB.mix=[];
    newB.Vessel=[];


newB.T = ((delHi*4.184)/(sum(addStuff{2})*oldB.CP(1)))+ oldB.T(1); %assumes constant CP interms of temperature and the change in the mixture.


newB.Amount = [sum(addStuff{1}) ; sum(addStuff{2})];

newB.V=oldB.V+addVol;
newB.Vvapor = oldB.Vvapor-addVol;
newB.P=oldB.P*(dev(2)+1);

newB.H(1) = (delHi+oldB.H(1)*oldB.Amount(1)*1000)/(newB.Amount(1)*1000);
newB.H(2) = (delHi+oldB.H(1)*oldB.Amount(1)*1000)/(newB.Amount(2)*1000);
newB.H(3) = newB.H(2)*(newB.F(2)/3.6);

for type=1:2
    fl = (addStuff{type}'*comp{type});
    newB.x(type,:) = fl./newB.Amount(type);
%     newB.x(2,:) = fl(2,:)./newB.Amount(2);
end

%Assign special properties

% newB.fail = augfailTP(newB.T,newB.P,newB.failadj,1*10^-6,.1*10^-6,newB.classFunc) ;
        

switch oldB.classFunc
    case 3
        if oldB.classVar(1)==1
            %Heater
            %change fail rate if it is outside the 40 degree difference
            
        else
            %Compressor
            
            
        end
    case 4
    newB.dHreac=oldB.dHreac+ newB.CP*(newB.T-25);
    
%     newB.level=oldB.level+ (addMass /(s.density(2)*1000))/oldB.size;
end

%change the chemical list to adjust for the new additions.
newB.compList = find(newB.x(1,:) > 10^-8);
newB.numComp = length(newB.compList);
tempList=[];
for j=1:size(newB.x,2)
    if newB.x(1,j)>10^-10
        tempList=[tempList,j];
    end
end

newB.compList=tempList;
newB.numComp=length(tempList);


end

