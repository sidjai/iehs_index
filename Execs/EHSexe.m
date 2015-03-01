clear all; close all;tic
%uses the output of the aspen program in terms of the stream list and
    %vessels
    %Sheet 1=Vessels, Sheet 2= streamlist
%Input the process safety information using the worksheet 
%Use WkstTemplete as a blank wkst to input data

    
%Input File names to be searched through
% InfoWkstFile = 'safetyInfoProofAspenwhy.xlsx';
% AspenStreams = 'ProofStreamsv2.4.xlsx';
% AspenPredictReport = 'proofofconcept  v2.4chem.txt';
% AspenBlocksName = 'proofofconcept  v2.4blocks.txt';
% CRWreport = 'proofofconcept crw.txt';
% outName = 'ProofEHS3_4.xlsx';

configFile = 'F:\EHSI\Aspen reports\DME\IEHSconfig.txt';
cfg = grabConfig(configFile);

%Global
cfg.parameters ={'Mobility','Fire','Acute Toxicity','Prob Runaway Decomp','Expected damage','Explosion','Irritation','Chronic Toxicity','Water effects','Air effects','Solid Effects','Bioaccumalation','Degredation'};

%load all the required files
[u , s ] = loadDesign(cfg); 
disp('Loaded Design')
[chem ] = loadChemInfo(cfg.AspenPredictReport, cfg.TESToutput);
disp('Loaded Chem files')


%make files if you don't have it.
[msg, go,u] = makeUserInput(u,chem,cfg);

if ~go
    msg{1}
    msg{2}
else
%only do the rest of the code if the safety wkst is filled and the Unit Alt is filled
%Load up variables

untLen = size(u,1);
stmLen = length(s);
chmLen = length(chem);

%Global variables

% justify{1} = {{'Gas Dispersion','Flashing liquid','Vaporizing liquid'}...Mobility
%     {'Flash Temp','LFL','LOC','SDS'}...Fire
%     {'TLC-C','IDLH','LC50 inhalation'}...Acute Tox
%     {'dTad','SDS','AIT','Self-reactivity, CRW'}...Prob Runaway
%     {'Reaction Energy','Energy of Decomposition','Combustion energy','Gibbs formation energy','SDS'}...Potential damage
%     {'Fudamental burning velocity','Energy of explosion','max Overpressure','Oxygen content','SDS'}...Explosion
%     {'DOT label','SDS','pH','LD50 dermal'}...Irritation
%     {'TLV-TWA','SDS'}...Chronic Tox
%     {'Fish toxicity','Small incects toxicity','Lower Growth of bacteria','SDS'}...Water effects
%     {'ERPG2','RfD','SDS','Chronic toxicity index'}...Air effects
%     {'Yes or no'}...Solid effects
%     {'BCF experimental','BCF prediction','Kow','SDS'}...Bioaccumalation
%     {'Percistency','Half-life','OCED','SDS'}}; %Degredation

% justify{2} = {{'MW'}...Mobility
%     {'Fuel+Oxi','LFL Burn mixed','LFL Chatlier'}...Fire
%     {'Exposure','Chatlier'}...Acute Tox
%     {'CRW hazards'}...Prob Runaway
%     {''}...Potential damage
%     {'Su mixed, Equivalence', 'BLEVE'}...Explosion
%     {''}...Irritation
%     {'Exposure','Chatlier'}...Chronic Tox
%     {''}...Water effects
%     {''}...Air effects
%     {''}...Solid effects
%     {''}...Bioaccumalation
%     {''}}; %Degredation

%vars intitalized


[chem ] = loadSafetyInfo(cfg.InfoWkstFile,cfg.CRWreport, chem);
disp('Loaded Safety Info')

%do the calculations that require some sort of interation over the index
[chem, u, s] = loadAddProp(chem,u,s,cfg.parameters);


for n=1:untLen
    
    disp(['Unit specfic indices:' num2str(n) '/' num2str(untLen)])
    
    %numCon is number of changing variables 
    %need to do it for all those different units

    for ch=1:(u(n).numCon)
        [u,chem] = calcIntVal(u,chem,n,ch,cfg.parameters,cfg.IOnotes);
    end
    
   
end

%make the HotSpot method
u =calcHotSpot(u,cfg.parameters);
for n=1:untLen
    [out, sum] = aggregate('unit',u,n,cfg.parameters);
    u(n).edp = out{1};
    u(n).mdp = out{2};
    u(n).cdp = out{3};
    u(n).udp = out{4};
    u(n).chemSmmy = sum;
end
%make the Unit unit interactions
paraLen=size(u(1).IntVal,1);

%find the streams that go to a vessel and can be used for UU int.
eligible = [];
for su=1:length(s)
    if s(su).to~=999 && ~s(su).flagPsu
        eligible=[eligible su];
    end
end

%Do the UUs for all the streams
k=1;
for sind=eligible
   disp(['Unit-Unit interactions: Stream ' num2str(k) '/' num2str(length(eligible))])
   s = calcUUvals(u, chem, s, sind,cfg.parameters,cfg.IOnotes);
   k=k+1;
end

%add vessels together
[out, sum] = aggregate('final',u,s,cfg.parameters);

fin.idp = out;
fin.unitsum = sum;
fin.avgUUs = mean([s.UUavgdev],2);

disp('wooooooooooooooooooooooo')

[u,s,fin.design] = makePlotsI(u,s,cfg.plotFlag);

%Output results
if cfg.xsFlag
    [sheet names] = makeXlsSheets(u,chem,s,fin,cfg);
    
    for fe= 1: length(sheet)
        if fe<4
            xlswrite(cfg.outName,sheet{fe},fe)
        else
            xlswrite(cfg.outName,sheet{fe},names{fe})
        end
    end
end
end
toc
%made by Siddarta Jairam 2014