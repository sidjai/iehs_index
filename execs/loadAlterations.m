function [val] = loadAlterations(numAlt,des,par)
%loadAlterations Intialize GUI and capture the responses.

%[left bottom width height]
%des={'first','second','third','fourth','fifth','sixth'};
%numAlt=1;
but={};


vdist=1/11;
widEntry=.15;
% hdist=[.35/3, 2*.35/3, .1];
hdist{1} = [0 ,((1-widEntry*numAlt)/3)];
hdist{2} = [hdist{1}(2) ,2*((1-widEntry*numAlt)/3)];
for i=1:numAlt
    hdist{2+i} = [sum(hdist{1+i}) ,(widEntry)];
end

add= min(0.95,.4+widEntry*numAlt);
han = figure('Units','normalized','Position',[.05,.05,add,.85]);
B = get(han, 'color');
bot=1-vdist;
makeText('Unit Alterations',{[0 1]},vdist,bot,1,B,'title');
% uicontrol(han, 'Style', 'text', 'String','\fontsize{16}(Unit Alternatives)','Units','normalized',...
%     'Position', [0 bot 1 vdist],'foregroundcolor', [0 0 0], 'backgroundcolor', BackColor);

bot= bot-2*vdist;
% text(.2,1-2*vdist,des,...
%     'Units','normalized','HorizontalAlignment','center');
makeText(des,{[.05 .35]},2*vdist,bot,1,B,'big');
% uicontrol(han, 'Style', 'text', 'String',des,'Units','normalized',...
%     'Position', [.05 bot .35 3*vdist],'foregroundcolor', [0 0 0], 'backgroundcolor', BackColor);


bot= bot-vdist;
%Where the process effects the system
%   Max or Base

makeText('Where does the addition apply?',hdist,vdist,bot,1,0,'big');
s={'Max: protection against extreme conditions',...
   'Base: overall reliability'...
   'Detection: Ability for an event to be anticipated or alerted.',...
   'Consequences: mitigates the hazard after an event'...
   };
makeText(s,hdist,vdist,bot,2,0,'reg');

for a=1:numAlt
    but{1,a} = uicontrol(han, 'Style', 'listbox','max',4,...
    'String', {'Max','Base', 'Detection','Consequences'},...
    'Units','normalized','Position', [hdist{2+a}(1) bot hdist{2+a}(2) vdist]);
end


bot= bot-2*vdist;
%Type of mechanism
%   Passive or active->computerized or otherwise

makeText('Type of mechanism',hdist,2*vdist,bot,1,B,'big');
makeText('[image of flowchart goes here]',hdist,2*vdist,bot,2,B,'reg');
for a=1:numAlt
    but{2,a} = uicontrol(han, 'Style', 'edit','Units','normalized',...
    'Position', [hdist{2+a}(1) bot hdist{2+a}(2) vdist],...
    'foregroundcolor', [0 0 0], 'backgroundcolor', B);
end

bot= bot-vdist;
%Protection against what condition?
%   Temperature defense, Pressure defense or flow defense

makeText('Protection against what condition?',hdist,vdist,bot,1,0,'big');
makeText('Defends against excessive temperature...etc.',hdist,vdist,bot,2,0,'reg');
for a=1:numAlt
    but{3,a} = uicontrol(han, 'Style', 'listbox','max',2,...
    'String', {'T','P'},'Units','normalized',...
    'Position', [hdist{2+a}(1) bot hdist{2+a}(2) vdist]);
end


bot= bot-vdist;
%Strength of the system
%   How much it actually protects against the condition

makeText('Strength',hdist,vdist,bot,1,B,'big');
makeText({'Qualatative assesment of the protection the system offers',...
    '0 for only temporary relief or a notification',...
    '5 for state of the art control'},...
    hdist,vdist,bot,2,B,'reg');
for a=1:numAlt
    but{4,a} = uicontrol(han, 'Style', 'edit','Units','normalized',...
    'Position', [hdist{2+a}(1) bot hdist{2+a}(2) vdist],...
    'foregroundcolor', [0 0 0], 'backgroundcolor', B);
end

bot= bot-vdist;
%Relaibility of the system
%   confidence of the ability to withstand wear and tear.

makeText('Reliability',hdist,vdist,bot,1,0,'big');
makeText({'Best guess for the ability of the system to withstand wear and tear',...
    '0: MTTF ~ 1 month',...
    '5: for failing within 5 years.'},...
    hdist,vdist,bot,2,0,'reg');
for a=1:numAlt
    but{5,a} = uicontrol(han, 'Style', 'edit','Units','normalized',...
    'Position', [hdist{2+a}(1) bot hdist{2+a}(2) vdist]);
end
bot= bot-vdist;
%Maintainability
%   Ease of access and repairability

makeText('Maintainability ',hdist,vdist,bot,1,B,'big');
makeText({'Qualatative assesment of the ease of access and repairability in terms of shut down hours',...
    '0 shutdown time on the order of months with esoteric parts and crampoed quarters',...
    '5 little or no shutdown time with hardly any skill required to access or repair the system'},...
    hdist,vdist,bot,2,B,'reg');
for a=1:numAlt
    but{6,a} = uicontrol(han, 'Style', 'edit','Units','normalized',...
    'Position', [hdist{2+a}(1) bot hdist{2+a}(2) vdist],...
    'foregroundcolor', [0 0 0], 'backgroundcolor', B);
end

set(but{1},'Callback',{@firstbut,but,par});

%make the clear push button

bot= bot-.75*vdist;

but{7} = uicontrol(han,'Style','pushbutton','String','Done',...
    'Units','normalized','Position', [.5*(1- bot*hdist{3}(2)) bot hdist{3}(2) .5*vdist],...
    'Callback',{@donebut,but});
waitfor(but{7},'UserData')
val=get(but{7},'UserData');
close(han)
end

function makeText(string,hdist,vdist,bot, col,back,bigFlag)
t = uicontrol('Style', 'text', 'String',string,'Units','normalized',...
    'Position', [hdist{col}(1) bot hdist{col}(2) vdist]);
if back~=0
    set(t,'foregroundcolor', [0 0 0], 'backgroundcolor', back);
end

if strcmp(bigFlag,'title')
    set(t,'FontSize',20,'FontWeight','demi')
elseif strcmp(bigFlag,'big')
    set(t,'FontSize',14)
elseif strcmp(bigFlag,'big')
    set(t,'FontSize',13)
end

end
function firstbut(obj,event,but,par)
temp={};
if ~isempty(find(get(obj,'Value')<3,1))
    temp={'T','P'};
end
if ~isempty(find(get(obj,'Value')>2,1))
    te=1;
    for k=length(temp)+1:length(temp)+length(par)
        temp{k}= par{te};
        te=te+1;
    end
end
set(but{3},'String',temp)

end
function donebut(obj,event,but)
num=size(but,2);
for b=1:length(but)
    for a =1:num
        
        if strcmp(get(but{b,a},'Style'),'edit')
            val{b,a}=str2double(get(but{b,a},'String'));
        else
            val{b,a}=get(but{b,a},'Value');
        end
    end
        
    
    
    
end
set(obj,'UserData',val)
end