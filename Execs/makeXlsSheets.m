function [ sheet names ] = makeXlsSheets( u,chem,s, r, cfg )
%makeXlsSheets Transfers the data into easy to read excel files about the results.
%   Each tab is a new method output with the final two being the unit and
%   stream list along with the appendix.

idp=r.idp;
design=r.design;


hierarchy = {'IntVal';'edp'; 'mdp';'cdp';'udp';'idp'};
names=hierarchy;
names{7}='Stream Table';
names{8}='Unit Table';
names{9}='Appendix';
add = {'';'SecEffs';'mix';'Vessel'; 'HotSpot';'UUs'};

easy = [4,1,5];%ind of lowest mean and highest in summary

untLen = size(u,1);
parLen = length(cfg.parameters);
chmLen = length(chem);

% %Load up the parameters

%Get the chem names and put them here
for c=1:chmLen
    headings{1}{c}=chem(c).ID; 
    headings{2}{c}=chem(c).ID; 
end

headings{3} = {'Avg Agg','Best case','Worst case'};



%do IntVal sheet

%do EDP sheet

for hi=1:6
    
    %go through the units and output it in the cell
    uind=1;
    col=1;
    cadd=0;
    
    
    
    if hi==6;
        nU=1;
    else
        nU=untLen;
    end
    
    while uind<nU+1
        did=0;
        ete=0;
        if hi~=6
            
            %do the labels to the side of the vals
            sheet{hi}{1,col}=['Unit ' mat2str(uind)];
            sheet{hi}{4,col}=u(uind).name;
            sheet{hi}{5,col}=u(uind).type;
            sheet{hi}{6,col}=u(uind).Indes;
            sheet{hi}{7,col}=u(uind).Outdes;
        else
            sheet{hi}{1,col}='All Units';
        end
        
        
        sheet{hi}{2,col}=hierarchy{hi};
        
        if hi<3
            typecols=1:size(u(uind).IntVal,2)+1;
            trans= u(uind).compList;
        elseif hi~=6
            typecols=1:4;
        else
            typecols=1:length(idp)+1;
        end
        
        col=col+1;
        for type=typecols
            tind=type-1;
%             if type>1
%                 tind=1;
%                 if hi==6
%                     ind = idpsum{3,easy(type-1)};
%                 else
% %                     ind = u(uind).(['sum' hierarchy{hi}]){3,easy(type-1)};
%                     
%                 end
%             else
%                 tind=type-1;
                
%             end
            
            
            for row =1:(parLen+1)
                if type==1
                    if row >1
                        sheet{hi}{row,col} = cfg.parameters{row-1};
                    end
                else       
                    if row ==1
                        if hi<3
                            sheet{hi}{row,col} = headings{hi}{trans(tind)};
                        elseif hi~=6
                            sheet{hi}{row,col} = headings{3}{type-1};
                        else
                            sheet{hi}{row,col} = r.unitsum{tind}{3,1};
                        end
                    else
                        if hi<3
                            
                            sheet{hi}{row,col} = rund(u(uind).(hierarchy{hi})(row-1,tind));
                        elseif hi==6
                            sheet{hi}{row,col} = rund(idp{tind}(row-1));
                        else
                            sheet{hi}{row,col} = rund(u(uind).(hierarchy{hi})(row-1,tind));
                            
                        end
                    end
                        
                end
%                 sheet{hi}{row,col}
                
            end
            %put the avg value at the bottom if its needed
            endrow=row+3;
            if hi>2
            if type == 1 
                sheet{hi}{row+2,col}='Average';
                sheet{hi}{row+3,col}='std';
%                 sheet{hi}{row+4,col}='ind';
%                 sheet{hi}{row+5,col}='Agregation';
%                 endrow=row+3;
%                 if hi>3
%                     sheet{hi}{row+6,col}='Previous';
%                     endrow=row+6;
%                 end
                
            else
                if hi<6
                    sheet{hi}{row+2,col} = rund(mean(u(uind).(hierarchy{hi})(:,1)));
                    sheet{hi}{row+3,col} = rund(std(u(uind).(hierarchy{hi})(:,1)));
                else
                    sheet{hi}{row+2,col} = rund(mean(idp{1}(:,1)));
                    sheet{hi}{row+3,col} = rund(std(idp{1}(:,1)));
                end
                    
%                 if hi~=7
                    
%                     sheet{hi}{row+2,col}=rund(u(uind).(['sum' hierarchy{hi}]){2,easy(type-1)});
%                     sheet{hi}{row+3,col}=rund(u(uind).(['sum' hierarchy{hi}]){2,2});
%                     sheet{hi}{row+4,col}=ind;
%                     unit=u(uind).(['unit' hierarchy{hi}]);
%                 else
%                     sheet{hi}{row+2,col}=rund(idpsum{2,easy(type-1)});
%                     sheet{hi}{row+3,col}=rund(idpsum{2,2});
%                     sheet{hi}{row+4,col}=ind;
%                     unit=idpunit;
%                 end

%                 for er=1:size(unit,2)
%                     sheet{hi}{row+4+er,col}=unit{ind,er};
%                 end
                
%                 endrow=row+4+er;
            end
            else
                endrow=row;
                
            end
            
            %%%%%%%%%%%%%%
            %do justification
            %%%%%%%%%%%%%%
            
            nurow=1;
            
            bege=endrow+4;
            
            
            
            sheet{hi}{bege,1}='Justification';
            if hi==1
                sheet{hi}{bege+1,1}='Descriptor vals';
            else
                sheet{hi}{bege+1,1}=add{hi};
            end
            if did==0
                for row =(bege+1):(parLen+bege)

                    if type==1
                        sheet{hi}{row,col+cadd}=cfg.parameters{nurow};
                    else
                        if hi==1
                            sheet{hi}{row,col}=u(uind).IntValJust{nurow,tind};
                        elseif hi==2
                            sheet{hi}{row,col}=rund(u(uind).fate(nurow,tind));
%                                 chmLen=size(u(uind).fate,2);
%                                 for ete=1:chmLen
%                                     sheet{hi}{row,col+cadd+ete-1}=rund(u(uind).fate(nurow,ete));
%                                 end
%                                 ete=ete-1;

                        elseif hi==6
                            sheet{hi}{row,col}=rund(r.avgUUs(nurow,1));
                            did=1;
                        else
                            sheet{hi}{row,col}=rund(u(uind).(add{hi})(nurow,1));
                            did=1;
                            if hi==3
                                sheet{hi}{row,col+1}=u(uind).mixJust{nurow};
                            end
                        end
                        
                        

                    end
                    nurow=nurow+1;
                end
                if did==1 && hi==4
                    %Do the fail outputs
                    row=row+1;
                    sheet{4}{row+1,col}='Standard Conditions';
                    sheet{4}{row+1,col+1}='T Extreme Conditions';
                    sheet{4}{row+1,col+2}='P Extreme Conditions';
                    sheet{4}{row+2,col-1}='Base fail rate';
                    sheet{4}{row+3,col-1}='Model Alterations fail rate';
                    sheet{4}{row+4,col-1}='Final fail rate (process + model)';
                    sheet{4}{row+6,col-1}='dfT, process';
                    sheet{4}{row+7,col-1}='dfP, process';
                    
                    sheet{4}{row+2,col}=u(uind).failbase;
                    sheet{4}{row+2,col+1}=u(uind).failbaseTM;
                    sheet{4}{row+2,col+2}=u(uind).failbasePM;
                    
                    sheet{4}{row+3,col}=u(uind).failadj;
                    sheet{4}{row+3,col+1}=u(uind).failadjTM;
                    sheet{4}{row+3,col+2}=u(uind).failadjPM;
                    
                    sheet{4}{row+4,col}=u(uind).fail;
                    
                    sheet{4}{row+6,col}=u(uind).df(1);
                    sheet{4}{row+7,col}=u(uind).df(2);
                    



                end
            end
            
            cadd=cadd+ete;
                    
                    
                    
            
            col=col+1;
        end
        col=col+1;%skip column between units
        uind=uind+1;
    
    end
    
    if  hi==6
        sheet{6}{1,col}='Design Factor';
        sheet{6}{5,col+3}='Design = (Sec effs + mix + Vessel+ UU)/IntVal';
        for d=2:(length(design))
            temp=design(d-1);
            if isinf(temp)
                temp='Infinity';
            else
                temp=rund(temp); 
            end
            sheet{6}{d,col}=temp;
        end
        
        sheet{6}{d+2,col}=rund(design(end));
        
        
    end
end

sheet{7} = makeStreamTable(s, cfg.parameters,  {u(:,1).name});
sheet{8} = makeUnitTable( u, cfg.parameters  );
    
%Appendix

sheet{9}{1,1} = 'APPENDIX';

sheet{9}{3,1} = 'Version';
sheet{9}{3,2} = '0.0.1';
sheet{9}{4,1} = 'Run Time';
sheet{9}{4,2} = num2str(clock);
sheet{9}{6,1} = 'Initial parameters';
initial = fieldnames(cfg);

for it = 1:length(initial)
    fd = initial{it};
    if(~iscell(cfg.(fd)))
        sheet{9}{7+it,1} = fd;
        sheet{9}{7+it,3} = cfg.(fd);
    end
    
end




% sheet{6}=makeStreamTable( s );
%round the values off
% for cg=[6,7,10]
%     for ro=1:size(sheet{6},1)
%         if ~ischar(sheet{6}{ro,cg})&& ~isempty(sheet{6}{ro,cg})
%             sheet{6}{ro,cg}=rund(sheet{6}{ro,cg});
%         end
%     end
% end
     
end

function out = rund(in)
messy=floor(in*1000)/1000;
if messy==0
    if in==0
        out=0;
    else
        out=10^-10;
    end
else
    str=num2str(in);
    ow=textscan(str,'%5.3f');
    out=ow{1}(1);
end
% if out==0 && in~=0
%     out=10^-10;
% end
end

