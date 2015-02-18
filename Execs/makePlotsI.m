function [ u , s , val ] = makePlotsI( u, s, flag )
%makePlotsI  Make the various plots to show the where the risk is
%   Right now just the chemical dependent methods versus the design
%   dependent methods

untLen = size(u,1);
paraLen=size(u(1).IntVal,1);
%find the eligible streams
parameters ={'Mobility','Fire','Acute Toxicity','Prob Runaway Decomp','Expected damage','Explosion','Irritation','Chronic Toxicity','Water effects','Air effects','Solid Effects','Bioaccumalation','Degredation'};

eligible = [];
for su=1:length(s)
    if s(su).to~=999
        eligible=[eligible su];
    end
end

% Calculate the Unit's stress on other streams

for ue=1:untLen
    temp=zeros(paraLen,1);
    for oni=u(ue).out
        if ~isempty(find(eligible == oni,1))
            temp=temp + s(oni).UUs;
        end
    end
    u(ue).UUstressout = temp;
    
    temp=zeros(paraLen,1);
    for ono=u(ue).in
        if ~isempty(find(eligible == ono,1))
            temp= temp + s(ono).UUs;
        end
    end
    u(ue).UUstressin = temp;
    u(ue).UUstress =temp-u(ue).UUstressout;
    
end


ind=2;
xs=zeros(2*untLen+1,paraLen);
ys=xs;
xms=zeros(2*untLen+1,1);
yms=xms;
reMean = @(x)(mean(nonzeros(x)));
for ur=1:untLen
    if ur==1
        ys(ind,:) = (u(1).UUstressin)';
        yms(ind,1) = reMean(ys(ind,:));
        ind=ind+1;
    end
        
    
    xs(ind,:) = xs(ind-1,:)+(u(ur).IntValRaw)';
    xms(ind,1) = reMean(xs(ind,:));
    ys(ind,:) = ys(ind-1,:)+(u(ur).fateRaw)'+(u(ur).mix)'+(u(ur).Vessel)';
    yms(ind,1) = reMean(ys(ind,:));
    ind=ind+1;
    
    if ur~=untLen
        ys(ind,:)=ys(ind-1,:)+(u(ur).UUstressout)';
        xs(ind,:)=xs(ind-1,:);
        xms(ind,1) = reMean(xs(ind,:));
        yms(ind,1) = reMean(ys(ind,:));
        ind=ind+1;
        ad=3;
    else
        ad=2;
    end
    
    % do the labels
    
    
    for pe=1:paraLen
        labx(ur,pe) = xs(ind-ad,pe)+(xs(ind-ad+1,pe)-xs(ind-ad,pe))/2;
        
    end
        labxm(ur,1) = xms(ind-ad,1)+(xms(ind-ad+1,1)-xms(ind-ad,1))/2;
        labym(ur,1) = max(min([yms(ind-ad-1:ind-ad+1)])-.7,0);
%         if ind<untLen-1
%             labym(ur,1) = max(max([yms(ind-ad-1:ind-ad+3)])+5,0);
%         else
%             labym(ur,1) = max(max([yms(ind-ad-1:ind-ad)])+.4,0);
%         end
%         labym(ur,1) = max(max([yms(ind-ad-1:ind-ad)])+.6,0);
end

xs(isnan(xs))=0;
ys(isnan(ys))=0;




if flag
    
    
    plot(xms,yms,'MarkerSize', 5,'LineWidth',3,'Color',[80/255 0 0])
        xlim([-.1 (xms(end)+.3)])
        ylim([-.1 (xms(end)+.3)])
        set(gca,'FontWeight','demi');
        set(gca,'FontName','Arial');
        for ul=1:5:untLen
            text(labxm(ul),labym(ul),['\fontsize{14}' u(ul).name]);
        end
            
        xlabel('Chemical dependent values','FontName','Calibri','FontWeight','demi', 'FontSize', 20)
        ylabel('Design dependent values','FontName','Calibri','FontWeight','demi', 'FontSize', 20)
%         title([relArea{1,tt+1} ' Calibration curve'],'FontWeight','bold', 'FontSize', 20)
    figure(2) 
    set(gca,'FontWeight','demi');
    hold on;
    type={'','-','--','-.',':'};
    for pa=2:5
        plot(xs(:,pa),ys(:,pa),[type{pa} 'x'],'MarkerSize', 5,'LineWidth',1.5)
    end
    xlim([-.1 2])
    box('on')
    xlabel('Chemical dependent values','FontWeight','demi', 'FontSize', 20)
    ylabel('Design dependent values','FontWeight','demi', 'FontSize', 20)

    h=legend(parameters{2},parameters{3},parameters{4},parameters{5},'Location', 'Best');
    set(h, 'Box', 'off');
    set(h, 'Color', 'none');
%     for pa=1:paraLen;
        
       
    
    
end
for p =1:paraLen
    val(p,1)=ys(end,p)/xs(end,p);
end
val(p+1)=yms(end)/xms(end);
% val(p+2)=std(val(1:p));
% ys
% xs
% vals = ys(end,:)./xs(end,:);
% val = [vals, ];

end

