
for ba=1:6
    [~,p(ba,:),~,t{ba}]=ttest(squeeze(data.avgDataTB_toStudyRTsCorr(ba,:,:))); 
    tt(ba,:)=t{ba}.tstat; 
end

figure();
bar(tt)

set(gca,'FontSize',16)
set(gca,'XTickLabel',{'\delta','\theta','\alpha','\beta','l-\gamma','h-\gamma'})
ylabel('T-Val (accross channel r)')
set(gca,'LineWidth',1)
grid on
legend('all','post','pre','post-pre','location','southeast')
