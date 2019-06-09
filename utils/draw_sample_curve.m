function y=draw_sample_curve(str)
    load(str);
    fig=figure;
    x_list=1:7;
    errorbar(x_list,res_RLS(:,1),res_RLS(:,2), 'k--o','LineWidth',1);
    hold on;
    errorbar(x_list,res_LapRLS_pcg(:,1),res_LapRLS_pcg(:,2),'b-.^','LineWidth',1);
    errorbar(x_list,res_nystrom_pcg(:,1),res_nystrom_pcg(:,2),'r-x','LineWidth',1);

    max_level=max(res_RLS(:,1));
    min_level=min(res_nystrom_pcg(:,1));
    step=max_level-min_level;
    
    xticklabels({'1%', '2%', '4%', '8%' ,'16%', '32%', '64%'});
    grid on
    set(gca,'XLim',[0.5 7.5])
    set(gca,'YLim',[min_level-0.5*step max_level+1.2*step])
    legend({'RLS', 'LapRLS CG&PCG', 'Nystrom LapRLS CG&PCG'}, 'FontSize',12);
    ylabel('RMSE');
    xlabel('partition of sampled data');
    set(gca,'FontSize',20,'Fontname', 'Times New Roman');
    hold off;

    print(fig,str,'-depsc')
end