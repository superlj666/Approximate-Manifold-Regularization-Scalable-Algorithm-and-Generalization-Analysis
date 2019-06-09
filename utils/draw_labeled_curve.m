function y=draw_labeled_curve(dataset)
    str = ['./result/labeled_curve_', dataset];
    load(str);
    fig=figure('Position',  [232   246   500   420]);
    
    x_list=1:7;
    errorbar(x_list,res_RLS(:,1),res_RLS(:,2), 'k--o','Markersize', 10, 'LineWidth',2.5);
    hold on;
    errorbar(x_list,res_LapRLS_pcg(:,1),res_LapRLS_pcg(:,2),'b-.^','Markersize', 10, 'LineWidth',2.5);
    errorbar(x_list,res_nystrom_pcg(:,1),res_nystrom_pcg(:,2),'r-x','Markersize', 10, 'LineWidth',2.5);

    ax = gca;
    ax.ActivePositionProperty = 'outerposition';
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + 3.4*ti(1);
    bottom = outerpos(2) + 4*ti(1);
    ax_width = outerpos(3) - 3.5*ti(1);
    ax_height = outerpos(4) - 4.5*ti(1);
    ax.Position = [left bottom ax_width ax_height];
    
    %title(dataset, 'Interpreter', 'none');
    xticklabels({'1%', '', '4%', '' ,'16%', '', '64%'});
    %yticklabels({'0.95', '1.00', '1.05', '1.10' ,'1.15'});
    grid on
    set(gca,'XLim',[0.5 7.5])
    legend({'RLS', 'LapRLS', 'Nystrom'});
    ylabel('RMSE','FontSize',25);
    xlabel('partition of labeled data','FontSize',25);
    set(gca,'FontSize',25, 'Fontname', 'Times New Roman');
    hold off;
    
    print(fig,['.\result\labeled_curve_', dataset],'-depsc')
end