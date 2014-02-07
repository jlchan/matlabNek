load CFL_Test_DJ;
figure
I = length(CFLV);
J = length(ReV);
nplots = I*J;

maxres = 0;
maxrat = 0;
for k = res.keys
    R = res(k{1});
    maxres = max(maxres,max(R(:)));
    R = res_ratio(k{1});
    maxrat = max(maxrat,max(R(:)));    
end
resvec = {res; ref_res; res_ratio};
maxrvec = [maxres maxres maxrat];

% plot w/out convection
for i = 1:length(CFLV)
    for j = 1:length(ReV)
        str = ['CFL=',num2str(CFLV(i)),', ReScale=',num2str(ReV(j))];
        for figs = 1:3
            figure(figs)
            subplot(I,J,(i-1)*J + j)
            surf(resvec{figs}(str))
            set(gca,'Xticklabel',[])
            set(gca,'Yticklabel',[])
            zlim([0 maxrvec(figs)])
            caxis([0 maxrvec(figs)])
            title(str)
            if (i==length(CFLV))
                xlabel('Order')
            end
            if (j==1)
                ylabel('Elements')
            end
            
        end
    end
end
 print(figure(1),'-depsc','res_agmg_cd_DJ.eps')
 print(figure(2),'-depsc','ref_agmg_cd_DJ.eps')
 print(figure(3),'-depsc','ratio_agmg_cd_DJ.eps')