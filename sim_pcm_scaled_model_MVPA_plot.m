%%
close all

squeeze(mean(Acc,3));
MEAN = squeeze(mean(mean(Acc,3)));
SEM = squeeze(nansem(mean(Acc,3)));

figure('Name', 'PCM_MVPA', 'Position', [100, 100, 1500, 600], 'Color', [1 1 1]);

h = errorbar(repmat((1:NbSteps)',1,4)+repmat(0:.1:.3,NbSteps,1),MEAN',SEM');
set(h(1), 'color', 'k', 'LineWidth', 1.2)
set(h(2), 'color', 'k', 'Linestyle', '--', 'LineWidth', 1.2)
set(h(3), 'color', 'b', 'LineWidth', 1.2)
set(h(4), 'color', 'b', 'Linestyle', '--', 'LineWidth', 1.2)

set(gca, 'xtick', 1:NbSteps, 'xticklabel', theta(:,1)./theta(:,2))
ax = axis;
axis([0 NbSteps+1 ax(3) ax(4)])

ylabel('accuracy')
xlabel('theta 1 / theta 2')

legend({'No scaling','Img scaling: Z-score','Feat scaling: mean centering',...
    'Img scaling: Z-score ; Feat scaling: mean centering'}, 'Location','SouthEast')

text(NbSteps-3, ax(3)+(ax(4)-ax(3))*.5,...
    sprintf(' Nb vox = %i\n Nb subj = %i\n Var_{sig} ~ N(1,1)\n Var_{noise} = 1.5', NbVox, numSim))

print(gcf, fullfile(Fig_dir, ['PCM_MVPA_' datestr(now, 'yyyy_mm_dd_HH_MM') '.tif']), '-dtiff')

