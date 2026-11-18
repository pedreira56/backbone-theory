function out = runE4(quick)
%RUNE4  Driver for experiment E4 (Proposition P5).
%
%   OUT = RUNE4()      full run   (use this for the paper)
%   OUT = RUNE4(true)  quick run  (smoke test)
%
%   Layout: put every .m file in  ...\Exp4\Cod\  and results appear in
%   ...\Exp4\Result\ , created automatically.  See BBOUTDIR.
%
%   Outputs: E4.mat, E4_fig.pdf/.png, E4_table.txt
%
%   Panel (a) is a calibration plot: the empirical distribution of the
%   two-sample test's p-values against the uniform.  A test with power bends
%   above the diagonal.  P5 predicts that this one cannot, and the panel is
%   how the prediction is checked rather than asserted.
%   Panel (b) shows the one thing that does differ between the two
%   populations, and the design never records it.
%   Panel (c) shows the price: the reconstruction settles on the value the
%   tail probability takes under the product coupling even though the data
%   came from the dependent population.

if nargin < 1 || isempty(quick), quick = false; end

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
outDir = bbOutDir(mfilename('fullpath'));
fprintf('[E4] code   : %s\n[E4] results: %s\n\n', codeDir, outDir);

cfg = struct();
if quick
    cfg.nGrid = round(logspace(2,3.7,5)); cfg.nRecv = 4000;
    cfg.nTest = 150; cfg.nPerm = 100; cfg.nRep = 60; cfg.nOracle = 2e5;
else
    cfg.nGrid = round(logspace(2,5,8));   cfg.nRecv = 40000;
    cfg.nTest = 250; cfg.nPerm = 200; cfg.nRep = 400; cfg.nOracle = 2e6;
end

out = E4_nonident(cfg);
out.outDir = outDir;

disp(' ');
disp('=== E4: block-wise summaries under the two populations ===');
disp(out.summary);
fprintf('\n=== E4: rejection rates of the two-sample test ===\n');
fprintf('  %10s %14s %14s\n','nominal','block 1','block 2');
for i = 1:numel(out.levels)
    fprintf('  %10.2f %14.4f %14.4f\n',out.levels(i),out.rejection(i,1),out.rejection(i,2));
end
fprintf('\n=== E4: joint tail probability ===\n');
fprintf('  under P-perp (Cond_CI holds)  : %.5f\n',out.tailPerp);
fprintf('  under P-prime (Cond_CI fails) : %.5f   (ratio %.2f)\n', ...
        out.tailPrime, out.tailPrime/out.tailPerp);
fprintf('  reconstruction, data from P-prime, at n = %d : %.5f\n', ...
        out.nGrid(end), out.tailRec(end));
fprintf('\n=== E4: cross-block correlation of the imputed pair ===\n');
fprintf('  (unconditional, so the reference is the unconditional value,\n');
fprintf('   not the conditional parameter theta = %.2f)\n',out.thetaCond);
fprintf('  true under P-prime            : %.4f\n',out.corrPrime);
fprintf('  true under P-perp             : %.4f\n',out.corrPerp);
fprintf('  reconstruction, by n:\n');
for j = 1:numel(out.nGrid)
    fprintf('    n = %8d : %+.4f\n',out.nGrid(j),out.corrRec(j));
end

f = figure('Position',[100 100 1200 380]);

% ---- (a) calibration of the two-sample test ---------------------------
subplot(1,3,1); hold on; box on
plot([0 1],[0 1],'-','Color',[.5 .5 .5],'LineWidth',1.1, ...
     'DisplayName','uniform (no power)');
cols = [0.00 0.45 0.74; 0.85 0.33 0.10];
for m = 1:2
    x = sort(out.pvals(:,m));
    plot(x,(1:numel(x))'/numel(x),'-','Color',cols(m,:),'LineWidth',1.4, ...
         'DisplayName',sprintf('observed, block %d',m));
end
set(gca,'FontSize',9); axis square
xlabel('nominal level'); ylabel('empirical rejection rate');
xlim([0 1]); ylim([0 1]);
text(-0.02,1.06,'(a)','Units','normalized','FontSize',11,'FontWeight','bold');
legend('Location','southeast','FontSize',8,'Box','off');

% ---- (b) the conditional joint the design never sees ------------------
% Two overlaid scatters of nearly identical clouds read as a single cloud,
% and the tilt -- the only thing that differs between the two populations --
% becomes invisible.  The exact 90 per cent probability ellipses are drawn
% on top: both laws are Gaussian conditionally on the backbone, so the
% ellipses are available in closed form and carry the message that the point
% clouds cannot.
subplot(1,3,2); hold on; box on
cols = [0.00 0.45 0.74; 0.85 0.33 0.10];
nm   = {'Cond\_CI holds','Cond\_CI fails'};
for kk = 1:2
    scatter(out.cond(kk).q1, out.cond(kk).q2, 3, cols(kk,:), 'filled', ...
            'MarkerFaceAlpha',0.10, 'HandleVisibility','off');
end
th = linspace(0,2*pi,200);
for kk = 1:2
    L = chol(out.cond(kk).Sig,'lower');
    E = out.cond(kk).mu.' + sqrt(chi2inv(0.90,2)) * L * [cos(th); sin(th)];
    plot(E(1,:),E(2,:),'-','Color',cols(kk,:),'LineWidth',1.8, ...
         'DisplayName',nm{kk});
end
set(gca,'FontSize',9); axis square
xlabel('exclusive attribute, block 1');
ylabel('exclusive attribute, block 2');
text(-0.02,1.06,'(b)','Units','normalized','FontSize',11,'FontWeight','bold');
legend('Location','northwest','FontSize',8,'Box','off');

% ---- (c) the cost --------------------------------------------------
% Error bars are the binomial standard error of the estimated probability.
% Without them the smallest sample sizes look like failures to converge when
% part of what is seen there is estimation noise and part is a genuine
% neighbourhood effect; the bars separate the two.
subplot(1,3,3); hold on; box on
plot(out.nGrid([1 end]), out.tailPrime*[1 1],'--','Color',[0.85 0.33 0.10], ...
     'LineWidth',1.5,'DisplayName','truth, data-generating law');
plot(out.nGrid([1 end]), out.tailPerp*[1 1],'--','Color',[0.00 0.45 0.74], ...
     'LineWidth',1.5,'DisplayName','value under the product coupling');
errorbar(out.nGrid, out.tailRec, 2*out.tailSE,'o-','Color',[0.25 0.25 0.25], ...
     'MarkerSize',4.5,'LineWidth',1.2,'MarkerFaceColor','w','CapSize',3, ...
     'DisplayName','reconstruction');
set(gca,'XScale','log','FontSize',9); axis square
xlabel('n   (donor block size)');
ylabel('P(both attributes exceed their upper quintile)');
ylim([0, out.tailPrime*1.35]);
xlim([out.nGrid(1)/1.6, out.nGrid(end)*1.6]);
text(-0.02,1.06,'(c)','Units','normalized','FontSize',11,'FontWeight','bold');
legend('Location','southeast','FontSize',8,'Box','off');

exportgraphics(f, fullfile(outDir,'E4_fig.pdf'), 'ContentType','vector');
exportgraphics(f, fullfile(outDir,'E4_fig.png'), 'Resolution',200);

save(fullfile(outDir,'E4.mat'),'-struct','out');
fid = fopen(fullfile(outDir,'E4_table.txt'),'w');
fprintf(fid,'E4 - Proposition P5 made operational\n\n');
fprintf(fid,'Block-wise summaries under the two populations\n');
fprintf(fid,'%s\n',evalc('disp(out.summary)'));
fprintf(fid,'\nRejection rates of the two-sample energy test (%d replications)\n',out.cfg.nRep);
fprintf(fid,'  %10s %14s %14s\n','nominal','block 1','block 2');
for i = 1:numel(out.levels)
    fprintf(fid,'  %10.2f %14.4f %14.4f\n',out.levels(i),out.rejection(i,1),out.rejection(i,2));
end
fprintf(fid,'  mean p-value: block 1 = %.4f, block 2 = %.4f  (uniform: 0.5)\n', ...
        mean(out.pvals(:,1)),mean(out.pvals(:,2)));
fprintf(fid,'\nJoint tail probability\n');
fprintf(fid,'  under P-perp                 : %.5f\n',out.tailPerp);
fprintf(fid,'  under P-prime                : %.5f  (ratio %.3f)\n', ...
        out.tailPrime,out.tailPrime/out.tailPerp);
fprintf(fid,'  reconstruction, by n (with binomial standard error):\n');
for j = 1:numel(out.nGrid)
    fprintf(fid,'    n = %8d : %.5f  (s.e. %.5f)\n', ...
            out.nGrid(j),out.tailRec(j),out.tailSE(j));
end
fprintf(fid,'\nCross-block correlation of the imputed pair\n');
fprintf(fid,'  (unconditional; conditional parameter is theta = %.2f)\n',out.thetaCond);
fprintf(fid,'  true under P-prime = %.4f ; true under P-perp = %.4f\n', ...
        out.corrPrime,out.corrPerp);
for j = 1:numel(out.nGrid)
    fprintf(fid,'    n = %8d : %+.4f\n',out.nGrid(j),out.corrRec(j));
end
fclose(fid);

fprintf('\nWrote E4.mat, E4_fig.pdf, E4_fig.png, E4_table.txt to:\n  %s\n', outDir);
end
