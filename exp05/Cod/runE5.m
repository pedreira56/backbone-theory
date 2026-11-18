function out = runE5(quick)
%RUNE5  Driver for experiment E5 (Theorem T4).
%
%   OUT = RUNE5()      full run   (use this for the paper)
%   OUT = RUNE5(true)  quick run  (smoke test)
%
%   Layout: put every .m file in  ...\Exp5\Cod\  and results appear in
%   ...\Exp5\Result\ , created automatically.  See BBOUTDIR.
%
%   Outputs: E5.mat, E5_fig.pdf/.png, E5_table.txt
%
%   THE TWO FAMILIES OF REGIMES
%   The first varies the number of blocks with all anchoring constants equal,
%   and tests the rate.  The second fixes four blocks and redistributes the
%   anchoring constants among them, holding or changing their squared sum.
%   Two entries of that family are chosen to have nearly the same sum through
%   very different arrangements: if the error depends on the constants only
%   through the sum, as Theorem T4 asserts, those two must land on top of
%   each other.  That is a sharper test than merely observing that the error
%   grows with the number of blocks, which almost any bound would predict.

if nargin < 1 || isempty(quick), quick = false; end

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
outDir = bbOutDir(mfilename('fullpath'));
fprintf('[E5] code   : %s\n[E5] results: %s\n\n', codeDir, outDir);

base = {'d',2,'b',2,'pm',3,'meanType','smooth','beta',1, ...
        'sigma0',0.5,'hetero',0,'rhoWithin',0.3,'theta',0};

R = struct('name',{},'params',{});
% family 1: number of blocks, equal constants
for M = [2 3 4 6]
    R(end+1) = struct('name',sprintf('M = %d, equal',M), ...
        'params',{[base {'M',M,'blockScale',ones(1,M)}]});                %#ok<AGROW>
end
% family 2: four blocks, constants redistributed
sc = {[1 1 1 1], [1 2 2 2], [1 1 1 4], [1 3 1 2], [1 .5 .5 .5]};
for i = 1:numel(sc)
    R(end+1) = struct('name',sprintf('M = 4, scales (%g %g %g)', ...
        sc{i}(2),sc{i}(3),sc{i}(4)), ...
        'params',{[base {'M',4,'blockScale',sc{i}}]});                    %#ok<AGROW>
end

cfg = struct();
if quick
    cfg.nGrid = round(logspace(2,3.7,6)); cfg.repMax = 10; cfg.nRecv = 1000;
else
    cfg.nGrid = round(logspace(2,5,10));  cfg.repMax = 40; cfg.nRecv = 2000;
end

out = E5_coupling(R,cfg);
out.outDir = outDir;

disp(' ');
disp('=== E5: rate and constant structure of Theorem T4 ===');
disp('    discr_per_Lambda2 must agree across the four-block regimes');
disp('    ratio_to_bound must not exceed one');
disp(out.tbl);

f = figure('Position',[100 100 1000 400]);
isFam1 = contains(string({R.name}),'equal') | ~contains(string({R.name}),'scales');

% ---- (a) the rate, one curve per number of blocks ---------------------
subplot(1,2,1); hold on; box on
g1   = find(isFam1);
cols = lines(numel(g1));
for i = 1:numel(g1)
    g = g1(i);
    plot(out.nGrid, out.W2sum(g,:), 'o-','Color',cols(i,:),'MarkerSize',4, ...
         'LineWidth',1.2,'MarkerFaceColor','w', ...
         'DisplayName',sprintf('M = %d   (slope %.2f)',out.M(g),out.slope(g)));
end
set(gca,'XScale','log','YScale','log','FontSize',9);
xlabel('n   (size of each donor block)');
ylabel('W_2^2 to the product coupling');
text(-0.02,1.05,'(a)','Units','normalized','FontSize',11,'FontWeight','bold');
legend('Location','southwest','FontSize',8,'Box','off');

% ---- (b) the constant structure ---------------------------------------
% Every regime, from both families, on one plot: the discrepancy at the
% largest sample size against the sum of squared anchoring constants.  T4
% predicts a straight line through the origin, and in particular predicts
% that two regimes with the same abscissa coincide however differently their
% constants are distributed among blocks.
subplot(1,2,2); hold on; box on
xx = out.lam2; yy = out.W2sum(:,end);
sl = xx\yy;
xl = [0, max(xx)*1.08];
plot(xl, sl*xl, '-','Color',[.45 .45 .45],'LineWidth',1.2, ...
     'DisplayName','proportional fit through the origin');
plot(xx(g1), yy(g1), 'o','Color',[0.00 0.45 0.74],'MarkerSize',7, ...
     'LineWidth',1.4,'MarkerFaceColor','w', ...
     'DisplayName','varying the number of blocks');
g2 = setdiff(1:numel(R), g1);
plot(xx(g2), yy(g2), 's','Color',[0.85 0.33 0.10],'MarkerSize',7, ...
     'LineWidth',1.4,'MarkerFaceColor','w', ...
     'DisplayName','four blocks, constants redistributed');
set(gca,'FontSize',9);
xlabel('\Sigma_l \Lambda_l^2   (sum over donor blocks)');
ylabel(sprintf('W_2^2 to the product coupling  (n = %d)',out.nGrid(end)));
xlim(xl); ylim([0, max(yy)*1.10]);
text(-0.02,1.05,'(b)','Units','normalized','FontSize',11,'FontWeight','bold');
legend('Location','northwest','FontSize',8,'Box','off');

exportgraphics(f, fullfile(outDir,'E5_fig.pdf'), 'ContentType','vector');
exportgraphics(f, fullfile(outDir,'E5_fig.png'), 'Resolution',200);

fprintf('\n=== E5: constant structure at n = %d ===\n',out.nGrid(end));
fprintf('  proportional fit slope           = %.3e\n',sl);
fprintf('  relative spread of the ratio     = %.2f %%\n', ...
        100*std(out.ratio)/mean(out.ratio));
fprintf('  largest ratio to the T4 majorant = %.4f  (must be <= 1)\n', ...
        max(out.tight));

save(fullfile(outDir,'E5.mat'),'-struct','out');
fid = fopen(fullfile(outDir,'E5_table.txt'),'w');
fprintf(fid,'E5 - rate and constant structure of Theorem T4\n\n');
fprintf(fid,'%s\n',evalc('disp(out.tbl)'));
fprintf(fid,'\nConstant structure at n = %d\n',out.nGrid(end));
fprintf(fid,'  proportional fit slope           = %.4e\n',sl);
fprintf(fid,'  relative spread of the ratio     = %.2f %%\n', ...
        100*std(out.ratio)/mean(out.ratio));
fprintf(fid,'  largest ratio to the T4 majorant = %.4f\n',max(out.tight));
fclose(fid);

fprintf('\nWrote E5.mat, E5_fig.pdf, E5_fig.png, E5_table.txt to:\n  %s\n', outDir);
end
