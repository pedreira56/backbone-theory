function figE3(hint)
%FIGE3  Redraw the figure of experiment E3 with the curve of Proposition P8.
%
%   FIGE3() finds E3.mat (see BBFINDMAT) and writes E3_fig_v2.pdf,
%   E3_fig_v2.png and E3_P8_table.txt next to it.  Nothing is simulated.
%   FIGE3(FOLDER) searches FOLDER first.
%
%   Needs only this file and bbFindMat.m.
%
%   WHAT IS ADDED
%   Proposition P8 gives the trade-off of panel (b) in closed form.  For the
%   average of the k nearest donors,
%       individual risk      ->  s2 * (1 + 1/k)
%       W2^2 to the true law ->  tau2 * (1 - 1/sqrt(k))^2      (Gaussian case)
%   where s2 is the mean conditional variance of the attribute and tau2 the
%   mean trace of the conditional covariance of the block.  The curve is
%   drawn for continuous k and the simulated points are placed on top of it,
%   so the agreement is seen rather than asserted.  Panel (a) is unchanged.

if nargin < 1, hint = ''; end
addpath(fileparts(mfilename('fullpath')));
[matFile, outDir] = bbFindMat('E3', hint);
out = load(matFile);

gA = 1;                       % regime shown in the figure
jL = numel(out.nGrid);        % largest sample size

f = figure('Position',[100 100 1020 420]);

% ---- (a) unchanged -----------------------------------------------------
subplot(1,2,1); hold on; box on
wV = squeeze(out.w2K(gA,:,1));
plot(out.nGrid([1 end]), out.tau2E(gA)*[1 1], '-', ...
     'Color',[0.74 0.74 0.74],'LineWidth',5, ...
     'DisplayName','floor  E[\tau^2(X_B)]   (Proposition P4)');
plot(out.nGrid, wV, 'o-','Color',[0.85 0.33 0.10], ...
     'MarkerSize',4.5,'LineWidth',1.3,'MarkerFaceColor','w', ...
     'DisplayName',sprintf('value transfer   (slope %.2f / %.2f)', ...
                   out.slopeW(gA), -2*out.beta(gA)/out.dB(gA)));
plot(out.nGrid, out.w2M(gA,:), 'o','Color',[0.00 0.45 0.74], ...
     'MarkerSize',5.5,'LineWidth',1.4,'MarkerFaceColor','w', ...
     'DisplayName','mean transfer');
set(gca,'XScale','log','YScale','log','FontSize',9);
xlabel('n   (donor block size)');
ylabel('W_2^2 to the true conditional law');
ylim([min(wV)/4, out.tau2E(gA)*6]);
xlim([out.nGrid(1)/1.6, out.nGrid(end)*1.6]);
text(-0.02,1.05,'(a)','Units','normalized','FontSize',11, ...
     'FontWeight','bold','HorizontalAlignment','left');
legend('Location','southwest','FontSize',8,'Box','off');

% ---- (b) simulated points on the curve of Proposition P8 ---------------
subplot(1,2,2); hold on; box on
rk = reshape(out.riskK(gA,jL,:),[],1);
wk = reshape(out.w2K(gA,jL,:),[],1);
rM = out.riskM(gA,jL);  wM = out.w2M(gA,jL);
s2   = rM;                    % mean conditional variance of the attribute
tau2 = out.tau2E(gA);         % mean trace of the conditional covariance
kc   = logspace(0,4,400);
rTh  = s2*(1 + 1./kc);
wTh  = tau2*(1 - 1./sqrt(kc)).^2;
spanX = max([rk; rM]) - min([rk; rM]);
spanY = max([wk; wM]) - min([wk; wM]);
plot(rTh, wTh, '--','Color',[0.35 0.35 0.35],'LineWidth',1.3, ...
     'DisplayName','Proposition P8');
plot(rk, wk, 'o','Color',[0.49 0.18 0.56],'MarkerSize',5.5,'LineWidth',1.3, ...
     'MarkerFaceColor','w','DisplayName','average of k nearest donors');
plot(rM, wM, 's','Color',[0.00 0.45 0.74],'MarkerSize',9,'LineWidth',1.5, ...
     'MarkerFaceColor','w', ...
     'DisplayName','oracle conditional mean  (k \rightarrow \infty)');
set(gca,'FontSize',9);
xlabel('individual quadratic risk');
ylabel('W_2^2 to the true conditional law');
xlim([min([rk; rM])-0.12*spanX, max([rk; rM])+0.14*spanX]);
ylim([min([wk; wM])-0.10*spanY, max([wk; wM])+0.10*spanY]);
for i = 1:numel(out.kList)
    if i == 1, ha = 'right'; dx = -0.015*spanX; else, ha = 'left'; dx = 0.015*spanX; end
    text(rk(i)+dx, wk(i), sprintf('k=%d',out.kList(i)), 'FontSize',7.5, ...
         'HorizontalAlignment',ha,'VerticalAlignment','middle');
end
text(-0.02,1.05,'(b)','Units','normalized','FontSize',11, ...
     'FontWeight','bold','HorizontalAlignment','left');
% The trade-off curve runs from the top left to the bottom right, so the
% upper-right corner is the only one no curve passes through.
legend('Location','northeast','FontSize',8,'Box','off');

exportgraphics(f, fullfile(outDir,'E3_fig_v2.pdf'),'ContentType','vector');
exportgraphics(f, fullfile(outDir,'E3_fig_v2.png'),'Resolution',200);

% ---- table: prediction against simulation ------------------------------
kL  = out.kList(:);
rPr = s2*(1 + 1./kL);
wPr = tau2*(1 - 1./sqrt(kL)).^2;
hdr = sprintf('  %6s %12s %12s %12s %12s\n','k','risk P8','risk sim','W2^2 P8','W2^2 sim');
fprintf('\n=== E3: Proposition P8 against simulation, n = %d ===\n',out.nGrid(jL));
fprintf('  s2 = %.6f   tau2 = %.6f\n%s',s2,tau2,hdr);
fid = fopen(fullfile(outDir,'E3_P8_table.txt'),'w');
fprintf(fid,'E3 - Proposition P8 against simulation, n = %d\n',out.nGrid(jL));
fprintf(fid,'  s2 = %.6f   tau2 = %.6f\n\n%s',s2,tau2,hdr);
for i = 1:numel(kL)
    row = sprintf('  %6d %12.4f %12.4f %12.4f %12.4f\n',kL(i),rPr(i),rk(i),wPr(i),wk(i));
    fprintf('%s',row); fprintf(fid,'%s',row);
end
fclose(fid);
fprintf('\nWrote E3_fig_v2.pdf, E3_fig_v2.png, E3_P8_table.txt to:\n  %s\n',outDir);
end
