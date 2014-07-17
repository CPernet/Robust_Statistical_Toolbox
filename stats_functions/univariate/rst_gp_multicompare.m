function [difference,CI,h]=rst_gp_multicompare(gp1,gp2,flag)

% Test for mean differences between groups with an adjustment for
% multiple comparisons using thr max under H0
%
% FORMAT [difference,CI,h]=rab_gp_multicompare(gp1,gp2,flag)
%
% INPUT gp1 and gp2 are 2 matrices n*p and m*p with p the conditions to
%       compare - 
%       flag - indicates to make a figure or not
%
% OUTPUT difference is the mean difference
%        CI the adjusted 95% CI
%        h hypothesis testing
%
% Cyril Pernet 

nboot = 1000;
alpha = 5/100;
if nargin<3
    flag = 1;
end

% multiple comparisons between 2 groups
boot_gp1 = zeros(size(gp1,1),size(gp1,2),nboot); B=1;
while B~=nboot+1
    tmp1 = randi(size(gp1,1),1,size(gp1,1))';
    if length(unique(tmp1)) ~= 1
        boot_gp1(:,:,B) = gp1(tmp1,:);
        B=B+1;
    end
end

boot_gp2 = zeros(size(gp2,1),size(gp2,2),nboot); B=1;
while B~=nboot+1
    tmp2 = randi(size(gp2,1),1,size(gp2,1))';
    if length(unique(tmp2)) ~= 1
        boot_gp2(:,:,B) = gp2(tmp2,:);
        B=B+1;
    end
end

means_gp1 = squeeze(nanmean(boot_gp1,1)); % means over resampled subjects
means_gp2 = squeeze(nanmean(boot_gp2,1));
boot_means1 = mean(nanmeans_gp1,2); % means over all bootstrap
boot_means2 = mean(nanmeans_gp2,2);
Tau1 = (1/(nboot-1))*sum((means_gp1-boot_means1(:,ones(1,nboot))).^2,2); % standard error of each condition
Tau2 = (1/(nboot-1))*sum((means_gp2-boot_means2(:,ones(1,nboot))).^2,2);
H = abs(means_gp1 - means_gp2 - repmat((boot_means1-boot_means2),1,nboot)) ...
    ./ repmat((sqrt(Tau1+Tau2)),1,nboot); % difference between conditions
Hmax = sort(max(H)); % maximum differences over all bootstraps across conditions

u = round((1-alpha)*nboot);
difference = nanmean(gp1,1)' - nanmean(gp2,1)';
bound = Hmax(u).*(sqrt(Tau1+Tau2)); % because of Hmax, this adjusts for multiple testing
CI = NaN(length(difference),2);
CI(:,1) = difference - bound;
CI(:,2) = difference + bound;
h = single(CI(:,1) > 0) + single(CI(:,2) < 0); 

% non corrected version would be
% boot_difference = sort((means_gp1 - means_gp2)');
% l = round((alpha*nboot)/2); u = nboot-l;
% CI(:,1) = boot_difference(l+1,:);
% CI(:,2) = boot_difference(u,:);

% adjustement for large N

if length(gp1)>80 && length(gp2) > 80
    disp('adjustment needed for large N - to complete Wilcox p314')
end

% plot
if flag == 1
    low = round((1-alpha/2)*nboot);
    high = nboot - low;
    figure('Name', 'Pair-wise comparisons')
    set(gcf,'Color','w')
    subplot(1,2,1); plot(nanmean(gp1),'LineWidth',2); hold on
    plot(nanmean(gp2),'r','LineWidth',2); grid on;
    means_gp1 = sort(means_gp1,2); means_gp2 = sort(means_gp2,2);
    fillhandle = patch([1:size(gp1,2) size(gp1,2):-1:1], [means_gp1(:,high); flipud(means_gp1(:,low))], [0 0 1]);
    set(fillhandle,'EdgeColor',[0 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    fillhandle = patch([1:size(gp1,2) size(gp1,2):-1:1], [means_gp2(:,high); flipud(means_gp2(:,low))], [1 0 0]);
    set(fillhandle,'EdgeColor',[1 0 0],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    title('Means','FontSize',14); subplot(1,2,2);
    bar([1:length(difference)],difference); hold on
    errorbar([1:length(difference)],difference,bound,bound,'r','LineWidth',2);
    title(['Difference between groups and ' num2str(100-alpha*100) '% CI'],'FontSize',16); grid on
    ylabel('Difference between groups','FontSize',14);
    xlabel('Conditions','FontSize',14);
    set(gca,'FontSize',12);
    set(gca,'xtick',[1:length(difference)]);
    v=axis; ymin = v(3) + (v(4)-v(3))/10;
end


