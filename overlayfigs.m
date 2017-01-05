%overlay two figures in one

%First load the figures:
fig1 = open('ifirstPRC_analyticplusnum.fig');
fig2 = open('iEPRC_analyticplusnum.fig');

%Get the axes objects from the figures
ax1 = get(fig1, 'Children');
ax2 = get(fig2, 'Children');

%Now copy the hangle graphics objects from ax2 to ax1. The loop isn't
%neccesary if your figures only have a single axes
for i = 1 : numel(ax2) 
   ax2Children = get(ax2(i),'Children');
   copyobj(ax2Children, ax1(i));
end