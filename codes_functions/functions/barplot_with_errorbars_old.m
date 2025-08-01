function barplot_with_errorbars(data,names,s,colors,nanflag)

%% inputs
% data --> data matrix with size n x m , n--> number of observations,
%          m --> no. of variables
% names --> cell containing the names you want to put on the x-axis
% s --> variable telling what to show in error bars
%         if s==1 --> std deviation (default)
%            s==2 --> standard error of the mean

% eg. names = {'Vehicle','LHRH'};

if nargin <3
    s = 1;
    nanflag = 1;
end
if nargin <4
    nanflag = 1;
end
if nargin <5
    nanflag = 1;
end
if nanflag == 1
    mean_ = squeeze(mean(data,'omitnan'));
else
    mean_ = squeeze(mean(data));
end

if s==1
    if nanflag == 1
        std_ = squeeze(std(data,'omitnan'));
    else
        std_ = squeeze(std(data));
    end
elseif s==2

    if nanflag == 1
        std_ = squeeze(std(data,'omitnan')/sqrt(size(data,1)));
    else
        std_ = squeeze(std(data)/sqrt(size(data,1)));
    end

end

if size(size(data)) <= 2

%     figure;
if ~exist('colors',"var")
    bar(mean_)
else
    b = bar(mean_);
    b(1).FaceColor = colors;
end
    hold on
    errorbar(1:length(mean_),mean_,std_,std_, 'LineWidth',2, 'MarkerSize',1,LineStyle = 'none');
%     er.LineStyle = 'none';

    if exist("names","var")
        set(gca,'xtick',1:length(mean_),'xticklabel',names)
    end

else
    b = bar(mean_, 'grouped');
   
    hold on
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(mean_);
    % Get the x coordinate of the bars
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
    % Plot the errorbars
    errorbar(x',mean_,std_,std_,'k', 'LineWidth',2, 'MarkerSize',1,LineStyle = 'none');
    hold off
    if exist("names","var")
        set(gca,'xtick',1:length(mean_),'xticklabel',names)
    end
    
    if exist('colors',"var")
        for bi = 1:length(b)
            b(bi).FaceColor = colors(bi,:);
        end
    end
end
