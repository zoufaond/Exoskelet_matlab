function [time_training, data_training,time_simulation,data_simulation] = motion_polyfit(motion_full,plot_polyfits,numdata_training,numdata_simulation)
time = motion_full(:,1)-motion_full(1,1);
time_training = linspace(0,time(end,1),numdata_training);
time_simulation = linspace(0,time(end,1),numdata_simulation);
mot_eul = motion_full(:,8:end)*pi/180;
data_training = zeros(numdata_training,11);
data_simulation = zeros(numdata_simulation,11);
for jnt = 1:length(mot_eul(1,:))
    % p = polyfit(time,mot_eul(:,jnt),4);
    % polyvals(:,jnt) = p;
    % data(:,jnt) = polyval(p,time_new);
    data_training(:,jnt) = interp1(time,mot_eul(:,jnt),time_training,'spline');
    data_simulation(:,jnt) = interp1(time_training,data_training(:,jnt),time_simulation,'spline');
    if strcmp(plot_polyfits,'plot_polyfits')
        figure
        plot(time,mot_eul(:,jnt),'o')
        hold on
        plot(time_simulation,data_simulation(:,jnt))
        hold on 
        plot(time_training,data_training(:,jnt))
        legend('osim','simulation','training')
    end
end

end