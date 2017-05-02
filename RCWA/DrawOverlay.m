function DrawOverlay(varargin)
%ABSORPTION_SPECTRUM Summary of this function goes here
%   Detailed explanation goes here

loc = DefaultLoc;
loc = ApplyVarargin(loc,varargin);
shift = 0*loc.nmLx/4;

nmLz = loc.nmLp+loc.nmLi+loc.nmLn;

    line([-loc.nmLx+shift,loc.nmLx+shift],[nmLz+loc.nmLw,nmLz+loc.nmLw],'color',[0,0,0],'LineWidth',1.2); % Window/Air
    line([-loc.nmLx+shift,loc.nmLx+shift],[nmLz,nmLz],'color',[0,0,0],'LineWidth',1.2); % Junction/Window
    line([-loc.nmLx+shift,loc.nmLx+shift],[0,0],'color',[0,0,0],'LineWidth',1.2); % Grating region/Junction
    
    if(loc.type == 0)
    line([-loc.nmLx+shift,-loc.zeta*loc.nmLx/2+shift],[-loc.nmda+loc.nmLg,-loc.nmda+loc.nmLg],'color',[0,0,0],'LineWidth',1.2); % Top of metal
    line([-loc.zeta*loc.nmLx/2+shift,-loc.zeta*loc.nmLx/2+shift],[-loc.nmda+loc.nmLg,-loc.nmda],'color',[0,0,0],'LineWidth',1.2); % Side of metal
    line([-loc.zeta*loc.nmLx/2+shift,loc.zeta*loc.nmLx/2+shift],[-loc.nmda+loc.nmLg,-loc.nmda+loc.nmLg],'color',[0,0,0],'LineWidth',1.2); % Top of metal
    line([loc.zeta*loc.nmLx/2+shift,loc.zeta*loc.nmLx/2+shift],[-loc.nmda+loc.nmLg,-loc.nmda],'color',[0,0,0],'LineWidth',1.2);
    else
       x = -loc.nmLx:1:loc.nmLx;
       hold on;
       plot(x,gx(x,loc)-loc.nmda,'-black');
       hold off;
    end
    line([-loc.nmLx/2+shift,loc.nmLx/2+shift],[-loc.nmda,-loc.nmda],'color',[0,0,0],'LineWidth',1.2);
    line([-loc.nmLx+shift,loc.nmLx+shift],[-loc.nmda-loc.nmLm,-loc.nmda-loc.nmLm],'color',[0,0,0],'LineWidth',1.2);    


end

