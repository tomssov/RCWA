for i = 1:1:loc.nlambda
    
    contourf(nmx,nmz - nmgtop, sqrt(squeeze(normE2l(:,:,i))),200,'linestyle','none');
    %pause(0.001);
    axis equal;
    xlim([nmx(1),nmx(end)]);
    caxis([0,100]);
    colorbar;
    title(strcat('\lambda = ',num2str(nmlambda(i))))
    DrawOverlay(loc);
    
    frame = getframe(2);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if i == 1;
          imwrite(imind,cm,'test.gif','gif', 'Loopcount',inf,'DelayTime',0);
      else
          imwrite(imind,cm,'test.gif','gif','WriteMode','append','DelayTime',0);
      end
    
end

