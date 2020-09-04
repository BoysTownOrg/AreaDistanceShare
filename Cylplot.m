function Cylplot(TubeArea,z,area,xlim,xtick,ms,co)
    z0=[xlim(1),xtick(1)];
    plot(z0,repmat(TubeArea,1,2),':k');
    hold on;
    N=length(z);
    if isempty(co)
      cp='k';
    else
      cp=co;
    end
    plot(repmat(z(1),1,2),[TubeArea,area(1)],'LineStyle',':','Color',cp);
    for kk=1:(N-1)
      dz=[z(kk),z(kk+1)];
      plot(dz,repmat(area(kk),1,2),'LineStyle','-','Color',cp);
      plot(repmat(z(kk+1),1,2),[area(kk),area(kk+1)],...
        'LineStyle','-','Color',cp);
      plot(z(kk),area(kk),[cp,'o'],'MarkerSize',ms);
    end
    plot(z(N),area(N),[cp,'o'],'MarkerSize',ms);
  end