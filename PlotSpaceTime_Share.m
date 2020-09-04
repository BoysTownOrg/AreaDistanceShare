% PlotSpaceTime_Share, 8/18/20.  Calls external files drawaxis.m and 
%   arrow.m to overcome to place axes and provide arrows.
close all, clear all
N=6;
x=0:N;
y=-1:(N-1);
x1=1:(N-1);
y1=0:(N-2);
yamax=10; % max time in samples
d=0.30;
xlim=[0,N-d];
ylim=[-1,yamax+0.7];
xtick=1:N;
ytick=0:yamax;
xstra={'1','2','3','4','5'};
xstr={'0','1','2','3','4',''};
ystr={'0','1/2','1','3/2','2','5/2','3','7/2','4','9/2','5'};

hf=figure(1);
set(hf,'Color','white','DefaultLineMarkerSize',6,...
  'Position',[100,100,580,960]);
plot(N+1,N+1); % dummy plot so XLim and YLim can be set before arrow
ha=gca;
set(ha,'DataAspect',[1,1,1],'XLim',xlim,'XTick',xtick,'XTickLabel',xstr,...
  'XColor','w','YColor','w',...
  'YLim',ylim,'YTick',ytick,'YTickLabel',ystr);%,'Box','off');
set(ha,'XGrid','on','YGrid','on');
hold on;
% move axes first before gray patch
drawaxis(ha,'x',0,'movelabel',0); % move x axis to x==0
drawaxis(ha,'y',1,'movelabel',1); % move y axis to y==1
% extra code to handle effects of drawaxis 
mya=gca;
hd=findobj(mya,'type','Text'); % find all text labels on axes
delete(hd); % discard all text in x, y axis labels
% init for gray fill
x0=0.8; % between 0,1
y0=x0-1;
xp=[x0,x0,0,0];
cg=repmat(0.9,1,3);
% gray fill only for TM boundary with measured RF at probe tip
  % draw gray background, patch with slopes of \infty, -1, 1
  yp=[y0,-2*y0+2*(N-2),N-2,-y0];
xp(3)=xp(2)-(yp(3)-yp(2));
xp(4)=xp(3)+(yp(4)-yp(3));
hpa=patch(xp,yp,cg);
hpa.LineStyle='none';
set(gca,'children',flipud(get(gca,'children'))); % invert front to back graphics objects
plot([0,x0],[0,0],'Color','w','LineWidth',1); % paint out x axis in [0,1] include label 0
plot([x0,1],[0,0],'Color',cg,'LineWidth',1); % paint out x axis in [0,1] include label 0
% erase axes before new origin
plot([1,1],[-1,0],'Color','w','LineWidth',1); % paint out x axis in [0,1] include label 0
plot(1,0,'ko','MarkerFaceColor','k'); % redraw bottom circle on y axis at 0
% Create x, y axis labels
fst=12;
text(2,-0.5,'Distance  (units of {\it cT/2} )',...
  'FontSize',fst); 
temp=-0.30;
text(temp,N-1.4,'Time  (units of {\it T} )','Rotation',90,...
  'FontSize',fst); 
text(1.4:(N-1+0.4),repmat(-0.27,1,N-1),xstra); % add x-axis labels
text(0.6+zeros(1,yamax+1),0:yamax,ystr); % fix y axis at 0 time

fsp=9;
sp='g'; % test symbol for variable at each grid point
plotStart(sp,d,x,y,fsp); % for incident diagonal
sw=1;
for ii=1:(N-1)
  yii=2*ii+y1;
  rng=find(yii<=yamax);
  plot1(sp,d,x1(rng),yii(rng),fsp,N);
end
for ii=2:(yamax-1)
  plot2(N,yamax,d,ii,ii-1,ii-2);
end
yvert=[0,ylim(2)];
ytt=0.3;
for jj=2:(N-1)
  plot([jj,jj],yvert,'k:');
  if jj<N-1
    text(jj+0.04,ytt,['m=',int2str(jj)]);
  else
    plot([jj+0.04,jj+0.04],yvert,'k:');
    text(jj+0.08,ytt,['m=',int2str(jj),'=M']);
    text(jj+0.08,ytt+0.3,'Near TM');
  end
end
dt=0.29;
xt=dt;
fst=10;
text(xt+0.2,-1+0.95*dt,'\it Input:  g_r^+[0,0]=1',...
  'HorizontalAlignment','center','FontSize',fst);
text(xt,1-dt,'\it h[0]','HorizontalAlignment','right','FontSize',fst);
text(xt,3-dt,'\it h[1]','HorizontalAlignment','right','FontSize',fst);
text(xt,5-dt,'\it h[2]','HorizontalAlignment','right','FontSize',fst);
text(xt,7-dt,'\it h[3]','HorizontalAlignment','right','FontSize',fst);
text(xt,9-dt,'\it h[4]','HorizontalAlignment','right','FontSize',fst);
text(xt,11-dt,'\it h[5]','HorizontalAlignment','right','FontSize',fst);
% Increase font sizes by a factor of 1.1
h=findall(gcf, '-property', 'fontsize');
hFontSize = cell2mat(get(h,'FontSize'));
newFontSize = hFontSize * 1.1;
set(h,{'FontSize'}, num2cell(newFontSize))
h=findall(gcf, '-property', 'LineWidth');
lw=1;
set(h,{'LineWidth'}, num2cell(lw))

function plotStart(sp,d,x,y,fsp)
m=length(x)-1;
for jj=1:m
  pmid=[x(jj)+2*x(jj+1),y(jj)+2*y(jj+1)]/3;
  p1=[x(jj)+d,y(jj)+d];
  if jj==1
    arrow(p1,pmid,'Type','Line','Length',8,'Width',0);
    DotextInput(x(jj+1),y(jj+1));
  else
    plot(x(jj),y(jj),'ko','MarkerFaceColor','k');
    hold on;
    if jj<m
      arrow(p1,pmid,'Type','Line','Length',8,'Width',0);
      if jj<m-1
        DotextInc(x(jj+1),y(jj+1));
      else
        xn=x(jj+1);
        yn=y(jj+1);
        ss=ssMake(xn-1,yn);
        text(xn,yn,['\it ',sp,'_r^+',ss],'FontSize',fsp,...
          'HorizontalAlignment','right','VerticalAlignment','top');
        text(xn,yn,['\it ',sp,'_r^-',ss],'FontSize',fsp,...
          'HorizontalAlignment','right','VerticalAlignment','bottom');
      end
    else
    end
  end
end

  function DotextInput(xn,yn)
    ss=ssMake(xn,yn);
    text(xn,yn,['\it ',sp,'_l^+',ss],'FontSize',fsp,...
      'HorizontalAlignment','left','VerticalAlignment','bottom');
  end

  function DotextInc(xn,yn)
    ss=ssMake(xn-1,yn);
    text(xn,yn,['\it ',sp,'_r^+',ss],'FontSize',fsp,...
      'HorizontalAlignment','right','VerticalAlignment','top');
    text(xn,yn,['\it ',sp,'_r^-',ss],'FontSize',fsp,...
      'HorizontalAlignment','right','VerticalAlignment','bottom');
    ss=ssMake(xn,yn);
    text(xn,yn,['\it ',sp,'_l^+',ss],'FontSize',fsp,...
      'HorizontalAlignment','left','VerticalAlignment','bottom');
  end
end

function ss=ssMake(xn,yn)
  if rem(yn,2)
    ss=['[\it ',int2str(xn),',',int2str(yn),'/2]'];
  else
    yn2=floor(yn/2);
    ss=['[\it ',int2str(xn),',',int2str(yn2),']'];
  end
end

function plot1(sp,d,x,y,fsp,N)
m=length(x);
for jj=1:m
  plot(x(jj),y(jj),'ko','MarkerFaceColor','k');
  hold on;
  pmid=[x(jj)+2*(x(jj)+1),y(jj)+2*(y(jj)+1)]/3;
  p1=[x(jj)+d,y(jj)+d];
  if jj<m 
    arrow(p1,pmid,'Type','Line','Length',8,'Width',0);
  else
    if m<N-1
      arrow(p1,pmid,'Type','Line','Length',8,'Width',0);
    else
    end
  end
  if jj>1
    if jj==m
      DotextAnechoic(sp,x(jj),y(jj),fsp);
    else
      Dotext(sp,x(jj),y(jj),fsp);
    end
  else
    DotextProbe(sp,x(jj),y(jj),fsp);
  end
end

  function Dotext(sp,xn,yn,fsp)
    ss=ssMake(xn-1,yn);
    text(xn,yn,['\it ',sp,'_r^+',ss],'FontSize',fsp,...
      'HorizontalAlignment','right','VerticalAlignment','top');
    text(xn,yn,['\it ',sp,'_r^-',ss],'FontSize',fsp,...
      'HorizontalAlignment','right','VerticalAlignment','bottom');
    ss=ssMake(xn,yn);
    text(xn,yn,['\it ',sp,'_l^+',ss],'FontSize',fsp,...
      'HorizontalAlignment','left','VerticalAlignment','bottom');
    text(xn,yn,['\it ',sp,'_l^-',ss],'FontSize',fsp,...
      'HorizontalAlignment','left','VerticalAlignment','top');
  end

  function DotextAnechoic(sp,xn,yn,fsp)
    ss=ssMake(xn-1,yn);
    text(xn,yn,['\it ',sp,'_r^+',ss],'FontSize',fsp,...
      'HorizontalAlignment','right','VerticalAlignment','top');
    text(xn,yn,['\it ',sp,'_r^-',ss],'FontSize',fsp,...
      'HorizontalAlignment','right','VerticalAlignment','bottom');
  end
end

function DotextProbe(sp,xn,yn,fsp)
ss=ssMake(xn,yn);
text(xn,yn,['\it ',sp,'_l^+',ss],'FontSize',fsp,...
  'HorizontalAlignment','left','VerticalAlignment','bottom');
text(xn,yn,['\it ',sp,'_l^-',ss],'FontSize',fsp,...
  'HorizontalAlignment','left','VerticalAlignment','top');
end

function plot2(N,yamax,d,ii,xii,yii)
if ii<=N
  for jj=ii:-1:2
    pmid=[xii+2*(xii-1),yii+2*(yii+1)]/3;
    p1=[xii-d,yii+d];
    arrow(p1,pmid,'Type','Line','Length',8,'Width',0);
    hold on;
    xii=xii-1;
    yii=yii+1;
  end
else
  xii=xii-(ii-N); % set max fixed boundary at TM
  yii=yii+(ii-N);
  jj=N;
  while jj>=2 && yii<=yamax
    pmid=[xii+2*(xii-1),yii+2*(yii+1)]/3;
    p1=[xii-d,yii+d];
    arrow(p1,pmid,'Type','Line','Length',8,'Width',0);
    hold on;
    jj=jj-1;
    xii=xii-1;
    yii=yii+1;
  end
end
end


