function simpleplot(fn, mode, tit)
  % This function provides a simple command for plotting figure
  if nargin<3
    tit='';
  end
  
  [xk, emesh, dos]=read_data(fn);
  
  fid=fopen('klabels.dat');
  xtic_data=textscan(fid, ' %f %s %s');
  xtics=xtic_data{1};
  nn=length(xtics);
  
  xticlabels=cell(nn+1,1);
  for ii=1:nn
    if (ii==1)
      xticlabels{ii}=xtic_data{2}{ii};
    else
      if (xtic_data{2}{ii}==xtic_data{3}{ii-1})
        xticlabels{ii}=xtic_data{2}{ii};
      else
        xticlabels{ii}=strcat(xtic_data{3}{ii-1},'/',xtic_data{2}{ii});
      end
    end
  end
  
  xticlabels{nn+1}=xtic_data{3}{nn};
  xtics(nn+1)=max(xk);
  
  fig=figure();
  set(fig, 'PaperUnits', 'centimeters');
  set(fig, 'PaperPosition', [0 0 8 6]);
  set(fig, 'PaperPositionMode', 'Manual');
  
%  [xx, yy]=meshgrid(xk, emesh);
%  xkf=zeros(length(xk)*2-1);
%  for ii=1:length(xk)*2-1
%      if (mod(ii,2)==1)
%          xkf(ii)=xk((ii+1)/2);
%      else
%          xkf(ii)=(xk(ii/2)+xk(ii/2+1))/2;
%      end
%  end
%  emeshf=zeros(length(emesh)*2-1);
%  for ii=1:length(emesh)*2-1
%      if (mod(ii,2)==1)
%          emeshf(ii)=emesh((ii+1)/2);
%      else
%          emeshf(ii)=(emesh(ii/2)+emesh(ii/2+1))/2;
%      end
%  end
%  [xxf, yyf]=meshgrid(xkf, emeshf);
%  dosf=interp2(xx, yy, dos, xxf, yyf, 'spline');
  dosmax=max(max(dos));

  if (mode==0)
    pcolor(xk, emesh, dos);
  else
    pcolor(xk, emesh, log(abs(dos)));
  end
  dy=round((max(emesh)-min(emesh))/5.0,1);
  ymax=floor(max(emesh)/dy)*dy;
  ymin=ceil(min(emesh)/dy)*dy;
  
  r=linspace(1.0,0.0,256);
  g=linspace(1.0,0.0,256);
  b=linspace(1.0,1.0,256);
  ctmp=[r;g;b];
  cmap=transpose(ctmp);
  colormap(cmap);
  caxis([0 dosmax/2]);
  shading interp;
  hold on;
  title(tit);
  set(gca, 'FontName', 'Times');
  set(gca, 'FontSize', 10);
  set(gca, 'Xtick', xtics, 'Ytick', [ymin:dy:ymax],'XTickLabel', xticlabels);
  plot([min(xk) max(xk)], [0 0], 'Color', [0.5,0.0,0.0]);
%  plot([min(xk) max(xk)], [-0.2 -0.2], 'Color', [0.5,0.0,0.0], 'LineStyle', '--');
%  plot([min(xk) max(xk)], [-0.24 -0.24], 'Color', [0.5,0.0,0.0], 'LineStyle', ':');
%s  plot([min(xk) max(xk)], [-0.28 -0.28], 'Color', [0.5,0.0,0.0], 'LineStyle', '-.');
  for ii=2:nn
    plot([xtics(ii) xtics(ii)], [min(emesh) max(emesh)], 'Color', [0.5,0.0,0.0]);
  end
  xlabel('K-points');
  ylabel('Energy (eV)');
%  colorbar;
  hold off;
  print(gcf, 'plot.png', '-dpng', '-r600');
end
