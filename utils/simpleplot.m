function simpleplot(fn, mode)
  % This function provides a simple command for plotting figure
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
  set(fig, 'PaperPosition', [0 0 12 8]);
  set(fig, 'PaperPositionMode', 'Manual');
  if (mode==0)
    pcolor(xk, emesh, dos);
  else
    pcolor(xk, emesh, log(dos));
  end
  colormap('copper');
  shading interp;
  hold on;
  set(gca, 'FontName', 'Helvetica');
  set(gca, 'FontSize', 8);
  set(gca, 'Xtick', xtics, 'Ytick', [-1:0.2:1],'XTickLabel', xticlabels);
  plot([min(xk) max(xk)], [0 0], 'Color', [0.5,0.0,0.0]);
  for ii=2:nn
    plot([xtics(ii) xtics(ii)], [min(emesh) max(emesh)], 'Color', [0.5,0.0,0.0]);
  end
  xlabel('K-points');
  ylabel('Energy (eV)');
%  colorbar;
  hold off;
  saveas(gcf, 'plot.tif', 'tiffn');
end
