function [ xxx ] = fillsc( xx, sc )
%fillsc fills supercell defined as sc(3,3) with unit cell 
%   Detailed explanation goes here
  nn=det(sc);    % Number of unit cells in this super cell
  ii=0;
  xxlat=zeros(nn,3);
  for ix=0:nn-1
    for iy=0:nn-1
      for iz=0:nn-1
        xlat=[ix,iy,iz];
        tlat=xlat/sc;
        tlat=tlat-floor(tlat);
        if (ii==0)
          ii=ii+1;
          xxlat(ii,:)=tlat(:);
        else
          isold=0;
          for jj=1:ii
            dx(1)=xxlat(jj,1)-tlat(1);
            dx(2)=xxlat(jj,2)-tlat(2);
            dx(3)=xxlat(jj,3)-tlat(3);
            if (norm(dx)<1e-4)
              isold=1;
            end
          end
          if (isold==0)
            ii=ii+1;
            xxlat(ii,:)=tlat;
          end
        end  % if ii==0
      end  % iz
    end  %iy
  end  %ix

  nxx=size(xx,1);
  nxxx=nxx*nn;
  
  xxx=zeros(nxxx,3);
  for ii=1:nn
    for jj=1:nxx
      xxx((ii-1)*nxx+jj,:)=xx(jj,:)/sc+xxlat(ii,:);
      xxx((ii-1)*nxx+jj,:)=xxx((ii-1)*nxx+jj,:)-floor(xxx((ii-1)*nxx+jj,:));
    end  % jj
  end  % ii

end

