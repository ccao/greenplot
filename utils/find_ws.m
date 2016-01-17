function [ rvec, wt ] = find_ws( slat, rmax )
%find_ws returns the rvecs and weights for a specific cell
%   Detailed explanation goes here
  metric=slat*transpose(slat);
  for loop=0:1
    nnrpt=0;
    for ir1=-rmax(1):rmax(1)
      for ir2=-rmax(2):rmax(2)
        for ir3=-rmax(3):rmax(3)
          ir=[ir1, ir2, ir3];
          %
          dist=zeros(1,125);
          ii=0;
          for iir1=-2:2
            for iir2=-2:2
              for iir3=-2:2
                ii=ii+1;
                iir=[iir1, iir2, iir3];
                dx=ir-iir.*transpose(rmax);
                dist(ii)=dx*metric*transpose(dx);
              end
            end
          end
          %
          if ((dist(63)-min(dist))<1e-7)
            nnrpt=nnrpt+1;
            if (loop==1)
              for ii=1:125
                if (abs(dist(ii)-dist(63))<1e-7)
                  wt(nnrpt)=wt(nnrpt)+1;
                end
              end
              rvec(nnrpt,:)=ir;
            end  % loop=1
          end
          %
        end  % ir3
      end % ir2
    end % ir1
    if (loop==0)
      rvec=zeros(nnrpt,3);
      wt=zeros(nnrpt,1);
    end

  end % loop
  
  checksum=sum(1.0./wt);
  if (abs(checksum-rmax(1)*rmax(2)*rmax(3))>1E-8)
    fprintf('!!! FATAL ERROR in finding W-S: %12.9f', abs(checksum-rmax(1)*rmax(2)*rmax(3)));
  end

end

