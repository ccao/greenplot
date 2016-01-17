function [ revlist ] = revindex( rvec )
%revindex generates a reverse list of rvecs
%   Detailed explanation goes here
  nrpt=size(rvec, 1);
  rmax=max(rvec);
  nn=2*rmax+1;
  revlist=zeros(nn(1),nn(2),nn(3));
  for ii=1:nrpt
    revlist(rvec(ii,1)+rmax(1)+1, rvec(ii,2)+rmax(2)+1, rvec(ii,3)+rmax(3)+1)=ii;
  end
%
end

