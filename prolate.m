%
%  Prolate filter design
%
%  function v=prolate(N,sigma);
%  N: filter order, should be even
%  sigma: width of the main lobe of the window response
%
function v=prolate(N,sigma);
M=N/2;
R=zeros(M,M);
R(1,1)=0.5-sigma;
k=1
l=[1:M];
R(k,l+1)=-sqrt(2)*sigma*sinc(2*sigma*l);
l=1;
k=[1:M];
R(k+1,l)=-sqrt(2)*sigma*sinc(2*sigma*k);
for k=1:M
    R(k+1,k+1)=0.5-sigma*(1+sinc(4*sigma*k));
end
for k=1:M
    for l=1:M
        if k~=l
            kpl=k+l;
            kml=k-l;
            R(k+1,l+1)=-sigma*sinc(2*sigma*kpl)-sigma*sinc(2*sigma*kml);
        end
    end
end
[V,D]=eig(R);
[a,indx]=min(diag(D));
v=V(:,indx);
v=[flipud(v(2:end)); v];
v=sign(v(1+M))*v;
v(1+M)=v(1+M)*sqrt(2);
v=v/v(1+M);

