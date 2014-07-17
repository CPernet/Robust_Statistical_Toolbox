function  [MIxy RMIxy1 RMIxy2 RMIxy3]=kernlden2d(x,y,pnts)

%__________________________________________________________________________
% Usage: calculations of two dimensions kernel density and mutual
% information by kernels.
 
% Inputs: 
%   y  and x are  vertical vectors of time series.
%   pnts is number of points that kernels should be evaluated.
 
% Output:
%   MIxy is Mutual information between x and y. RMIxy1 to RMIxy3 are 
%   rescaled form of MIxy.
%   A 3D graph will be drawn by running the code.
 
% Keywords: Bivariate kernel density plot, Mutual Information Estimation,
% Nonlinear Correlation.
 
% Copyright(c) Shapour Mohammadi, University of Tehran, 2009
% shmohammadi@gmail.com
%__________________________________________________________________________

d=2;
n=length(x);

endp1=ceil(pnts/10);
endp2=ceil(pnts/20);

minx=min(x);maxx=max(x);grx=(maxx-minx)/(pnts-endp1);
miny=min(y);maxy=max(y);gry=(maxy-miny)/(pnts-endp1);

h1x=(4/(3*n))^(1/5)*std(x);
h1y=(4/(3*n))^(1/5)*std(y);

for k=1:pnts
xi(k,1)=minx+grx*(k-endp2);
yi(k,1)=miny+gry*(k-endp2);

fx(k,1)=(1/((2*pi)^0.5*n*h1x))*sum(exp(-((xi(k,1)-x).^2)/(2*h1x^2)));
fy(k,1)=(1/((2*pi)^0.5*n*h1y))*sum(exp(-((yi(k,1)-y).^2)/(2*h1y^2)));

px(k,1)=(1/((2*pi)^0.5*n*h1x))*sum(exp(-((xi(k,1)-x).^2)/(2*h1x^2)))*grx;
py(k,1)=(1/((2*pi)^0.5*n*h1y))*sum(exp(-((yi(k,1)-y).^2)/(2*h1y^2)))*gry;

end

[gx gy]=meshgrid(xi,yi);

sigma=((n*var(x)+n*var(y))/(n+n))^0.5;
h=sigma*(4/(d+2))^(1/(d+4))*(n^(-1/(d+4)));
tic
for i=1:pnts
    for j=1:pnts
       
       fxy(i,j)=(1/(2*pi*n*h^2))*sum(exp(-((gx(i,j)-x).^2+...
           (gy(i,j)-y).^2)/(2*h^2)));
       pxy(i,j)=(1/(2*pi*n*h^2))*sum(exp(-((gx(i,j)-x).^2+...
           (gy(i,j)-y).^2)/(2*h^2)))*grx*gry;
       I1xy(i,j)= pxy(i,j)*log(pxy(i,j)/(px(i)*py(j)));
    end
end
toc 

imagesc(fxy); xlabel('data1'); ylabel('data2');
%figure;set(gcf,'Color','w');
%surf(gx,gy,fxy); xlabel('data1'); ylabel('data2'); zlabel('density')

Hx=-(px'*log(px));
Hy=-(py'*log(py));

MIxy=(sum(sum(I1xy)));
RMIxy1=2*MIxy/(Hx+Hy);
RMIxy2=MIxy/(Hx*Hy)^0.5;
RMIxy3=MIxy/min(Hx,Hy);
