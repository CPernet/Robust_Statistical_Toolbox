function dhat=akerd(x,pts,plotit)
%dhat=akerd(x,pts,plotit)
%
% Computes an adaptive kernel density estimate for univariate data,
% using expected frequency as initial estimate of density (See Silverman, 1986).
%
% Matlab adaptation of R functions from Rand Wilcox:
% akerd, rdplot & near
% see Wilcox, R. R. (2012). Introduction to robust estimation and
% hypothesis testing (3rd ed.). Amsterdam, Boston: Academic Press.
%
% See also MADR, IDEALF, WINVAR, DISKER

% -----------------------------------------------------------
% Guillaume Rousselet, University of Glasgow, 16 October 2012
% -----------------------------------------------------------

if nargin<3 || isempty(plotit)
    plotit=0;
end
if nargin<2 
    pts=[];
end

fr=.8;aval=.5;
x=x(~isnan(x));
x=x(:);
x=sort(x);
n=length(x);
m=madr(x);
if m==0
    [ql,qu]=idealf(x);
    m=(qu-ql)./1.34898;%(qnorm(.75)-qnorm(.25))
end
if m==0
    m=sqrt(winvar(x,20)/.4129);
end
if m==0
    error('in akerd: all measures of dispersion are equal to 0')
end

% rdplot: Expected frequency curve --------------------------------------
rmd=zeros(length(x),1);
for p=1:length(x) % near: determine which values in x are near pt based on fr * mad
    rmd(p)=sum(abs(x-x(p))<=(fr.*m));
end
if m~=0
    fhat=(rmd./(2.*fr.*m))./n;
end
% -----------------------------------------------------------------------

if m>0
    fhat=fhat./(2.*fr.*m);
end

% compute hval
sig=sqrt(var(x));
[ql,qu]=idealf(x);
iq=(qu-ql)/1.34;
A=min([sig iq]);
if A==0
    A=sqrt(winvar(x,20))/.64;
end
hval=1.06.*A./length(x)^(.2);
% See Silverman, 1986, pp. 47-48

gm=exp(mean(log(fhat(fhat>0)))); %fprintf('%.7f',gm)
alam=(fhat./gm).^(0-aval);

if isempty(pts);pts=sort(x);end
dhat=zeros(length(pts),1);
pts=sort(pts);
for p=1:length(pts)
    temp=(pts(p)-x)./(hval.*alam);
    epan=zeros(n,1);
    epan(abs(temp)<sqrt(5))=.75*(1-.2*temp(abs(temp)<sqrt(5)).^2)./sqrt(5);
    dhat(p)=mean(epan./(alam.*hval));
end

if plotit==1
    figure('Color','w','Name','Adaptive kernel density estimate','NumberTitle','off')
    hold on
    plot(pts,dhat,'ko')
    plot(pts,dhat,'k-')
end
