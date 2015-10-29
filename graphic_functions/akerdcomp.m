function akerdcomp(a,b)

% GAR, Glasgow University, 26 Nov 2012

da=akerd(a,a,0);
db=akerd(b,b,0);
figure('Color','w','NumberTitle','off');hold on
plot(sort(a),da,'k','LineWidth',2)
plot(sort(b),db,'k--','LineWidth',2)
box on
set(gca,'FontSize',14)
