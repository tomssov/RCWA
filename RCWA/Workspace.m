
nodes = loc.Nt+1.8;
hold off
data = max(squeeze(abs(El(1,1,3,:,mat_cat==1))));

plot((nmz(mat_cat==1)-150), data,'r')
hold on


for(i=1:10:30)
plot((nmz(mat_cat==1)-150),squeeze(abs(El(1,1,3,i,mat_cat==1))));
end

[val,locs]=findpeaks(data)



scatter(nmz(locs)-150, val)
hold off

diff=locs(2:end)-locs(1:end-1);
numb=12;
width=180/(loc.Ng);
%width=1/(15*loc.Nt^2);


sample = data(locs(numb)-ceil(diff(1)/4):locs(numb)+ceil(diff(1)/4));
m=(sample(1)-sample(end))/length(sample);
x=linspace(-length(sample)/2,length(sample)/2,length(sample));
sample2= sample-(-x*m + sample(1));
maxp=max(sample2);
sample2=sample2/maxp;
plot(x,sample2,x,exp(-abs(x)/15),x,exp(-(x).^2/(10^2)),x,exp(-sqrt(abs(x).^3)/40),x,1./(width*x.^2+1))