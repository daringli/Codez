data=dlmread("gamma_addback.tsv",'\t',1,0);

data_gamma_x=data(:,2)(find(data(:,3)<=10));
data_gamma_y=data(:,3)(find(data(:,3)<=10));
data_gamma_events=data(:,1)(find(data(:,3)<=10));

%plot 2D histogram of gamma mult at different E* and E_gamma
[counts,xbin,ybin] = hist2d([data_gamma_x,data_gamma_y],30);
figure(1)
clf
 colormap ("hot");
imagesc(xbin,ybin,counts)
colorbar()
set(gca,'YDir','normal')

%gather event-wise number of identified gamma
n_events=data_gamma_events(end);
Eex=zeros(1,n_events);
n_gamma_id=zeros(1,n_events);
for i=1:n_events
indices=find(data(:,1)==i)
Eex(i)=data(indices(1),2)
n_gamma_id(i)=sum(data(indices',3)<=10);
endfor

%gather event-wise number of fired gamma
n_gamma=zeros(1,n_events);
event_nr=0;
fid=fopen("test.gun",'r');
while (line=fgetl(fid))!=-1
  if (line(1)=='*')
    ++event_nr;
  endif
  if(length(line)>=3)
    if(line(3)=='g')
      n_gamma(event_nr)=n_gamma(event_nr)+1;
    endif
  endif
endwhile

%next we average over datapoints with the same energy
n_gamma_avg=zeros(1,length(unique(Eex)));
n_gamma_id_avg=zeros(1,length(unique(Eex)));

i=1;
for E=unique(Eex)
  indices=find(Eex==E);
  n_gamma_avg(i)=sum(n_gamma(indices));
  n_gamma_id_avg(i)=sum(n_gamma_id(indices));
  i++;
endfor

%take ratio
firedOverIDed=n_gamma_avg./n_gamma_id_avg;
figure(2)
plot(unique(Eex),firedOverIDed);


%{
%get gamma mult for each event
%id:ed
n_events=data_gamma_events(end);
n_gammas_ided=zeros(1,n_events);
Eex=zeros(1,n_events);

for i=1:n_events
  n_gammas_ided(i)=length(find(data_gamma_events==i));
  Eex(i)=data(find(data(:,1)==i)(1),2);
endfor


figure(2)
clf
colormap ("hot");
hold on
[counts,x,y]=hist2d([Eex',n_gammas_ided'],30);
imagesc(x,y,counts);
set(gca,'YDir','normal')

%get the actual gamma mult by parsing the gunfile
n_gammas=zeros(1,length(unique(data(:,1))));
event_nr=0;
fid=fopen("test.gun",'r');
while (line=fgetl(fid))!=-1
  if (line(1)=='*')
    ++event_nr;
  endif
  if(length(line)>=3)
    if(line(3)=='g')
      n_gammas(event_nr)=n_gammas(event_nr)+1;
    endif
  endif

endwhile
figure(3)
clf
colormap ("hot");
hold on
[counts,x,y]=hist2d([Eex',n_gammas'],30);
imagesc(x,y,counts);
set(gca,'YDir','normal')
%}