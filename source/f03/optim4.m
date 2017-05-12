clear
close

str{1,1}='Maize';str{1,2}='Soy';str{1,3}='Millet';,str{1,4}='Swheat';
str{1,5}='Wwheat';str{1,6}='Sbarley';str{1,7}='Wbarley';str{1,8}='Sorghum';

str{2,1}=10;str{2,2}=5;str{2,3}=5;str{2,4}=6;str{2,5}=10;str{2,6}=8;str{2,7}=8;str{2,8}=6;


ress=200;

dirr{1}='/home/sm1epk/SDGVM/SOUP2/source/f03';
dirr{2}='/data/sm1epk/SDGVM_runs/tempoutputopt/opt';
dirr{3}='/data/sm1epk/SDGVM_runs/tempoutputopt/';

%Reads the parameter file which holds the crop ID,the number of batches and the number 
%of gridcells to be used in the optimization
eval(['cd ',dirr{1}]);
a=dlmread('nopts.dat');
crop=a(1);
nopts=a(2);
npoints=a(3);

mm=[];
for i=1:nopts
  eval(['cd ',dirr{2},num2str(i)]);
  load('opt_par.mat')
  mm=[mm m];
end
mm=mm(:,:);

eval(['cd ',dirr{1}]);
%Probability Density Figure
figure
ecornerplot_me(m,'ks',true,'color',[.5 .5 .5],...
   'names',{'Par1' 'Par2' 'log(\sigma)'},'scatter',false)
eval(['cd ',dirr{3}]);
print(gcf,['pdf_all'],'-dpng',['-r',num2str(ress)])
close(gcf)
  
eval(['cd ',dirr{1}]);
%MAP estimate of parameters
m2=mm;
[count edges mid loc]=histcn(m2',10,10,10);%'
[C,ia,ic]=unique(loc,'rows');
ix=C(mode(ic),:);
sub=mid{1};ma(1)=sub(ix(1));sub=mid{2};ma(2)=sub(ix(2));sub=mid{3};ma(3)=sub(ix(3));

ma_all=ma;

hist3(mm(1:2,:)',[20 20])%'  
view(225,32)
set(gcf,'renderer','opengl')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
eval(['cd ',dirr{3}]);  
print(gcf,['hist_all'],'-dpng',['-r',num2str(ress)])
close(gcf)
 


ref_all=[];
sim_all=[];

for i=1:nopts
  %Simulated values for the batch
  sim=zeros(npoints,2);


  eval(['cd ',dirr{2},num2str(i)]);
  load('opt_par.mat')

  
  %Auto-correlation Figure
  eval(['cd ',dirr{1}]);
  figure 
  [C,lags,ESS]=eacorr_me(m);
  plot(lags,C,'.-',lags([1 end]),[0 0],'k');
  grid on
  xlabel('lags')
  ylabel('autocorrelation');
  text(lags(1),0,sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','left')
  title('Markov Chain Auto Correlation')
  legend('Par1','Par2','\sigma')
  eval(['cd ',dirr{2},num2str(i)]);
  print(gcf,['auto_',num2str(i)],'-dpng',['-r',num2str(ress)])
  close(gcf)

  
  %Probability Density Figure
  eval(['cd ',dirr{1}]);
  figure
  ecornerplot_me(m,'ks',true,'color',[.5 .5 .5],...
     'names',{'Par1' 'Par2' 'log(\sigma)'},'scatter',false)
  eval(['cd ',dirr{2},num2str(i)]);  
  print(gcf,['pdf_',num2str(i)],'-dpng',['-r',num2str(ress)])
  close(gcf)

  
  %Reads the reference yield data for the gridcells in the batch
  eval(['cd ',dirr{2},num2str(i)]);    
  ref=dlmread('yie_dat.dat');
  ref_all=[ref_all;ref];  


  %Gets the MAP for the specific batch
  eval(['cd ',dirr{1}]);
  %MAP estimate of parameters
  m2=m(:,:);
  [count edges mid loc]=histcn(m2',5,5,5);%'
  [C,ia,ic]=unique(loc,'rows');
  ix=C(mode(ic),:);
  sub=mid{1};ma(1)=sub(ix(1));sub=mid{2};ma(2)=sub(ix(2));sub=mid{3};ma(3)=sub(ix(3));

  
  %Does run for each batch with original values
  or=[450.0;0.0];
  eval(['cd ',dirr{2},num2str(i)]);    
  dlmwrite('opt_par.dat',or,'delimiter','\t','precision',5)
  system(['./sdgvm.exe ./test_two_crops_opt_r_',num2str(i),'.dat /data/sm1epk/SDGVM_runs/tempoutputopt/opt',num2str(i),'/res']);        
  eval(['cd ',dirr{2},num2str(i),'/res']);
  a=dlmread([str{1,crop},'ryield.dat']);
  sim(:,1)=(1/0.45)*10/1000*mean(a(:,end-11:end-1));


  %Does run with MAP values obtained from this batch
  or=[ma(1);ma(2)];
  eval(['cd ',dirr{2},num2str(i)]);  
  dlmwrite('opt_par.dat',or,'delimiter','\t','precision',5)
  system(['./sdgvm.exe ./test_two_crops_opt_r_',num2str(i),'.dat /data/sm1epk/SDGVM_runs/tempoutputopt/opt',num2str(i),'/res']);        
  eval(['cd ',dirr{2},num2str(i),'/res']);
  a=dlmread([str{1,crop},'ryield.dat']);
  sim(:,2)=(1/0.45)*10/1000*mean(a(:,end-11:end-1));

  %Does run with MAP values obtained from all batches together
  or=[ma_all(1);ma_all(2)];
  eval(['cd ',dirr{2},num2str(i)]);  
  dlmwrite('opt_par.dat',or,'delimiter','\t','precision',5)
  system(['./sdgvm.exe ./test_two_crops_opt_r_',num2str(i),'.dat /data/sm1epk/SDGVM_runs/tempoutputopt/opt',num2str(i),'/res']);        
  eval(['cd ',dirr{2},num2str(i),'/res']);
  a=dlmread([str{1,crop},'ryield.dat']);
  sim(:,3)=(1/0.45)*10/1000*mean(a(:,end-11:end-1));


  sim_all=[sim_all;sim];


  %Does the scatter plot for each batch with original parameters and MAP for the batch
  figure
  sub=sim(:,1);
  ind=ref>0 & sub>0;
  scatter(ref(ind),sub(ind))
  hold
  sub=sim(:,2);
  ind=ref>0 & sub>0;
  scatter(ref(ind),sub(ind))  
  plot([0 str{2,crop}],[0 str{2,crop}])
  hold off
  eval(['cd ',dirr{2},num2str(i)]);    
  legend('Orig','Optim','Location','SouthEast')
  set(gca,'XLim',[0 str{2,crop}],'YLim',[0 str{2,crop}]);
  print(gcf,['sc_',num2str(i)],'-dpng',['-r',num2str(ress)])
  close(gcf)

end

%Does the scatter plot for all points,pooled vs batch MAP
eval(['cd ',dirr{3}]); 
figure
scatter(ref_all,sim_all(:,2),'filled','SizeData',10)
hold
scatter(ref_all,sim_all(:,3),'filled','SizeData',10)
plot(0:1:20,0:1:20,'black')
set(gca,'XLim',[0 str{2,crop}],'YLim',[0 str{2,crop}])
xlabel('Data Yield tn/hc')
ylabel('SDGVM Yield tn/hc')
legend('Individually','Pooled');
print(gcf,['sc_all'],'-dpng',['-r',num2str(ress)])
close(gcf)

eval(['cd ',dirr{1}]);
