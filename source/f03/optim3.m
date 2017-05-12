clear

str{1,1}='Maize';str{1,2}='Soy';str{1,3}='Millet';,str{1,4}='Swheat';
str{1,5}='Wwheat';str{1,6}='Sbarley';str{1,7}='Wbarley';str{1,8}='Sorghum';

%Reads the parameter file which holds the crop ID,the number of batches and the number 
%of gridcells to be used in the optimization
a=dlmread('nopts.dat');
crop=a(1);
nopts=a(2);
npoints=a(3);


s=pwd;
%Reads the reference yield data for the gridcells in the batch
ref=dlmread('yie_dat.dat');

%Here it does a single run with some initial values for the parameters.
%It does this to get an estimate of the variance between simulations and data
%which will be used by the MCMC

%Initial values of the parameters and writes them in file to be read by SDGVM
or=[450.0;0.0];
eval(['cd ',s]);
dlmwrite('opt_par.dat',or,'delimiter','\t','precision',5)
%Run SDGVM for this batch
system(['./sdgvm.exe ./test_two_crops_opt_r_',s(42:end),'.dat /data/sm1epk/SDGVM_runs/tempoutputopt/opt',s(42:end),'/res']);

%Goes into the folder that holds the runs for the batch and reads yield for the specific crop  
eval(['cd ',s,'/res']);
a=dlmread([str{1,crop},'ryield.dat']);
sim=(1/0.45)*10/1000*mean(a(:,end-11:end-1));

%Finds the indexes where both simulated and data yields are >0.Only these will
%be used for the optimization
ind=ref>0 & sim>0;

%Calculates sigma
sigma=std(ref(ind)-sim(ind));  



lognormpdf=@(x,mu,sigma)-0.5*((x-mu)./sigma).^2  -log(sqrt(2*pi).*sigma);

logLike=@(m)sum(lognormpdf(sim(ind),sim_y_opt(m,str,crop,s,ind),exp(m(3))));

nn=8;
logprior = @(m)(m(1)>50)&(m(1)<1000)&(m(2)>0)&(m(2)<2.);
mball(1,:)=unifrnd(50,1000,nn,1);
mball(2,:)=unifrnd(0.,2,nn,1);
mball(3,:)=log(sigma)*(rand(nn,1));
  
eval(['cd ',s])
m=gwmcmc(mball,{logprior logLike},1000,'burnin',.2,'stepsize',2);

save('opt_par.mat','m')




