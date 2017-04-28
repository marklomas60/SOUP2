clear

dirr{1}='/home/sm1epk/SDGVM/SOUP2/source/f03';

str{1,1}='Maize';
%Lats/lons of the gridcells to be considered in optimization
str{2,1}={'35.25 -78.25','46.75 -98.25','37.75 -121.25','28.25 -106.25','14.75 -88.25','21.25 -100.75','18.25 -98.75','-34.25 -70.75','-21.75 -47.75','-32.75 -63.75','-35.75 -59.75','41.75 -4.25','39.25 -6.25','46.75 0.25','45.75 4.75','45.25 8.75','40.25 22.25','49.25 36.75','28.75 30.75','54.75 38.25','11.75 6.25','-14.75 27.25','-1.25 37.25','31.25 72.75','17.25 81.25','44.75 86.75','-7.75 111.25','44.75 125.25','34.75 116.75','22.75 109.25','30.25 105.25','34.25 119.75','25.25 113.25'};
str{2,1}={'35.25 -78.25','46.75 -98.25'};
str{5,1}='maize';
str{7,1}='Maize';
str{9,1}=0.8;

str{1,2}='Soy';
%Lats/lons of the gridcells to be considered in optimization
str{2,2}={'35.25 22.45';'21.75 13.25';'31.75 42.25';'1.75 3.25';'27.75 33.25'};
str{5,2}='soybean';
str{7,2}='Soy';
str{9,2}=0.85;

for crop=1:1


  cd /data/sm1epk/SDGVM_runs  
  %Parameter file template
  fin=fopen('test_two_crops_opt.dat','r');
  %Parameter file.Same as above but it will add the lat/lons to be considered by replacins the REPLACE STRING
  fout=fopen('test_two_crops_opt_r.dat','w');
  coun=0;
  %Holds the lat/lon of the gridcells
  sub=cellstr(str{2,crop});
  %Read every line of the template input file till the end of the file
  while ~feof(fin)
    %Read each line as a string
    s=fgetl(fin);
    %If the string REPLACE is in the line,replace with the lats/lons to be considered in optimization
    %The amount of REPLACE in the template file must equal the amount of gridcells you need
    if(strcmp(s,'REPLACE'))
      coun=coun+1;
      s=regexprep(s,'REPLACE',sub(coun));
    end
    %Print the new line
    fprintf(fout,'%s\n',s);
  end
  fclose(fin);
  %Now you have the parameter file with the lats/lons you consider in the optimization
  fclose(fout);


  %We have the lats/lons we will consider as strings in str{2,*}.Here we get them
  %as numbers and pixel image space
  lat=89.75:-0.5:-89.75;
  lon=-179.75:0.5:179.75;
  %Lat/lons Im considering in optimization
  sub=str{2,crop};
  latlon=zeros(size(sub,2),2);
  latloni=zeros(size(sub,2),2);
  for i=1:size(sub,2)
    subb=sub{i};
    subb=sub{i};
    subb=strsplit(subb);
    latlon(i,1)=str2num(subb{1});latlon(i,2)=str2num(subb{2});
    [x,latloni(i,1)]=min(abs(lat-str2num(subb{1})));
    [x,latloni(i,2)]=min(abs(lon-str2num(subb{2})));
  end
  %Lat/lons in numbers
  str{3,crop}=latlon;
  %Lat/lons in pixel space
  str{4,crop}=latloni;


  %Reads the yield data I will calibrate against.I know the pixel space where my gridcells are
  cd /data/sm1epk/crop_sets/SAGE
  load('yield_data.mat');
  eval(['sub=y_data2.',str{5,crop},'_yie;']);
  %Pixel space of lat/lons of the gridcells
  latloni=str{4,crop};
  %Gets the yield of the gridcells I will calibrate against.Multiplies to get dry yield
  str{6,crop}=str{9,crop}*sub(sub2ind(size(sub),latloni(:,1),latloni(:,2)));

  %Compiles
  eval(['cd ',dirr{1}]);
  system('make clean');system('make');
  eval(['cd ',dirr{1},'/bin']);
  system(['cp ./sdgvm.exe ',dirr{1}]);
  eval(['cd ',dirr{1}]);


  or=[450.0;0.0];
  eval(['cd ',dirr{1}]);
  dlmwrite('opt_par.dat',or,'delimiter','\t','precision',5)
  


  
  eval(['cd ',dirr{1}]);
  system(['./sdgvm.exe /data/sm1epk/SDGVM_runs/test_two_crops_opt_r.dat']);

  eval(['cd /data/sm1epk/SDGVM_runs/tempoutputr'])
  a=dlmread([str{7,crop},'ryield.dat']);
  str{8,crop}=(1/0.45)*10/1000*a(:,end-7);
  
  nn=5;
  sigma=std(str{6,crop}-str{8,crop});  
  lognormpdf=@(x,mu,sigma)-0.5*((x-mu)./sigma).^2  -log(sqrt(2*pi).*sigma);


  logLike=@(m)sum(lognormpdf(str{6,crop},sim_y_opt(m,str,crop),exp(m(3))));


  logprior = @(m)(m(1)>50)&(m(1)<1000)&(m(2)>0)&(m(2)<2.);
  mball(1,:)=unifrnd(50,1000,nn,1);
  mball(2,:)=unifrnd(0.,2,nn,1);
  mball(3,:)=log(sigma)*(rand(nn,1))*1;
  
  eval(['cd ',dirr{1}]);
  m=gwmcmc(mball,{logprior logLike},300,'burnin',.2,'stepsize',2);

end

cd /home/sm1epk/SDGVM/SOUP2/source/f03





