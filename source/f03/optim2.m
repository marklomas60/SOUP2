clear

dirr{1}='/home/sm1epk/SDGVM/SOUP2/source/f03';

str{1,1}=29;str{2,1}='maize';  str{3,1}=0.80;str{4,1}=9;
str{1,2}=30;str{2,2}='soybean';str{3,2}=0.85;str{4,2}=4;
str{1,3}=31;str{2,3}='millet'; str{3,3}=0.85;str{4,3}=4;
str{1,4}=32;str{2,4}='wheat';  str{3,4}=0.85;str{4,4}=5;
str{1,5}=33;str{2,5}='wheat';  str{3,5}=0.85;str{4,5}=9;
str{1,6}=34;str{2,6}='barley'; str{3,6}=0.85;str{4,6}=7;
str{1,7}=35;str{2,7}='barley'; str{3,7}=0.85;str{4,7}=7;
str{1,8}=36;str{2,8}='sorghum';str{3,8}=0.85;str{4,8}=5;


%Reads the parameter file which holds the crop ID,the number of batches and the number 
%of gridcells to be used in the optimization
eval(['cd ',dirr{1}]);
a=dlmread('nopts.dat');
crop=a(1);
nopts=a(2);
npoints=a(3);

%Reads the cover I will use to run SDGVM for each crop.We need this because we want to
%make sure that the gridcells we will pick to optimize will have crop cover in SDGVM
cd /home/sm1epk/SDGVM/SDGVM_misc/SDGVM/data/land_use/global/hyde_crops
a=dlmread(['cont_lu-',num2str(str{1,crop}),'-2001.dat']);
a=reshape(a,[],360);
a=a';%'
a(a>100)=0;

%Reads the yield data I will calibrate against.Also multiply to get dry yield
cd /data/sm1epk/crop_sets/SAGE
load('yield_data.mat');
eval(['sub=y_data2.',str{2,crop},'_yie;']);
sub=str{3,crop}*sub;

lon=-179.75:0.5:179.75;
lat=89.75:-0.5:-89.75;

%Creates a number of intervals same as the number of gridcells I will optimize
%of yield from 0 to a high observed value.I will pick one gridcell from each interval
inter=linspace(0,str{4,crop},npoints);
inter=[inter str{4,crop}+1];

ind=zeros(nopts,npoints,2);
yie_dat=zeros(nopts,npoints);

%For each point
for i=1:size(inter,2)-1
    %Find the gridcells in the dataset that fall in the interval
    k=find(sub>inter(i) & sub<inter(i+1));	
    %For each of the batches
    for j=1:nopts
        cov=0;
        %Keep sampling random gridcells until the cover is positive
        while(cov==0)
            k0=datasample(k,1);
            cov=a(k0);
        end
        %Make the index in image space
        [k1,k2]=ind2sub(size(sub),k0);
        ind(j,i,1)=k1;ind(j,i,2)=k2;
        %Store as lat and lon to be written
        ind(j,i,1)=lat(k1);ind(j,i,2)=lon(k2);
        %Get the yield values of the dataset
        yie_dat(j,i)=sub(k1,k2);
    end
end

%For each batch write in the relevant batch folder the latslons used and the reference yield data
for i=1:nopts
    eval(['cd /data/sm1epk/SDGVM_runs/tempoutputopt/opt',num2str(i)]);
    dlmwrite('latslons.dat',squeeze(ind(i,:,:)),'delimiter','\t','precision',5)
    dlmwrite('yie_dat.dat',squeeze(yie_dat(i,:))','delimiter','\t','precision',5)%'
end


%Compiles
eval(['cd ',dirr{1}]);
system('module load compilers/gcc/5.2')
system('make clean');system('make');
eval(['cd ',dirr{1},'/bin']);
system(['cp ./sdgvm.exe ',dirr{1}]);
eval(['cd ',dirr{1}]);


%For each batch
for i=1:nopts
    cd /data/sm1epk/SDGVM_runs  
    %Parameter file template
    fin=fopen('test_two_crops_opt.dat','r');
    %Parameter file.Same as above but it will add the lat/lons to be considered by replacing the REPLACE STRING
    eval(['cd /data/sm1epk/SDGVM_runs/tempoutputopt/opt',num2str(i)]);
    fout=fopen(['test_two_crops_opt_r_',num2str(i),'.dat'],'w');

    coun=0;
    %Read every line of the template input file till the end of the file
    while ~feof(fin)
        %Read each line as a string
        s=fgetl(fin);
        %If the string REPLACE is in the line,replace with the lats/lons to be considered in optimization
        %The amount of REPLACE in the template file must equal the amount of gridcells you need
        if(strcmp(s,'REPLACE'))
            coun=coun+1;
            s=regexprep(s,'REPLACE',[num2str(squeeze(ind(i,coun,1))),' ', num2str(squeeze(ind(i,coun,2)))]);
        end
        %Print the new line
        fprintf(fout,'%s\n',s);
    end
    fclose(fin);
    %Now you have the parameter file with the lats/lons you consider in the optimization
    fclose(fout);

    %Copy in each batch folder the necessary files
    eval(['cd ',dirr{1},'/bin']);
    system(['cp ./sdgvm.exe ','/data/sm1epk/SDGVM_runs/tempoutputopt/opt',num2str(i)]);

    eval(['cd ',dirr{1}]);
    system(['cp ./misc_params.dat ','/data/sm1epk/SDGVM_runs/tempoutputopt/opt',num2str(i)]);
  
    eval(['cd ',dirr{1}]);
    system(['cp -R ./inc ','/data/sm1epk/SDGVM_runs/tempoutputopt/opt',num2str(i)]);

    eval(['cd ',dirr{1}]);
    system(['cp ./optim3.m ','/data/sm1epk/SDGVM_runs/tempoutputopt/opt',num2str(i)]);
  
    eval(['cd ',dirr{1}]);
    system(['cp ./nopts.dat ','/data/sm1epk/SDGVM_runs/tempoutputopt/opt',num2str(i)]);

    eval(['cd ',dirr{1}]);
    system(['cp ./sim_y_opt.m ','/data/sm1epk/SDGVM_runs/tempoutputopt/opt',num2str(i)]);
  
    eval(['cd ',dirr{1}]);
    system(['cp ./gwmcmc.m ','/data/sm1epk/SDGVM_runs/tempoutputopt/opt',num2str(i)]);

end


cd /home/sm1epk/SDGVM/SOUP2/source/f03
