function f = sim_y_opt(x,str,crop,s,ind)
  
  or=x(1:end-1);
  dlmwrite('opt_par.dat',or,'delimiter','\t','precision',5)

  system(['./sdgvm.exe ./test_two_crops_opt_r_',s(42:end),'.dat /data/sm1epk/SDGVM_runs/tempoutputopt/opt',s(42:end),'/res']);
    
  eval(['cd ',s,'/res']);

  a=dlmread([str{1,crop},'ryield.dat']);

  eval(['cd ',s]);  

  f=(1/0.45)*10/1000*mean(a(:,end-11:end-1));
  f=f(ind);
    
end
