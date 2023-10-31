clear all
close all
clc 
%Problem formulation
%objective =@ (x) flamingfunction(x);
objective =@(x) flamingfunction(x)
nv=2;
min=-4;
max=4;

%NSGA Parameter
Max_iter =100;
np=50;
pc=0.7;
nc=round(np*pc/2);
pm=0.2;
nm=round(np*pm);
mu=0.02;
sigma = 0.1;

%initialization
empty.crowdingdistance=[];
empty.rank=[];
empty.dominationset=[];
empty.dominationcount=[];

population=repmat(empty,np,1);
a= repmat(min,np,nv);
b=repmat(max,np,nv);

x=a+(rand(np,nv)).*(b-a);
z=flamingfunction(x,nv);

%non dominate sorting
[population,F]=nondominatedsorting(population,z);

%crowding distance
[population]=crowdingdistance(z,F,population);

%sort the population
[population, z, x, F]= sortpopulation(population, z,x);

%main loop
for it = 1:Max_iter
    % cross over operation
    popc = repmat(empty, 2*nc, 1);
    
    y1 = zeros(nc,nv);
    y2 = zeros(nc,nv);
    zrcrossy1 = zeros(nc,nv);
    zrcrossy2 = zeros(nc,nv);
    
    for k=1:nc
        i1 = randi(np);
        p1 = x(i1,:);
        ipk = 1;
        while ipk<(ipk+1)
            i2 = randi(np);
            if i2 ~= i1
                break;
            end
            ipk = ipk+1;
        end
        p2 = x(i2,:);
        
        [y1(nc,:), y2(nc,:)] = crossover(p1, p2,nv, min, max)
        zrcrossy1(k,:) = flamingfunction(y1(k,:), nv);
        zrcrossy2(k,:) = flamingfunction(y2(k,:), nv);
    end
    zcross = [zrcrossy1;zrcrossy2];
    ycross = [y1;y2];
    
    %mutation operation
    popm = repmat(empty, nm, 1);
    ym = zeros(nm,nv);
    zmutate = zeros(nm,nv);
    
    for k =1:nm
        i = randi(np);
        p = x(i,:);
        
        ym (k,:) = mutation(p, nv, sigma, min, max);
        
        zmutate(k,:) = flamingfunction(ym(k,:), nv);
    end
    
    %merge solution
    population_merge = [population
                        popc
                        popm];
    z_merge = [z
              zcross
              zmutate];
    x_merge =[x
             ycross
             ym];
   %non dominate sorting
   [population_nds,F_nds]=nondominatedsorting( population_merge,z_merge);

   %crowding distance
   [population_cds]=crowdingdistance(z_merge,F_nds,population_nds);

   %sort the population
   [population_sp, zs, xs, Fs]= sortpopulation(population_cds, z_merge,x_merge);   
   
   %survival selection
   population = population_sp(1:np)
   z=zs(1:np,:);
   x = xs(1:np,:);
   
   %non dominate sorting
   [population,F]=nondominatedsorting(population,z);

   %crowding distance
   [population]=crowdingdistance(z,F,population);

   %sort the population
   [population, z, x, F]= sortpopulation(population, z,x);

   %plot 
   figure(1);
   plot(z(:,1),z(:,2), 'o',col='r')
   
   xlabel('1st Obj')
   ylabel('2nd Obj')
   title('NSGA II')
   grid on
end
        
