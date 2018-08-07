

%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%
clear;

%Field Dimensions - x and y maximum (in meters)
xm=300;
ym=300;

%x and y Coordinates of the Sink
sink.x=0.5*xm;
sink.y=0.5*ym;

%Number of Nodes in the field
n=100

%Optimal Election Probability of a node to become cluster head
p=0.1;

%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.5;

%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;

%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;

%Data Aggregation Energy
EDA=5*0.000000001;

%Values for Hetereogeneity
%Percentage of nodes than are advanced
m=0.2;

%\alpha
a=3;

%maximum number of rounds
rmax=5000;

%Computation of do
do=sqrt(Efs/Emp);

%%%%%%%%%%%%%%%%%%%%%%%%% LEACH %%%%%%%%%%%%%%%%%%%%%%%%
%Creation of the random Sensor Network
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    S(i).yd=rand(1,1)*ym;
    S(i).G=0;
    S(i).type='N';
    S(i).E=Eo;  
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
    
        
%First Iteration
figure(1);

%?????
cluster=1;

first_dead_l=0;
half_dead_l=0;
all_dead_l=0;

tic
for r=0:1:rmax
    r

 %Operation for epoch
  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        S(i).G=0;
    end
  end
    
	%Number of dead nodes
	dead=0;

	for i=1:1:n
	    %checking if there is a dead node
	    if (S(i).E<=0)
	       dead=dead+1;
	    else
	        S(i).type='N';
	    end
    end
    
    if(first_dead_l==0)
        if(dead~=0)
            first_dead_l=r;
        end
    end

    if(half_dead_l==0)
        if(dead>0.5*n)
            half_dead_l=r;
            for i=1:1:n	    
                if (S(i).E<=0)
                    plot(S(i).xd,S(i).yd,'red x');
                else                        
                    plot(S(i).xd,S(i).yd,'blue o');
                end
            end
            plot(S(n+1).xd, S(n+1).yd, 'black +');
        end
    end
    
    if(all_dead_l==0)
        if(dead==n)
            all_dead_l=r;
        end
    end
    
	T(r+1)=toc;
	ns=n-dead;

	array_X4(r+1)= r+1;
	array_Y4(r+1)= ns;


    cluster=1;
    for i=1:1:n
       if(S(i).E>0)
           temp_rand=rand;     
               if ( (S(i).G)<=0)
                   %Election of Cluster Heads
                   if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
                       S(i).type='C';
		               S(i).G=round(1/p)-1;
		               C(cluster).xd=S(i).xd;
		               C(cluster).yd=S(i).yd;
                       
                       distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
                       C(cluster).distance=distance;
                       C(cluster).id=i;
                       cluster=cluster+1;
            
                       %Calculation of Energy dissipated
			            if (distance>do)
			                S(i).E=S(i).E- ( (ETX)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
			            end
			            if (distance<=do)
			                S(i).E=S(i).E- ( (ETX)*(4000)  + Efs*4000*( distance * distance )); 
			            end
                   end     
    
              end
        end 
    end


	%Election of Associated Cluster Head for Normal Nodes
	for i=1:1:n
	   if ( S(i).type=='N' && S(i).E>0 )
	     if(cluster-1>=1)
	       min_dis=sqrt( (S(i).xd- C(1).xd)^2 + (S(i).yd- C(1).yd)^2 );
	       min_dis_cluster=1;
	       for c=2:1:cluster-1
	           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
	           if ( temp<min_dis )
	               min_dis=temp;
	               min_dis_cluster=c;
	           end
	       end
       
	       %Energy dissipated by associated Cluster Head
	        
	            if (min_dis>do)
	                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
	            end
	            if (min_dis<=do)
	                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
	            end
	        %Energy dissipated
	        if(min_dis>0)
	          S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
	      
	        end

           
	     end
	   end
	end
    hold on;

end


%%%%%%%%%%%%%%%%%%%%%%%%% standard SEP %%%%%%%%%%%%%%%%%%%%%%%%
%Creation of the random Sensor Network
figure(2);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    S(i).yd=rand(1,1)*ym;
    S(i).CH=0;
    S(i).type='N';
    %Random Election of Normal Nodes
    if (i>=m*n+1)
       S(i).E=Eo;
	   S(i).G=0;
	   S(i).ENERGY='Nor';
       Sinitial(i)=Eo;
       S(i).p=p/(1+m*a);  
	else
       S(i).E=Eo*(1+a);
       S(i).p=p/(1+m*a)*(1+a);
	   S(i).G=0;
       Sinitial(i)=Eo*(1+a);
	   S(i).ENERGY='Adv';
	end
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
    
    
Maxdisb = 0;
for i=1:1:n
    disb=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
	if(disb>Maxdisb)
		Maxdisb=disb;
	end
    S(i).disb=disb;
    if(S(i).disb <= do)
        S(i).range='internal';
    else
        S(i).range='external';
    end
end
%First Iteration
figure(2);

%?????
cluster=1;

first_dead_s=0;
half_dead_s=0;
all_dead_s=0;

tic
for r=0:1:rmax
    r

 %Operation for epoch
  if(mod(r, round((1+a*m)/p) )==0)
    for i=1:1:n
		if(S(i).ENERGY=='Nor')
        	S(i).G=0;
		end
    end
  end

  if(mod(r, round((1+a*m)/(p*(1+a))) )==0)
    for i=1:1:n
		if(S(i).ENERGY=='Adv')
        	S(i).G=0;
		end
    end
  end

	%Number of dead nodes
	dead=0;

	TOTAL_E=0;
	for i=1:1:n
	    %checking if there is a dead node
	    if (S(i).E<=0)
	       dead=dead+1;
	    else
	        S(i).type='N';
			TOTAL_E=TOTAL_E+S(i).E;
	    end
    end
    
    if(first_dead_s==0)
        if(dead~=0)
            first_dead_s=r;
        end
    end

    if(half_dead_s==0)
        if(dead>0.5*n)
            half_dead_s=r;
            for i=1:1:n	    
                if (S(i).E<=0)
                    plot(S(i).xd,S(i).yd,'red x');
                else  
                    if(S(i).ENERGY=='Adv')
                        plot(S(i).xd,S(i).yd,'blue *');
                    else
                        plot(S(i).xd,S(i).yd,'blue o');
                    end
                end
            end
            plot(S(n+1).xd, S(n+1).yd, 'black +');
        end
    end
    
    if(all_dead_s==0)
        if(dead==n)
            all_dead_s=r;
        end
    end

	T(r+1)=toc;
	ns=n-dead;
	AVG_E=TOTAL_E/ns;
	array_X1(r+1)= r+1;
	array_Y1(r+1)= ns;


    cluster=1;
    for i=1:1:n
       if(S(i).E>0)
           temp_rand=rand;     
               if ( (S(i).G)<=0)
                   %Election of Cluster Heads
                   if(temp_rand<= ((S(i).p/(1-S(i).p*mod(r,round(1/S(i).p))))))
                       S(i).type='C';
		               S(i).G=round(1/p)-1;
                       S(i).CH=S(i).CH+1;
		               C(cluster).xd=S(i).xd;
		               C(cluster).yd=S(i).yd;
                       
                       distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
                       C(cluster).distance=distance;
                       C(cluster).id=i;
                       cluster=cluster+1;
            
                       %Calculation of Energy dissipated
			            if (distance>do)
			                S(i).E=S(i).E- ( (ETX)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
			            end
			            if (distance<=do)
			                S(i).E=S(i).E- ( (ETX)*(4000)  + Efs*4000*( distance * distance )); 
			            end
                   end     
    
              end
        end 
    end


	%Election of Associated Cluster Head for Normal Nodes
	for i=1:1:n
	   if ( S(i).type=='N' && S(i).E>0 )
	     if(cluster-1>=1)
	       min_dis=sqrt( (S(i).xd- C(1).xd)^2 + (S(i).yd- C(1).yd)^2 );
	       min_dis_cluster=1;
	       for c=2:1:cluster-1
	           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
	           if ( temp<min_dis )
	               min_dis=temp;
	               min_dis_cluster=c;
	           end
	       end
       
	       %Energy dissipated by associated Cluster Head
	        
	            if (min_dis>do)
	                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
	            end
	            if (min_dis<=do)
	                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
	            end
	        %Energy dissipated
	        if(min_dis>0)
	          S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
	      
	        end

	     end
	   end
	end
    hold on;

end


%%%%%%%%%%%%%%%%%%%%% multi-hop SEP  %%%%%%%%%%%%%%%%%%%%%%%%%%
%Creation of the random Sensor Network
figure(3);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    S(i).yd=rand(1,1)*ym;
    S(i).G=0;
    S(i).type='N';
	S(i).CH=0;
    %Random Election of Normal Nodes
    if (i>=m*n+1)
       S(i).E=Eo;
       S(i).p=p/(1+m*a);  
	   S(i).ENERGY='Nor';
	   Sinitial(i)=Eo;
	else
       S(i).E=Eo*(1+a);
       S(i).p=p/(1+m*a)*(1+a);
	   S(i).ENERGY='Adv';
	   Sinitial(i)=Eo*(1+a);
	end
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
    
for i=1:1:n
    disb=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
    S(i).disb=disb;
    if(S(i).disb <= do)
        S(i).range='internal';
    else
        S(i).range='external';
    end
end
%First Iteration
figure(3);

%?????
cluster=1;

first_dead_m=0;
half_dead_m=0;
all_dead_m=0;
tic
for r=0:1:rmax
    r

 %Operation for epoch
	if(mod(r, round((1+a*m)/p) )==0)
	    for i=1:1:n
			if(S(i).ENERGY=='Nor')
	        	S(i).G=0;
			end
	    end
	  end

	  if(mod(r, round((1+a*m)/(p*(1+a))) )==0)
	    for i=1:1:n
			if(S(i).ENERGY=='Adv')
	        	S(i).G=0;
			end
	    end
	  end


	%Number of dead nodes
	dead=0;
	TOTAL_E=0;
	for i=1:1:n
	    %checking if there is a dead node
	    if (S(i).E<=0)
	       dead=dead+1;
	    else
	        S(i).type='N';
			TOTAL_E=TOTAL_E+S(i).E;
	    end
    end

    if(first_dead_m==0)
        if(dead~=0)
            first_dead_m=r;
        end
    end

    if(half_dead_m==0)
        if(dead>0.5*n)
            half_dead_m=r;
            for i=1:1:n	    
                if (S(i).E<=0)
                    plot(S(i).xd,S(i).yd,'red x');
                else  
                    if(S(i).ENERGY=='Adv')
                        plot(S(i).xd,S(i).yd,'blue *');
                    else
                        plot(S(i).xd,S(i).yd,'blue o');
                    end
                end
            end
            plot(S(n+1).xd, S(n+1).yd, 'black +');
        end
    end
    
    if(all_dead_m==0)
        if(dead==n)
            all_dead_m=r;
        end
    end
	T(r+1)=toc;
	ns=n-dead;
	AVG_E=TOTAL_E/ns;

	array_X3(r+1)= r+1;
	array_Y3(r+1)= ns;


    cluster=0;
    for i=1:1:n
       if(S(i).E>0)
           temp_rand=rand;     
               if ( (S(i).G)<=0)
                   %Election of Cluster Heads
                   if(temp_rand<= ((S(i).p/(1-S(i).p*mod(r,round(1/S(i).p))))*((S(i).E/Sinitial(i))+(1-S(i).E/Sinitial(i))*(S(i).p/(1+S(i).CH)))))
                       S(i).type='C';
		               S(i).G=round(1/p)-1;
					   S(i).CH=S(i).CH+1;
                       cluster=cluster+1;
		               C(cluster).xd=S(i).xd;
		               C(cluster).yd=S(i).yd;
                       distance = S(i).disb;
                       C(cluster).disb=S(i).disb;
                       C(cluster).id=i;
					   C(cluster).range = S(i).range;
            
                       %Calculation of Energy dissipated
			           % if (distance>do)
			            %    S(i).E=S(i).E- ( (ETX)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
			            %end
			            if (distance<=do)
			                S(i).E=S(i).E- ( (ETX)*(4000)  + Efs*4000*( distance * distance )); 
			            end
                   end     
    
              end
        end 
    end

    clusterin = 0;
	clusterout = 0;
    if(cluster > 0)
      for i=1:1:cluster
	     if(C(i).range=='internal')
		    clusterin = clusterin + 1;
		    CIN(clusterin).xd = C(i).xd;
			CIN(clusterin).yd = C(i).yd;
			CIN(clusterin).disb = C(i).disb;
			CIN(clusterin).id = C(i).id;
	     else
		    clusterout = clusterout + 1;
		    COUT(clusterout).xd = C(i).xd;
			COUT(clusterout).yd = C(i).yd;
			COUT(clusterout).disb = C(i).disb;
			COUT(clusterout).id = C(i).id;
		 end
      end
    end

	if(clusterout==0)
		continue;
	end

	if(clusterin==0)
		for i=1:1:clusterout
			S(COUT(i).id).E=S(COUT(i).id).E- ( (ETX)*(4000) + Emp* 4000 * (  S(COUT(i).id).disb *S(COUT(i).id).disb * S(COUT(i).id).disb * S(COUT(i).id).disb)); 
        end
	else
	%????
	    for cout=1:1:clusterout
		    min_out_in=COUT(cout).disb;
			min_out_in_cluster = 0;
            for cin=1:1:clusterin
				    temp = min(min_out_in,sqrt( (COUT(cout).xd-CIN(cin).xd)^2+ (COUT(cout).xd-CIN(cin).xd)^2));
					if(temp<min_out_in)
					    min_out_in=temp;
						min_out_in_cluster=cin;
					end
		    end
			if(min_out_in_cluster==0)
				S(COUT(cout).id).E=S(COUT(cout).id).E- ( (ETX)*(4000) + Emp* 4000 * (  S(COUT(cout).id).disb *S(COUT(cout).id).disb * S(COUT(cout).id).disb * S(COUT(cout).id).disb)); 
			else
                S(CIN(min_out_in_cluster).id).E = S(CIN(min_out_in_cluster).id).E- ( (ERX + EDA)*4000 );
				if(min_out_in>do)
					S(COUT(cout).id).E=S(COUT(cout).id).E- ( (ETX)*(4000) + Emp* 4000 * min_out_in* min_out_in*min_out_in*min_out_in); 
				else
					S(COUT(cout).id).E=S(COUT(cout).id).E- ( (ETX)*(4000) + Emp* 4000 * min_out_in * min_out_in);
				end
			end
		end
	end


	%Election of Associated Cluster Head for Normal Nodes
	for i=1:1:n
	   if ( S(i).type=='N' && S(i).E>0 )
	     if(S(i).range=='internal')
		    S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( S(i).disb * S(i).disb)); 
	     else
	   		 %???????????????ch
		     min_dis=sqrt( (S(i).xd-COUT(1).xd)^2 + (S(i).yd-COUT(1).yd)^2 ) ;
   	         min_dis_clusterout=1;
			 if(clusterout>1)
			 %???
			    for c=2:1:clusterout
	   	           temp=min(min_dis,sqrt( (S(i).xd-COUT(c).xd)^2 + (S(i).yd-COUT(c).yd)^2 ) );
	   	           if ( temp<min_dis )
	   	               min_dis=temp;
	   	               min_dis_clusterout=c;
	   	           end
				end
			 else
			 %???
                 if (min_dis>do)
                   S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
	             end
	             if (min_dis<=do)
	                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
	             end
		 
                 S(COUT(min_dis_clusterout).id).E = S(COUT(min_dis_clusterout).id).E- ( (ERX + EDA)*4000 );
		     end
		  end
		 end
	end
		 
		 
       
	   
	   
	   
    hold on;

end



%%plot(array_X1, array_Y1, 'blue -',array_X1, array_Y3,'red -',array_X1, array_Y4,'black -')