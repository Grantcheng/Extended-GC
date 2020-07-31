clc
clear
epstot=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16]
Leng_IE=length(epstot);
Case_Type_Choose=input('Case_Type_Choose is : 0 - 280282290;  1 - 280290371;  2 - 280290371_2;  3 - 280290371_Q; 4 - 280282290_Q; 5 - 280290371_2_Q; ');
Taken_dim=zeros(Leng_IE,1);
for ie=1:Leng_IE
       
    C=epstot(ie);
    Couple_ie=C;
    display(['The Strength is : ' num2str(C)]);
        
   if ie==1
       if Case_Type_Choose==0
       eigvalue=load('HH_Three_280/Couple=0/eig_21_10000000.txt');
       elseif Case_Type_Choose==1
       eigvalue=load('HH_Three_280290371/Couple=0/eig_21_10000000.txt');   
       elseif Case_Type_Choose==2
       eigvalue=load('HH_Three_280290371_2/Couple=0/eig_21_10000000.txt');   
       elseif Case_Type_Choose==3
       eigvalue=load('HH_Three_280290371_Q/Couple=0/eig_21_10000000.txt');  
       elseif Case_Type_Choose==4
       eigvalue=load('NeuroType_Case/Case_1/Couple=0/eig_21_10000000.txt');
       elseif Case_Type_Choose==5
       eigvalue=load('NeuroType_Case/Case_2/Couple=0/eig_21_10000000.txt');    
       end
   elseif ie==2
       if Case_Type_Choose==0
       eigvalue=load('HH_Three_280/Couple=0.02/eig_21_10000000.txt');
       elseif Case_Type_Choose==1
       eigvalue=load('HH_Three_280290371/Couple=0.02/eig_21_10000000.txt');   
       elseif Case_Type_Choose==2
       eigvalue=load('HH_Three_280290371_2/Couple=0.02/eig_21_10000000.txt');   
       elseif Case_Type_Choose==3
       eigvalue=load('HH_Three_280290371_Q/Couple=0.02/eig_21_10000000.txt'); 
       elseif Case_Type_Choose==4
       eigvalue=load('NeuroType_Case/Case_1/Couple=0.02/eig_21_10000000.txt');
       elseif Case_Type_Choose==5
       eigvalue=load('NeuroType_Case/Case_2/Couple=0.02/eig_21_10000000.txt');
       end
    elseif ie==3
       if Case_Type_Choose==0
       eigvalue=load('HH_Three_280/Couple=0.04/eig_21_10000000.txt');
       elseif Case_Type_Choose==1
       eigvalue=load('HH_Three_280290371/Couple=0.04/eig_21_10000000.txt');   
       elseif Case_Type_Choose==2
       eigvalue=load('HH_Three_280290371_2/Couple=0.04/eig_21_10000000.txt');   
       elseif Case_Type_Choose==3
       eigvalue=load('HH_Three_280290371_Q/Couple=0.04/eig_21_10000000.txt'); 
       elseif Case_Type_Choose==4
       eigvalue=load('NeuroType_Case/Case_1/Couple=0.04/eig_21_10000000.txt');
       elseif Case_Type_Choose==5
       eigvalue=load('NeuroType_Case/Case_2/Couple=0.04/eig_21_10000000.txt');
       end
    
    elseif ie==4
       if Case_Type_Choose==0
       eigvalue=load('HH_Three_280/Couple=0.06/eig_21_10000000.txt');
       elseif Case_Type_Choose==1
       eigvalue=load('HH_Three_280290371/Couple=0.06/eig_21_10000000.txt');   
       elseif Case_Type_Choose==2
       eigvalue=load('HH_Three_280290371_2/Couple=0.06/eig_21_10000000.txt');   
       elseif Case_Type_Choose==3
       eigvalue=load('HH_Three_280290371_Q/Couple=0.06/eig_21_10000000.txt'); 
       elseif Case_Type_Choose==4
       eigvalue=load('NeuroType_Case/Case_1/Couple=0.06/eig_21_10000000.txt');
       elseif Case_Type_Choose==5
       eigvalue=load('NeuroType_Case/Case_2/Couple=0.06/eig_21_10000000.txt');
       end
    
    elseif ie==5
      if Case_Type_Choose==0
       eigvalue=load('HH_Three_280/Couple=0.08/eig_21_10000000.txt');
       elseif Case_Type_Choose==1
       eigvalue=load('HH_Three_280290371/Couple=0.08/eig_21_10000000.txt');   
       elseif Case_Type_Choose==2
       eigvalue=load('HH_Three_280290371_2/Couple=0.08/eig_21_10000000.txt');   
       elseif Case_Type_Choose==3
       eigvalue=load('HH_Three_280290371_Q/Couple=0.08/eig_21_10000000.txt'); 
       elseif Case_Type_Choose==4
       eigvalue=load('NeuroType_Case/Case_1/Couple=0.08/eig_21_10000000.txt');
       elseif Case_Type_Choose==5
       eigvalue=load('NeuroType_Case/Case_2/Couple=0.08/eig_21_10000000.txt');
       end
    
    elseif ie==6
      if Case_Type_Choose==0
       eigvalue=load('HH_Three_280/Couple=0.1/eig_21_10000000.txt');
       elseif Case_Type_Choose==1
       eigvalue=load('HH_Three_280290371/Couple=0.1/eig_21_10000000.txt');   
       elseif Case_Type_Choose==2
       eigvalue=load('HH_Three_280290371_2/Couple=0.1/eig_21_10000000.txt');   
       elseif Case_Type_Choose==3
       eigvalue=load('HH_Three_280290371_Q/Couple=0.10/eig_21_10000000.txt'); 
       elseif Case_Type_Choose==4
       eigvalue=load('NeuroType_Case/Case_1/Couple=0.10/eig_21_10000000.txt');
       elseif Case_Type_Choose==5
       eigvalue=load('NeuroType_Case/Case_2/Couple=0.10/eig_21_10000000.txt');
       end
     
    elseif ie==7
      if Case_Type_Choose==0
       eigvalue=load('HH_Three_280/Couple=0.11/eig_21_10000000.txt');
       elseif Case_Type_Choose==1
       eigvalue=load('HH_Three_280290371/Couple=0.11/eig_21_10000000.txt');   
       elseif Case_Type_Choose==2
       eigvalue=load('HH_Three_280290371_2/Couple=0.11/eig_21_10000000.txt');   
       elseif Case_Type_Choose==3
       eigvalue=load('HH_Three_280290371_Q/Couple=0.11/eig_21_10000000.txt'); 
       elseif Case_Type_Choose==4
       eigvalue=load('NeuroType_Case/Case_1/Couple=0.11/eig_21_10000000.txt');
       elseif Case_Type_Choose==5
       eigvalue=load('NeuroType_Case/Case_2/Couple=0.11/eig_21_10000000.txt');
       end
     
    elseif ie==8
       if Case_Type_Choose==0
       eigvalue=load('HH_Three_280/Couple=0.12/eig_21_10000000.txt');
       elseif Case_Type_Choose==1
       eigvalue=load('HH_Three_280290371/Couple=0.12/eig_21_10000000.txt');   
       elseif Case_Type_Choose==2
       eigvalue=load('HH_Three_280290371_2/Couple=0.12/eig_21_10000000.txt');   
       elseif Case_Type_Choose==3
       eigvalue=load('HH_Three_280290371_Q/Couple=0.12/eig_21_10000000.txt'); 
       elseif Case_Type_Choose==4
       eigvalue=load('NeuroType_Case/Case_1/Couple=0.12/eig_21_10000000.txt');
       elseif Case_Type_Choose==5
       eigvalue=load('NeuroType_Case/Case_2/Couple=0.12/eig_21_10000000.txt');
       end
        
    elseif ie==9
      if Case_Type_Choose==0
       eigvalue=load('HH_Three_280/Couple=0.13/eig_21_10000000.txt');
       elseif Case_Type_Choose==1
       eigvalue=load('HH_Three_280290371/Couple=0.13/eig_21_10000000.txt');   
       elseif Case_Type_Choose==2
       eigvalue=load('HH_Three_280290371_2/Couple=0.13/eig_21_10000000.txt');   
       elseif Case_Type_Choose==3
       eigvalue=load('HH_Three_280290371_Q/Couple=0.13/eig_21_10000000.txt'); 
       elseif Case_Type_Choose==4
       eigvalue=load('NeuroType_Case/Case_1/Couple=0.14/eig_21_10000000.txt');
       elseif Case_Type_Choose==5
       eigvalue=load('NeuroType_Case/Case_2/Couple=0.14/eig_21_10000000.txt');
       end
    
    elseif ie==10
      if Case_Type_Choose==0
       eigvalue=load('HH_Three_280/Couple=0.14/eig_21_10000000.txt');
       elseif Case_Type_Choose==1
       eigvalue=load('HH_Three_280290371/Couple=0.14/eig_21_10000000.txt');   
       elseif Case_Type_Choose==2
       eigvalue=load('HH_Three_280290371_2/Couple=0.14/eig_21_10000000.txt');   
       elseif Case_Type_Choose==3
       eigvalue=load('HH_Three_280290371_Q/Couple=0.14/eig_21_10000000.txt');
       elseif Case_Type_Choose==4
       eigvalue=load('NeuroType_Case/Case_1/Couple=0.14/eig_21_10000000.txt');
       elseif Case_Type_Choose==5
       eigvalue=load('NeuroType_Case/Case_2/Couple=0.14/eig_21_10000000.txt');
       end
        
    elseif ie==11
       if Case_Type_Choose==0
       eigvalue=load('HH_Three_280/Couple=0.15/eig_21_10000000.txt');
       elseif Case_Type_Choose==1
       eigvalue=load('HH_Three_280290371/Couple=0.15/eig_21_10000000.txt');   
       elseif Case_Type_Choose==2
       eigvalue=load('HH_Three_280290371_2/Couple=0.15/eig_21_10000000.txt');   
       elseif Case_Type_Choose==3
       eigvalue=load('HH_Three_280290371_Q/Couple=0.15/eig_21_10000000.txt');
       elseif Case_Type_Choose==4
       eigvalue=load('NeuroType_Case/Case_1/Couple=0.15/eig_21_10000000.txt');
       elseif Case_Type_Choose==5
       eigvalue=load('NeuroType_Case/Case_2/Couple=0.15/eig_21_10000000.txt');
       end
       
       elseif ie==12
       if Case_Type_Choose==0
       eigvalue=load('HH_Three_280/Couple=0.16/eig_21_10000000.txt');
       elseif Case_Type_Choose==1
       eigvalue=load('HH_Three_280290371/Couple=0.16/eig_21_10000000.txt');   
       elseif Case_Type_Choose==2
       eigvalue=load('HH_Three_280290371_2/Couple=0.16/eig_21_10000000.txt');   
       elseif Case_Type_Choose==3
       eigvalue=load('HH_Three_280290371_Q/Couple=0.16/eig_21_10000000.txt');
       elseif Case_Type_Choose==4
       eigvalue=load('NeuroType_Case/Case_1/Couple=0.16/eig_21_10000000.txt');
       elseif Case_Type_Choose==5
       eigvalue=load('NeuroType_Case/Case_2/Couple=0.16/eig_21_10000000.txt');
       end
   end
   
   information_dim=0;
   sum_le=0;
   for j=1:length(eigvalue)
      sum_le_pre=sum_le;
      sum_le=sum_le+eigvalue(j);
      if sum_le<0 && sum_le_pre>0
          information_dim=(j-1)+sum_le_pre/abs(eigvalue(j));
      end    
   end
   
   Taken_dim(ie,1)=2*information_dim+1;
   
end
       
       
       
       
       
       
       