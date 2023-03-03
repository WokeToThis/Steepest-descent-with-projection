syms x y 
f = 1/3*x^2 + 3*y^2

fsurf(f,[-3 3])
%fcontour(f,[-1 1])



g = gradient(f,[x,y])

[X, Y] = meshgrid(-1:.1:1,-1:.1:1);
G1 = subs(g(1),[x y],{X,Y});
G2 = subs(g(2),[x y],{X,Y});
fcontour(f,[-1 1])
hold on;
quiver(X,Y,G1,G2)
figure;


% [xk,yk,k] = Steep(f,1,1,0.1,1,0.001)
% [xxk,yyk,kk] = Steep(f,1,1,0.3,1,0.001)
% [xxxk,yyyk,kkk] = Steep(f,1,1,3,1,0.001)
% [x1k,y1k,k1] = Steep(f,1,1,5,1,0.001)
% [x1k,y1k,k1] = Steep(f,5,-5,0.5,1,0.01)


% [x1k,y1k,k1] = SteepProj(f,5,-5,0.5,1,0.01,5)
% figure;
% [x1kk,y1kk,kk1] = SteepProj(f,-5,10,0.1,1,0.01,15)
% figure;
% [x1kkk,y1kkk,kkk1] = SteepProj(f,8,-10,0.2,1,0.01,0.1)
%% 
% Steepest Descent

function [xk,yk,k] = Steep(f,xk,yk,gk1,mode,epsilon)
%mode is about gk 1<-const 2<-minimizing f(x - gk*grad(f(x))) 3<-Armijo
k=0;
max_it = 100;

syms x y gk
    
    
   while(k<max_it)
      df = gradient(f);
      dfx = df(1,1);
      dfy = df(2,1);
%       dfx = diff(f,x);
%       dfy = diff(f,y);
        

    gradx = double(subs(dfx,[x y],[xk yk]));
    grady = double(subs(dfy,[x y],[xk yk]));
    if(mode==1)
        gk=gk1;
    end
    if(mode==2)
        g = @(gk)  subs(f,[x y],[xk - gk*gradx,yk - gk*grady]);
        gk = GoldenSec(g,-5,5,0.01,0.01,200,0.618);
       

     
    end
    if(mode == 3)
        gk = Armijo(f,xk,yk,gradx,grady,-df,df,1,0.007,0.5);
    end

    


       if(sqrt(gradx^2+grady^2)<epsilon) 
           break;
       end
    

        xk =double( xk - gk*gradx);
        yk = double(yk - gk*grady);
        
        k = k+1;
       kvec(k) = k;
       fvec(k) = double(subs(f,[x y],[xk yk]));

       plot(kvec,fvec);
       xlabel('Iterations')
       ylabel('F(xk,yk)')
       title('SteepestDescent')
   end

end
%% 
% Steepest Descent with Projection

function [xk,yk,k] = SteepProj(f,xk,yk,gk1,mode,epsilon,sk)

k=0;
max_it = 100;

syms x y gk
    
    
   while(k<max_it)
      df = gradient(f);
      dfx = df(1,1);
      dfy = df(2,1);
%       dfx = diff(f,x);
%       dfy = diff(f,y);
        

    gradx = double(subs(dfx,[x y],[xk yk]));
    grady = double(subs(dfy,[x y],[xk yk]));
    if(mode==1)
        gk=gk1;
    end
    if(mode==2)
        g = @(gk)  subs(f,[x y],[xk - gk*gradx,yk - gk*grady]);
        gk = GoldenSec(g,-5,5,0.01,0.01,200,0.618);
       

     
    end
    if(mode == 3)
        gk = Armijo(f,xk,yk,gradx,grady,-df,df,1,0.007,0.5);
    end

    


       if(sqrt(gradx^2+grady^2)<epsilon) 
           break;
       end
    

       [xxk,yyk] = Prx(f,xk,yk,-10,5,-8,12,sk);
       xk =double( xk + gk*(xxk - xk));
       yk = double(yk + gk*(yyk -yk));
        
        
        k = k+1;
       kvec(k) = k;
       fvec(k) = double(subs(f,[x y],[xk yk]));
       xvec(k) = xk;
       yvec(k) = yk;

       plot(kvec,fvec);
       xlabel('Iterations')
       ylabel('F(xk,yk)')
       title('SteepestDescent')

       
   end
   
       figure;
       fcontour(f,[-10 10])
       hold on;
       plot(xvec,yvec)

    

end
%% 
% Projection

function [xk,yk] = Prx(f,xk,yk,a1,b1,a2,b2,sk)

syms x y;

 df = gradient(f);
 dfx = df(1,1);
 dfy = df(2,1);
%       dfx = diff(f,x);
%       dfy = diff(f,y);
        

  gradx = double(subs(dfx,[x y],[xk yk]));
  grady = double(subs(dfy,[x y],[xk yk]));

   tmpxk =double(xk - sk*gradx);
   tmpyk = double(yk - sk*grady);

   
   
   if(tmpxk <= a1)
       xk = a1;
   
   elseif(tmpxk >= b1)
           xk = b1;
   else
       xk = tmpxk;
       
   end

   if(tmpyk <= a2)
       yk = a2;
   
   elseif(tmpyk >= b2)
           yk = b2;
   else
       yk = tmpyk;
       
   end

   

  



end