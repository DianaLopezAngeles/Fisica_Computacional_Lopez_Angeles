
__precompile__()

module herramientas

#Cargamos nuestros programas de la tarea 4:
export metodo_newton_arbitrario
function metodo_newton_arbitrario(f,df,x0,t)
    x = x0
    for i in 1:100                             
        x=x-(f(x)/df(x))                       
    end
return x
end;

export metodo_newton_intervalo
function metodo_newton_intervalo(f,df,intervalo)    
    list=zeros(length(intervalo))                   
    x = 0                                           
    for i in 1:length(intervalo)                    
        x = intervalo[i]                            
        for n in 1:200                              
            x = x-(f(x)/df(x))                      
            end                                     
    list[i]=x;                                      
    end
    list                                            
end;

export metodo_newton_epsilon
function metodo_newton_epsilon(f,df,intervalo)          
    t = []                                          
    epsilon = 10.0^-8                               
    list = metodo_newton_intervalo(f,df,intervalo)  
    push!(t,list[1])                                
    for i in 1:length(t)                            
        for n in 1:length(list)                     
            if abs(t[i]-list[n]) > epsilon          
                push!(t,list[n])                    
                end                               
        end
        return t                                            
    end
end

#Cargamos los programas de la tarea 5: Metodos de integracion
export sumas_Riemann
function sumas_Riemann(f,a,b,subintervalos)         
    N = subintervalos                                 
    A = a
    B = b
    datos = []                                      
    suma = 0                                          
    area_bajo_la_curva = 0
    for i in 1:N                                       
        b = ((B-A)/N)*(i) - A                          
        a = A + ((B-A)/N)*(i-1)
        area_bajo_la_curva = (b-a)*(f((b+a)/2))        
        push!(datos,area_bajo_la_curva)              
    end
    suma = sum(datos,1:N)                             
return suma
end;

export metodo_trapecio
function metodo_trapecio(f,a,b,subintervalos)            
    N = subintervalos                                  
    A = a
    B = b
    datos = []                                         
    suma = 0                                            
    area_bajo_la_curva = 0

    for i in 1:N
        b = ((B-A)/N)*(i) - A
        a = A + ((B-A)/N)*(i-1)
        area_bajo_la_curva = (b-a)*((f(b)+f(a))/2) 
        push!(datos,area_bajo_la_curva)
    end
    suma = sum(datos,1:N)
return suma
end;

#Cargamos los programas de la tarea 7:
export metodo_simpson
function metodo_simpson(f,a,b,subintervalos)           
    N = subintervalos                                  
    A = a
    B = b
    datos = []                                         
    suma = 0                                           
    area_bajo_la_curva = 0

    for i in 1:N
        b = ((B-A)/N)*(i) - A
        a = A + ((B-A)/N)*(i-1)
        area_bajo_la_curva = (b-a)/6*((f(a)+4f((a+b)/2))+f(b))
        push!(datos,area_bajo_la_curva)
    end
    suma = sum(datos,1:N)
return suma
end;

#Cargamos programas de la tarea 8:

export derivada_numerica
function derivada_numerica(f,x,h)
    df_ = (f(h+x)-f(x))/h 
    return df_ 
end;

#Cargamos programas de la tarea 11:

export metodo_euler
function metodo_euler(f,x0,t0,tf,h) 
    n=round((tf-t0)/h)+1                 
    listt=linspace(t0,tf,n)             
    listx=zeros(n)         
    listx[1]=x0             
    for i in 1:length(listx)-1   
        listx[i+1]=listx[i]+h*f(listx[i],listt[i])
    end
    return listt, listx
end;

export metodo_euler_dimensiones
function metodo_euler_dimensiones(f,list,x0)
     x = x0
     h = list[2]-list[1]
     listx=[]
     push!(listx,x)
     for i in 2:length(list)
        t = i*h
        x = x + f(x,t)*h
        push!(listx,x) 
     end
     return listx
    end;

#Cargamos programas de la tarea 12:
export euler_implicito
function euler_implicito(f,df,x0,t0,tf,h)
    n=round((tf-t0)/h)+1                 
    listt=linspace(t0,tf,n)
    listx = zeros(length(listt))
    x = x0 
    listx[1] = x0
    for i in 2:length(listt)
        g(x) = x - listx[i-1] - h*f(x,listt[i])
        dg(x) = 1 - h*df(x,listt[i])
        x = metodo_newton_arbitrario(g,dg,listx[i-1],listt[i])
        listx[i] = x
    end
    return (listt,listx)
end;

export euler_punto_medio
function euler_punto_medio(f,x0,t0,tf,h)
    n=round((tf-t0)/h)+1                 
    listt=linspace(t0,tf,n)             
    listx=zeros(n)
    listx[1]=x0                              
    for i in 1:length(listx)-1               
        listx[i+1]=listx[i]+h*f((listx[i]+(h/2)f(listx[i],listt[i])),listt[i+1])
    end
    return listt, listx
end;

export runge_kutta_4
function runge_kutta_4(f,x0,t0,tf,h)
    n=round((tf-t0)/h)+1               
    listt=linspace(t0,tf,n) 
    listx=zeros(n)  
    listx[1] = x0 
    for i in 1:length(listx)-1
        k1 = f(listx[i], listt[i])
        k2 = f(listx[i] + h*(k1)/2, listt[i+1])
        k3 = f(listx[i] + h*(k2)/2, listt[i+1])
        k4 = f(listx[i] + h*(k3), listt[i],)
        listx[i+1] = listx[i] + h/6*(k1 + 2*(k2) + 2*(k3) + k4)
    end
    return listt,listx
end;

export runge_kutta_dimensiones
function runge_kutta_dimensiones(f,list,x0)
     x = x0
     h = list[2]-list[1]
     listx=[]
     push!(listx,x)
     for i in 2:length(list)
        t = i*h
        k1 = f(x,t);
        k2 = f(x+(h/2)*k1,t+(h/2));
        k3 = f(x+(h/2)*k2, t+(h/2));
        k4 = f(x+h*k3, t+h);
        x = x+(h/6)*(k1+2*k2+2*k3+k4);
        push!(listx,x) 
     end
     return listx
end;
end


