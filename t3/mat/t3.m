format long;

R1 = 120000;
C = 0.007;
n = 10; %razao de espiras no transformador
eta = 1; %constante do material - usada na formula dos diodos
freq = 50;
I_s = 1e-14; %constante para calcular a tensao nos diodos
num_diodos = 22; %numero de diodos
Vs = 230; %tensao inicial no source
Vt = 0.025; %constante para o calculo da corrente nos díodos, thermal voltage
w = 2*pi*freq;

T = 1/(2*freq);
t_off = (1/4)*T;

R2 = 678100;

for i = 1:20
  func = (Vs/n)*C*w*sin(w*t_off) - (1/R1)*(Vs/n)*cos(w*t_off) - I_s*(exp(12/(eta*Vt*num_diodos))-1);
  dfunc = (Vs/n)*C*(w^2)*cos(w*t_off)+(1/R1)*(Vs/n)*w*sin(w*t_off);
  t_off = t_off - (func/dfunc);
endfor

t_on = (3/4)*T;

Req = 1/((1/R1)+(1/R2)); 



for i = 1:20
  func = (Vs/n)*cos(w*t_on)+(Vs/n)*cos(w*t_off)*exp(-(1/(Req*C))*(t_on-t_off));
  dfunc = -w*(Vs/n)*sin(w*t_on)-(Vs/n)*cos(w*t_off)*(1/(Req*C))*exp(-(1/(Req*C))*(t_on-t_off));
  t_on = t_on - (func/dfunc);
endfor

t = 0:(1e-4):0.2;

l = length(t);

v0envelope = ones(1,l);

for i = 1:l
  if t(i)<=t_off
    v0envelope(i) = abs((230/n)*cos(w*t(i)));
  else
    if t(i)<=t_on
      v0envelope(i) = (230/n)*abs(cos(w*t_off))*exp(-(1/(Req*C))*(t(i)-t_off));
    else
      t_off = t_off + T;
      t_on = t_on + T;
      if t(i)<=t_off
        v0envelope(i) = abs((230/n)*cos(w*t(i)));
      else
        if t(i)<=t_on
          v0envelope(i) = (230/n)*abs(cos(w*t_off))*exp(-(1/(Req*C))*(t(i)-t_off));
        endif
      endif
    endif
  endif    
endfor



v0envelope_dc = mean(v0envelope);
ripple_v0envelope = max(v0envelope) - min(v0envelope)
v0envelope_centro = (ripple_v0envelope/2) + min(v0envelope); 

rd = (eta*Vt)/(I_s*exp((12/num_diodos)/(eta*Vt)));
v0regulator_ac = ((num_diodos*rd)/(num_diodos*rd + R2))*(v0envelope - v0envelope_dc);

if v0envelope_centro >= 12
    v0regulator_dc = 12;
else
    v0regulator_dc = v0envelope_centro; 
endif

v0regulator = v0regulator_ac + v0regulator_dc;

aux2 = v0regulator - 12;
average = mean(aux2);
DCregulator = mean (v0regulator)
ripple = max(v0regulator) - min(v0regulator);

Cost = (R1+R2)/1000 + C*(10^6) + 0.1*(num_diodos + 4);
Merit = 1/(Cost*(ripple+abs(average) + 10^(-6)))

aux=abs((230/n)*cos(w*t));

f1=figure();
plot(t*1000,v0regulator,"linewidth",4,";V_{Reg};")
hold on
plot(t*1000,v0envelope,"linewidth",4,";V_{Env};")
hold on
plot(t*1000,aux,"linewidth",4,";V_{Trans};")
xlabel("t(ms)");
ylabel("V(Volts)");
legend('Location','northeast');
print(f1,"graph1.eps","-depsc");

f2=figure();
plot(t*1000,v0envelope -12,"linewidth",2,";V_{Env}-12;")
hold on
plot(t*1000,v0regulator -12,"linewidth",2,";V_{Reg}-12;")
xlabel("t(ms)");
ylabel("V(Volts)");
legend('Location','northeast');
print(f2,"graph2.eps","-depsc");

printf("Envelope_TAB\n");
printf("Ripple Envelope = %e V\n",ripple_v0envelope);
printf("Average Envelope = %e V\n",v0envelope_dc);
printf("Envelope_END\n");

printf("Regulator_TAB\n");
printf("Ripple Regulator = %e V\n",ripple);
printf("Average Regulator = %e V\n",DCregulator);
printf("Regulator_END\n");

