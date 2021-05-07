format long;

R1 = 98437;
C = 0.0073476945;
n = 10; %numero de espiras no transformador
eta = 1; %constante do material - usada na formula dos diodos
freq = 50;
I_s = 1e-14; %constante para calcular a tensao nos diodos
k = 22; %numero de diodos
%Vd = 0.6;
Vs = 230; %tensao inicial no source
Vt = 0.025; %constante para o calculo da corrente nos díodos, thermal voltage
w = 2*pi*freq;

T = 1/(2*freq);
t_off = (1/4)*T;

R2 = 677410;

for i = 1:20
  f = (Vs/n)*C*w*sin(w*t_off) - (1/R1)*(Vs/n)*cos(w*t_off) - I_s*(exp(12/(eta*Vt*k))-1);
  fl = (Vs/n)*C*(w^2)*cos(w*t_off)+(1/R1)*(Vs/n)*w*sin(w*t_off);
  t_off = t_off - (f/fl);
endfor

t_on = (3/4)*T;

Req = 1/((1/R1)+(1/R2)); %desprezar resistencia dos diodos?

%duvidas: sinal menos/mais linha 33 // onde obter o valor de VT-temperatura

for i = 1:20
  f = (Vs/n)*cos(w*t_on)+(Vs/n)*cos(w*t_off)*exp(-(1/(Req*C))*(t_on-t_off));
  fl = -w*(Vs/n)*sin(w*t_on)-(Vs/n)*cos(w*t_off)*(1/(Req*C))*exp(-(1/(Req*C))*(t_on-t_off));
  t_on = t_on - (f/fl);
endfor

t = 0:(1e-4):0.2;

l = length(t);

v0_env = ones(1,l);

for i = 1:l
  if t(i)<=t_off
    v0_env(i) = abs((230/n)*cos(w*t(i)));
  else
    if t(i)<=t_on
      v0_env(i) = (230/n)*abs(cos(w*t_off))*exp(-(1/(Req*C))*(t(i)-t_off));
    else
      t_off = t_off + T;
      t_on = t_on + T;
      if t(i)<=t_off
        v0_env(i) = abs((230/n)*cos(w*t(i)));
      else
        if t(i)<=t_on
          v0_env(i) = (230/n)*abs(cos(w*t_off))*exp(-(1/(Req*C))*(t(i)-t_off));
        endif
      endif
    endif
  endif    
endfor



v0_env_dc = mean(v0_env);
ripple_v0_env = max(v0_env) - min(v0_env)
v0_env_centro = (ripple_v0_env/2) + min(v0_env); %centro AC?

rd = (eta*Vt)/(I_s*exp((12/k)/(eta*Vt)));
v0_reg_ac = ((k*rd)/(k*rd + R2))*(v0_env - v0_env_dc);

if v0_env_centro >= 12
    v0_reg_dc = 12;
else
    v0_reg_dc = v0_env_centro; %centro em vez de dc?
endif

v0_reg = v0_reg_ac + v0_reg_dc;

testar = v0_reg - 12;
average = mean(testar);
DC_level = mean (v0_reg)
ripple = max(v0_reg) - min(v0_reg);

Cost = (R1+R2)/1000 + C*(10^6) + 0.1*(k + 4);
Merit = 1/(Cost*(ripple+abs(average) + 10^(-6)))

aux=abs((230/n)*cos(w*t));

f1=figure();
plot(t*1000,v0_reg,"linewidth",4,";V_{Reg};")
hold on
plot(t*1000,v0_env,"linewidth",4,";V_{Env};")
hold on
plot(t*1000,aux,"linewidth",4,";V_{Trans};")
xlabel("t(ms)");
ylabel("V(Volts)");
legend('Location','northeast');
print(f1,"graph1.eps","-depsc");

f2=figure();
plot(t*1000,v0_env -12,"linewidth",2,";V_{Env}-12;")
hold on
plot(t*1000,v0_reg -12,"linewidth",2,";V_{Reg}-12;")
xlabel("t(ms)");
ylabel("V(Volts)");
legend('Location','northeast');
print(f2,"graph2.eps","-depsc");

printf("Envelope_TAB\n");
printf("Ripple Envelope = %e V\n",ripple_v0_env);
printf("Average Envelope = %e V\n",v0_env_dc);
printf("Envelope_END\n");

printf("Regulator_TAB\n");
printf("Ripple Regulator = %e V\n",ripple);
printf("Average Regulator = %e V\n",DC_level);
printf("Regulator_END\n");

