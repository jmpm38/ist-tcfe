R1 = 1000;
R2 = 1000;
R3 = 150000;
R4 = 1000;
C1 = 220e-9;
C2 = 110e-9;
i = sqrt(-1);

wL = 1/(R1*C1);
wH = 1/(R2*C2);

if wL>wH
	aux = wL;
	wL = wH;
	wH = aux;
endif

w0 = sqrt(wL*wH);

input_imp = R1 + 1/(i*w0*C1);
output_imp = R2/(1+i*w0*C2*R2);


T_w0 = ((R1*C1*w0*i)/(1+R1*C1*w0*i))*(1/(1+R2*C2*w0*i))*(1+(R3/R4));
Gain_w0 = abs(T_w0);
Gain_w0_dB = 20*log10(Gain_w0);

f = logspace(1,8,70);
w = 2*pi*f;
T = ones(1,length(w));
T_dB = ones(1,length(w));

for k = 1:length(w)
	T(k) = ((R1*C1*w(k)*i)/(1+R1*C1*w(k)*i))*(1/(1+R2*C2*w(k)*i))*(1+(R3/R4))
	T_dB(k) = 20*log10(abs(T(k)));
endfor

phi = (180/pi)*arg(T);

freq_deviation = w0/(2*pi) - 1000;
gain_deviation = Gain_w0_dB - 40;

Cost = R1/1000 + R2/1000 + 3*100 + R4/1000 + (3*220) + 13323;

Merit = 1/(Cost*(abs(freq_deviation) + abs(gain_deviation) + 1e-6));

teorica = figure ();
plot(log10(f),T_dB,"g");
legend("v_o(f)/v_I(f)");

xlabel ("Log10(Frequency [Hz])");
ylabel ("Circuit Gain");

print (teorica, "teoria", "-depsc");

fase = figure ();
plot(log10(f),phi,"g");
legend("arg(v_o(f)/v_I(f))");

xlabel ("Log10(Frequency [Hz])");
ylabel ("Phase Degrees");

print (fase, "fase", "-depsc");

printf ("resultados_TAB\n");
printf ("Central Frequency Deviation = %e Hz \n", freq_deviation);
printf ("Gain Deviation = %e dB \n", gain_deviation);
printf ("Cost = %e MU \n", Cost);
printf ("Merit = %e \n", Merit);
printf ("resultados_END\n\n");

printf ("frequencias_TAB\n");
printf ("Lower Cut Off Frequency = %e rad/s\n", wL);
printf ("Higher Cut Off Frequency = %e rad/s\n", wH);
printf ("Central Frequency = %e rad/s\n", w0);
printf ("frequencias_END\n\n");

printf ("comparacao_TAB\n");
printf ("Gain = %e dB \n", Gain_w0_dB);
printf ("Central Frequency = %e Hz\n", w0/(2*pi));
printf ("Gain Deviation = %e dB \n", gain_deviation);
printf ("Central Frequency Deviation = %e Hz \n", freq_deviation);
printf ("comparacao_END\n\n");

printf ("ganhos_TAB\n");
printf ("Gain = %e \n", Gain_w0);
printf ("Gain = %e dB \n", Gain_w0_dB);
printf ("ganhos_END\n\n");

printf ("impedancias_TAB\n");
printf ("Input Impedance = %e + %ej \n", real(input_imp), imag(input_imp));
printf ("Output Impedance = %e + %ej\n", real(output_imp), imag(output_imp));
printf ("impedancias_END\n\n");

printf ("componentes_TAB\n");
printf ("R1 = %e \n", R1);
printf ("R2 = %e \n", R2);
printf ("R3 = %e \n", R3);
printf ("R4 = %e \n", R4);
printf ("C1 = %e \n", C1);
printf ("C2 = %e \n", C2);
printf ("componentes_END\n\n");
