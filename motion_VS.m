float  S, V, Acc, Dec, Jrek, Tacc, Tdec, T1, T2, T3, C1, C2,C3,C4,C5,C6,J1,J2;
float motion1()   /*V > Acc * Tacc> Dec * Tdec*,S > C1,double Acc, double Tacc, double Dec, double Tdec, double V, double C1, double S 无匀速运动，*/
{
	T1 = (V - Acc*Tacc) / Acc;
	T2 = (V - Dec*Tdec) / Dec;
	C1 = V*Tacc + V*T1 / 2 + V*Tdec + V*T2 / 2;
	T3 = (S - C1) / V;
	return 0;
}
float motion2()  
{
	T1 = V / Acc - Tacc;
	T2 = V / Dec - Tdec;
	T3 = 0;
	return 0;
}
float motion3()  
{
	J1 = V / (Tacc*Tacc);
	T1 = 0;
	T2 =( J1*Tacc*Tacc*Tdec/V-Tdec)/(1-(J1*Tacc*Tacc)/V);
	T3 = 0;
	return 0;
}
float motion4()
{
	T1 = 0;
	T2 = 0;
	T3 = 0;
	Tdec = pow(S / (2 * Tdec) / (1 + Tacc / Tdec), 1 / 3);
	return 0;
}
float motion5()
{
	C4 = V*Tacc + V*Tdec + V*T2 / 2;
	T1 = 0;
	T2=(V-Tdec*Dec)/Dec;
	T3 = (S - C4) / V;
	return 0;
}
float motion6()
{
	J1 = V / (Tacc*Tacc);
	T1 = 0;
	T2 = (J1*Tacc*Tacc*Tdec / V - Tdec) / (1 - (J1*Tacc*Tacc) / V);
	T3 = 0;
	return 0;
}
float motion7()
{
	T1 = 0;
	T2 = 0;
	T3 = 0;
	return 0;
}
float motion8()
{
	T1 = 0;
	T2 = 0;
	C6 = V*Tacc + V*Tdec;
	T3 = (S-C6)/V;
	return 0;
}
float motion9()
{
	T1 = 0;
	T2 = 0;
	T3 = 0;
	return 0;
}
int main()
{
	printf("依次输入V,Acc,Tacc,Dec,Tdec,S\n");
	scanf("%f", &V);
	scanf("%f", &Acc);
	scanf("%f", &Tacc);
	scanf("%f", &Dec);
	scanf("%f", &Tdec);
	scanf("%f", &S);
	int flag=0;
	if (V > Acc * Tacc > Dec * Tdec)
		if (S > C1)
			flag = 1;
		else
			if (C1 >= S > C2)
				flag = 2;
			else
				if (C2 >= S > C3)
					flag = 3;
				else
					flag = 4;
	else
		if (Acc * Tacc >= V > Dec * Tdec)
			if (S > C4)
				flag = 5;
			else
				if (C4 >= S > C5)
					flag = 6;
				else
					flag = 7;
		else
			if (Acc * Tacc >= Dec * Tdec > V)
				if (S > C6)
					flag = 8;
				else
					flag = 9;
			else
				printf("输入参数错误");
	switch(flag)
	{
	case 1:motion1, printf("%f", T1), printf("%f", T2), printf("%f", T3);
		break;
	case 2:motion2, printf("%f", T1), printf("%f", T2), printf("%f", T3);
		break;
	case 3:motion3,printf("%f", T1), printf("%f", T2), printf("%f", T3);
		 break;
	case 4:motion4,printf("%f", T1), printf("%f", T2), printf("%f", T3);
		break;
	case 5:motion5,printf("%f", T1), printf("%f", T2), printf("%f", T3);
		break;
	case 6:motion6,printf("%f", T1), printf("%f", T2), printf("%f", T3);
		break;
	case 7:motion7,printf("%f", T1), printf("%f", T2), printf("%f", T3);
		break;
	case 8:motion8,printf("%f", T1), printf("%f", T2), printf("%f", T3);
		break;
	case 9:motion9,printf("%f", T1), printf("%f", T2), printf("%f", T3);
		break;
	default:
		break;
	}
	return 0;
}