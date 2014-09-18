(* ::Package:: *)

(* ::Input:: *)
R=1;
\[Gamma] = 150.0/0.510998928+1;
\[Beta] = Sqrt[1-1/\[Gamma]^2];
c = 299792458;
\[Omega] = (\[Beta] c)/R;
\[Phi]=(5 \[Pi])/180;
L = 0.02;
\[Sigma]=0.0007;
k=(R \[Phi] + L)/(\[Beta] c);
m=(R \[Phi])/(\[Beta] c);
\[CapitalDelta][x_]:=c (k-x)
dx=9.77517 10^-5;
dy = 0.000393701;
lx=0.10132819;
ux=0.11256963;
ly=-0.00076224908;
uy=0.011836183;

Efldx[x_,y_,t_]:=If[(c(k-t))^2-(x-R Sin[\[Omega] t])^2-(y-R+R Cos[\[Omega] t])^2>=0,(\[Beta]^2) /(R (dy/4 +\[CapitalDelta][t]-\[Beta] Cos[\[Omega] t](x-R Sin[\[Omega] t])-\[Beta] Sin[\[Omega] t](y-R+R Cos[\[Omega] t]))^3) ((-x+R Sin[t \[Omega]]) (-R+(R-y) Cos[t \[Omega]]+x Sin[t \[Omega]])-\[Beta] (-R+y+R Cos[t \[Omega]]) \[CapitalDelta][t]+Sin[t \[Omega]] \[CapitalDelta][t]^2),0]

Efldy[x_,y_,t_]:=If[(c(k-t))^2-(x-R Sin[\[Omega] t])^2-(y-R+R Cos[\[Omega] t])^2>=0,(\[Beta]^2) /(R (dy/4 +\[CapitalDelta][t]-\[Beta] Cos[\[Omega] t](x-R Sin[\[Omega] t])-\[Beta] Sin[\[Omega] t](y-R+R Cos[\[Omega] t]))^3) (-(-R+y+R Cos[t \[Omega]]) (-R+(R-y) Cos[t \[Omega]]+x Sin[t \[Omega]])+\[Beta] (x-R Sin[t \[Omega]]) \[CapitalDelta][t]-Cos[t \[Omega]] \[CapitalDelta][t]^2),0]

Q=Table[(Erf[(x+dx/2)/(Sqrt[2]\[Sigma])]-Erf[(x-dx/2)/(Sqrt[2]\[Sigma])])/2 (Erf[(y+dy/2)/(Sqrt[2]\[Sigma])]-Erf[(y-dy/2)/(Sqrt[2]\[Sigma])])/2,{y,-8 dy,7 dy,dy},{x,-29 dx,28 dx,dx}];

For[l=0;t=0;dt=m/100;Ex=Table[0,{y,ly,uy,dy},{x,lx,ux,dx}],t<=m,l++;t+=dt,Print[l];
	For[j=1,j<=16,j++,
		For[i=1,i<=56,i++,
			Ex=Ex+dt Q[[j,i]]Table[Efldx[x-(i-29) dx,y-(j-9) dy,t],{y,ly,uy,dy},{x,lx,ux,dx}]]]]

Min[Ex]
Max[Ex]

Export["Eacc-Field.csv",Ex]

For[l=0;t=0;dt=m/100;Ey=Table[0,{y,ly,uy,dy},{x,lx,ux,dx}],t<=m,l++;t+=dt,Print[l];
	For[j=1,j<=16,j++,
		For[i=1,i<=56,i++,
			Ey=Ey+dt Q[[j,i]]Table[Efldy[x-(i-29) dx,y-(j-9) dy,t],{y,ly,uy,dy},{x,lx,ux,dx}]]]]

Min[Ey]

Max[Ey]

Export["Eyacc-Field.csv",Ey]

xend=R Sin[\[Phi]]+L Cos[\[Phi]];
yend =R(1-Cos[\[Phi]])+L Sin[\[Phi]];

f=ListInterpolation[Ey,{{ly,uy},{lx,ux}}];

yonline=Table[f[Sin[\[Phi]]x+yend,Cos[\[Phi]]x+xend],{x,(-(ux-lx)+dx)/2,((ux-lx)-dx)/2,dx}];
Export["Eyacc-online.csv",yonline]


g=ListInterpolation[Ex,{{ly,uy},{lx,ux}}];

xonline=Table[g[Sin[\[Phi]]x+yend,Cos[\[Phi]]x+xend],{x,(-(ux-lx)+dx)/2,((ux-lx)-dx)/2,dx}];
Export["Exacc-online.csv",xonline]

Exonline=Table[Cos[\[Phi]]xonline[[i]]-Sin[\[Phi]]yonline[[i]],{i,1,114}];
