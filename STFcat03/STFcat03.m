clear all
close all
clc
m= load('ncl.dat')

lambdaHPP = 34.2
lambdaNHPP1 =[35.32,2.32*2*pi,-0.2]
lambdaNHPP2= [35.22,0.224,-0.16]
A=4
B=6
Delta=0.01

i = tabulate(m(:,5))
j = i (:, 2)
n=j'
t = (1:length(n))/4
tn=[t;n]
ncum = cumsum(n)
p=A:Delta:B
HP = lambdaHPP*p
NHP1 = lambdaNHPP1(1)*p-lambdaNHPP1(2)/(2*pi)*cos(2*pi*(p+lambdaNHPP1(3)))+lambdaNHPP1(2)/(2*pi)*cos(2*pi*lambdaNHPP1(3))
NHP2 = lambdaNHPP2(1)*p+lambdaNHPP2(2)*(0.5*p-1/(8*pi)*sin(4*pi*(p+lambdaNHPP2(3)))+1/(8*pi)*sin(4*pi*lambdaNHPP2(3)))

plot(t,ncum,'b', 'LineWidth',3)
axis([4 6 139 220])
xlabel('Years')
ylabel('Aggregate number of losses / Mean value function')
hold on
plot(p,HP,'r*:')
hold on
plot(p,NHP1,'k*--')
hold on
plot(p,NHP2,'g*-.')