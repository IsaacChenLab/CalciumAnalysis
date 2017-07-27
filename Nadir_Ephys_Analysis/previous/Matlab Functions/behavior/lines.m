

scrsz = get(0,'ScreenSize');
f1=figure('Position',scrsz,'color',[0 0 0]);
drawnow;
for a=-5:1:5
    x=[a a a a a a a a a a a];
    y=[0 1 2 3 4 5 6 7 8 9 10];
    hold on;
    plot(x,y,'w');
end
%axis square;
axis equal;
axis fill;
axis off;

scrsz = get(0,'ScreenSize');
f1=figure('Position',scrsz,'color',[0 0 0]);
drawnow;
m=0;
hold on;
for n=-5:5
for i=1:5
    x1(i)=i;
    y1(i)=m*x1(i)+n;
end
plot(y1,'w');
end
%axis square;
axis fill;
axis off;


H=line(m,n,style)
%ecuatia unei drepte este y=m*x+n
%functie cu parametri m,n,style
%pt m dat  (m=tg(alfa)) unghiul dreptei ;n=valoarea lui y(ordonata)
%cind x=0
%daca vreau mai rare n=-5:5 , mai dese n mai mare
%daca vreau o retea peste grid;
%cu acest program vom obtine paralelograme paralele umplute
%deci linii groase paralele (grosimea depinde de n)
scrsz = get(0,'ScreenSize');
f1=figure('Position',scrsz,'color',[0 0 0]);
drawnow;
m=0;
hold on;
for n=-11:2:11 %pas 2 ca sa fac un interval intre paralelograme.
for i=1:5
    x1(i)=i;
    y1(i)=m*x1(i)+n;
    %disp(x1);
    %disp(y1);
   % X(i)=x1(i);
   % Y(i)=y1(i);
end
%paralelogramul ABCD pe care-l umplu cu fill
    %A
    X(1)=x1(1);
    Y(1)=m*x1(1)+n-1;
    %B
    X(2)=x1(5);
    Y(2)=m*x1(5)+n-1;
    %C
    X(3)=x1(5);
    Y(3)=m*x1(5)+n;
    %D
    X(4)=x1(1);
    Y(4)=m*x1(1)+n;
    %disp(X);
    %disp(Y);
 H=plot(y1,'w');
%fill(X,Y,'g');
for j=1:4
   XX(j)=X(j);
   YY(j)=Y(j);
end
fill(XX,YY,'w');
end
%axis square;
axis fill;
axis off;
%pause(5);clf;


-------------------------------------------------

paralelogram ABCD (A=jos stinga; B=jos dreapta;C=sus dreapta;D=sus
stinga;)
Obtinem coordonatele celor 4 puncte XX,YY
ciudat!! cand am obtinut X,Y imi iese unul in plus de aceea am facut bucla de 4
XX=X YY=Y

%ecuatia unei drepte este y=m*x+n
%trebuie sa fac o functie cu parametrii m,n,style
%pt m dau (m=tg(alfa))deci unghiul dreptei ;n=valoarea lui y(ordonata)
%cind x=0
m=1;
%daca vreau mai rare n=-5:5 ,mai dese n mai mare
%daca vrei o retea peste pui grid;
hold on;

%for n=-1:1
for n=1
for i=1:5
    x1(i)=i;
    y1(i)=m*x1(i)+n;
    %disp(x1);
    %disp(y1);
   % X(i)=x1(i);
   % Y(i)=y1(i);
end
   %A
    X(1)=x1(1);
    Y(1)=m*x1(1)+n-1;
    %B
    X(2)=x1(5);
    Y(2)=m*x1(5)+n-1;
    %C
    X(3)=x1(5);
    Y(3)=m*x1(5)+n;
    %D
    X(4)=x1(1);
    Y(4)=m*x1(1)+n;
    disp(X);
    disp(Y);
    plot(y1,'r');
%fill(X,Y,'g');
end
for j=1:4
   XX(j)=X(j);
   YY(j)=Y(j);
end
fill(XX,YY,'g');
%
for n=3
for i=1:5
    x1(i)=i;
    y1(i)=m*x1(i)+n;
    %disp(x1);
    %disp(y1);
   % X(i)=x1(i);
   % Y(i)=y1(i);
end
    %A
    X(1)=x1(1);
    Y(1)=m*x1(1)+n-1;
    %B
    X(2)=x1(5);
    Y(2)=m*x1(5)+n-1;
    %C
    X(3)=x1(5);
    Y(3)=m*x1(5)+n;
    %D
    X(4)=x1(1);
    Y(4)=m*x1(1)+n;
    disp(X);
    disp(Y);
 plot(y1,'r');
%fill(X,Y,'g');
end
for j=1:4
   XX(j)=X(j);
   YY(j)=Y(j);
end
fill(XX,YY,'g');
%
for n=5
for i=1:5
    x1(i)=i;
 y1(i)=m*x1(i)+n;
    %disp(x1);
    %disp(y1);
   % X(i)=x1(i);
   % Y(i)=y1(i);
end
    %A
    X(1)=x1(1);
    Y(1)=m*x1(1)+n-1;
    %B
    X(2)=x1(5);
    Y(2)=m*x1(5)+n-1;
    %C
    X(3)=x1(5);
    Y(3)=m*x1(5)+n;
    %D
    X(4)=x1(1);

 Y(4)=m*x1(1)+n;
    disp(X);
    disp(Y);

    plot(y1,'r');
%fill(X,Y,'g');
end
for j=1:4
   XX(j)=X(j);
   YY(j)=Y(j);
end
fill(XX,YY,'g');
axis square;
%axis off;
%pause(5);clf;

