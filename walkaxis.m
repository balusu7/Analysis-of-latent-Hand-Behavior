x=[0:-2:-60];
y=(180+4*x)/-3;
Caxis = vertcat(x,y);
Caxis=Caxis';
plot( x,y,'g*');