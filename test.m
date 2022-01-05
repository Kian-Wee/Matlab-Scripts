A = rand(10,5);
cm = colormap(parula(size(A,1)));                           % Default Colormap
% c = A(:,5);
for i=1:10
    x1=A(i,1);
    y1=A(i,2);
    x2=A(i,3);
    y2=A(i,4);
    plot([x1,x2],[y1,y2],'Color', cm(i,:))
    hold on
end