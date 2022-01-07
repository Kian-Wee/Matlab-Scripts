% A = rand(10,5);
% cm = colormap(parula(size(A,1)));                           % Default Colormap
% % c = A(:,5);
% for i=1:10
%     x1=A(i,1);
%     y1=A(i,2);
%     x2=A(i,3);
%     y2=A(i,4);
%     plot([x1,x2],[y1,y2],'Color', cm(i,:))
%     hold on
% end

% %% Break loop if keypress to save to excel
% DlgH = figure;
% H = uicontrol('Style', 'PushButton', ...
%                     'String', 'Break', ...
%                     'Callback', 'delete(gcbf)');
% 
% %% Display data
% while (ishandle(H))
%     disp(cputime);
% end
% disp("end")

ButtonHandle = uicontrol('Style', 'PushButton', ...
                         'String', 'Stop loop', ...
                         'Callback', 'delete(gcbf)');
for k = 1:1e6
  disp(k)
  if ~ishandle(ButtonHandle)
    disp('Loop stopped by user');
    break;
  end
pause(0.01); % A NEW LINE
end
disp("end")

%https://www.mathworks.com/matlabcentral/answers/285872-shading-with-plot3#answer_223465
r = sqrt(x.^2 + y.^2 + z.^2);
g = patch('Vertices', [x(:), y(:),z(:); nan nan nan], 'Faces', (1:length(x)+1).', 'FaceVertexCData', [r(:); nan], 'EdgeColor', 'interp', 'Marker','*')