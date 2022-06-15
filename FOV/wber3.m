   function PP = wber3(A,L,col)          
       
% Computes the boundary of the field of values W(A) of A for arbitrary 
%     complex matrices A .  W(A) := { x' * A * x / x' * x } 
% Input:  A is a square matrix 
%         L is the number of points plotted on the curve(s)    
%                  (use L = 40, 60, 100 for example)
%         col = color, 'k' for black, 'r' for red (default)
% Output: PP: complex numbers on the boundary of the FOV W(A) and
%         VV: corresponding unit vectors that generate the PP points

if nargin == 1, L = 80; col = 'r'; end 
if nargin == 2, col = 'r'; end
[n n] = size(A); Vt = zeros(n,1); PP = NaN * ones(2*(L+1),1);    % definitions
PP = NaN * ones(2*(L+1),1); opts.disp = 0; opts.maxit = 1600; warning off
H = (A+A')/2; K = (A-A')/(2i);                         % hermitean and skewherm parts
th = [0:pi/L:pi]; ko = cos(th); si = sin(th);           % setting up angular values
for k = 1:L+1                                           % rotated evalue computations
    k
  Ath = ko(k) * H + si(k) * K; 
  if n < 170 
    [V D] = eig(Ath);    
    Vt = V(:,n); PP(k) = Vt'*A*Vt;                      % desired evalues/evectors
    Vt = V(:,1); PP(L+1+k) = Vt'*A*Vt;                  % (use max and min evalue/vector)
  else 
    [V D] = eigs(Ath,1,'Lr',opts); PP(k) = V'*A*V;
    [V D] = eigs(Ath,1,'Sr',opts); PP(L+1+k) = V'*A*V;
  end 
end,

%plot(real(PP),imag(PP),col); %hold on, % plot(0,0,'*'), % grid on
%set(gca,'PlotBoxAspectRatio',[1,1,1]);  set(gca,'DataAspectRatio',[1,1,1]);
% title(['Numerical range of a ',num2str(n),' by ',num2str(n),' matrix A'],...
%        'FontSize',14), 
%xlabel('real  axis', 'FontSize',15), ylabel('imaginary  axis', 'FontSize',15), 
%rectangle('Position',[x1 y1 x2-x1 y4-y1]); hold off;   

% Ar = H; Ai =  K;
% Atest = H; [V D] = eig(Atest);    
%     Vt = V(:,1); P1= Vt'*A*Vt; P1 = [real(P1),imag(P1)];
%     plot(P1(1),P1(2),'or'), plot([0,45],[0,0],'r'), plot(P1(1),0,'or')
%     mult = 2;  plot([0,mult*P1(1)],[0,P1(2)*mult],'--r'), plot([P1(1),P1(1)],[-3,2*P1(2)],':r')
% ra = [.6,(1-.36)^.5];
% Atest = ra(1)*H + ra(2)*K;   [V D] = eig(Atest);    
%      Vt = V(:,1); P2= Vt'*A*Vt; P2 = [real(P2),imag(P2)];
%      plot(P2(1),P2(2),'ob'), mult = 1.5; plot([0,mult*P2(1)],[0,mult*P2(2)],'--b') 
%      mult = 35; plot([0,mult*ra(1)],[0,mult*ra(2)],'b')
%      rProj = D(1)*ra; plot(rProj(1),rProj(2),'ob'), 
%      mult = 1; plot([rProj(1),2*P2(1)-rProj(1)],[rProj(2),2*P2(2)-rProj(2)],':b')
%      text(20.5,-.7,'P','FontSize',13), text(-1,-1,'0','FontSize',13), 
%      text(2,10,'Proj(P)','FontSize',13), 



