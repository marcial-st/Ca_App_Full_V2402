%%%%%% Thesis
%%%% MST,JMRC,RBR
%%%% Kalman estimator function  Predictor/Corrector
function [xk,Pk]=kalman_v1(A,H,Q,R,zk,x,P,flag)
[m,n]=size(Q);
I=eye(m,n);
    % Only Prediction (when no cell is found)
    xk=A*x;
    Pk=A*P*A'+Q;
if 1==flag
    % Correction
    Kk=Pk*H'/(H*Pk*H'+R);
    xk=xk+Kk*(zk-H*xk);
    Pk=(I-Kk*H)*Pk;
end