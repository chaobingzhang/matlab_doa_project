%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%日期：2021/11/17
%函数名称：通过极大值来判断方位
%参数作用：
%       1、accuracy          目标方位点后测量精度，以目标为中心，前后20°
%       2、deltaSigma        目标估计小数点精度 
%       3、Nsigma            待测量目标矩阵长度
%       4、X                 带测量矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ret]=Extremum(accuracy,deltaSigma,Nsigma,X)
    %创建待处理数组
    diff_value_accuracy=accuracy*2/deltaSigma;% 计算差值
    half_diff_value_accuracy=diff_value_accuracy/2;
    Dispose_X = zeros(1,Nsigma+diff_value_accuracy);
    for i=1:half_diff_value_accuracy
        Dispose_X(i)=X(Nsigma-half_diff_value_accuracy+i);
    end
    for i=half_diff_value_accuracy+1:Nsigma+half_diff_value_accuracy
        Dispose_X(i)=X(i-half_diff_value_accuracy);
    end
    
    for i=Nsigma+half_diff_value_accuracy+1:Nsigma+diff_value_accuracy
        Index=i-Nsigma-half_diff_value_accuracy;
        Dispose_X(i)=X(Index);
    end

    
    Aim = zeros(1,Nsigma);%一行表示方位
    for i=half_diff_value_accuracy+1:Nsigma+half_diff_value_accuracy
        flag=1;
        for j=1:half_diff_value_accuracy
            if Dispose_X(i) <= Dispose_X(i+j)||Dispose_X(i) <= Dispose_X(i-j)
                flag=0;
            end
        end
        if flag==1
            Aim(i-half_diff_value_accuracy)=Dispose_X(i);
        end
    end
    
    [~,Ret]=sort(Aim,'descend');
    
%     %寻找两目标个方位值
%     flag=1;%判断标志位
%     max_num=0;%记录两个角度的极大值
%     max_i=0; %声能流最强目标
%     sec_num=0;
%     sec_i=0; %声能流第二强目标
%     for i=half_diff_value_accuracy+1:Nsigma+half_diff_value_accuracy
%         flag=1;
%         for j=1:half_diff_value_accuracy
%             if Dispose_X(i) <= Dispose_X(i+j)||Dispose_X(i) <= Dispose_X(i-j)
%                 flag=0;
%             end
%         end
%         if flag==1
%             if Dispose_X(i)>max_num
%                 sec_num=max_num;
%                 sec_i=max_i;
%                 max_num=Dispose_X(i);
%                 max_i=i-half_diff_value_accuracy;
%             elseif Dispose_X(i)>sec_num && Dispose_X(i)~=max_num
%                 sec_num=Dispose_X(i);
%                 sec_i=i-half_diff_value_accuracy;
%             end
%         end
%     end
end