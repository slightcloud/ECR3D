%****************************************************************
%该子函数用高斯-赛德尔加超松弛迭代的方法求解三维空间中的电势%
%****************************************************************

function phi = GS_SOR_3D(rho,nx,ny,nz,tol,w,choice)
%***************************************************************
%rho代表的是每个格点上的电荷密度；
%nx代表的是系统x方向的格点总数；
%ny代表的是系统y方向的格点总数；
%nz代表的是系统z方向的格点总数；
%tol代表的是迭代结果的容差；
%w代表的是松弛因子
%choice的取值为1和2，
%          choice=1：表示两个边界的电势为0
%          choice=2：表示两个边界的电势对坐标的一阶导数为0
%***************************************************************
global dg EPS0    %声明全局变量
solver_it = 20000;  %迭代步数
phi = zeros(nx,ny,nz);  %给电势赋初值
%**********************choice = 1*******************************
if choice == 1
    for i = 1:solver_it
        g = 1/6*((rho(2:nx-1,2:ny-1,2:nz-1)/EPS0)*dg*dg+phi(3:nx,2:ny-1,2:nz-1)+phi(1:nx-2,2:ny-1,2:nz-1)+...
            phi(2:nx-1,3:ny,2:nz-1)+phi(2:nx-1,1:ny-2,2:nz-1)+phi(2:nx-1,2:ny-1,3:nz)+phi(2:nx-1,2:ny-1,1:nz-2));   %只求解2到nz-1的电势值，依据是泊松方程
        phi(2:nx-1,2:ny-1,2:nz-1) = phi(2:nx-1,2:ny-1,2:nz-1)+w*(g-phi(2:nx-1,2:ny-1,2:nz-1));                %做超松弛迭代
        if mod(i,25) == 0
            res = rho(2:nx-1,2:ny-1,2:nz-1)/EPS0+(phi(3:nx,2:ny-1,2:nz-1)+phi(1:nx-2,2:ny-1,2:nz-1)+...
            phi(2:nx-1,3:ny,2:nz-1)+phi(2:nx-1,1:ny-2,2:nz-1)+phi(2:nx-1,2:ny-1,3:nz)+phi(2:nx-1,2:ny-1,1:nz-2)-6*phi(2:nx-1,2:ny-1,2:nz-1))/(dg*dg); %求出迭代的残差
            sum_res = sum(sum(sum(res.^2)));          %求出所有残差的和
            judge = sqrt(sum_res/(nx*ny*nz));       %求出平均残差
            fprintf('It is the %d interp\tThe residual is %g\n',i,judge);
            if judge < tol                  %将平均残差与设定容差进行对比，若满足残差要求则停止循环
                break;
            end
        end
%         fprintf('It is the %d interp\n',i);
    end
    clc;
    if i == solver_it
        fprintf('The interation is not converged!!!\n');
    end
end

%**********************choice = 2*******************************
if choice == 2
    for i = 1:solver_it
        g = 1/6*((rho(2:nx-1,2:ny-1,2:nz-1)/EPS0)*dg*dg+phi(3:nx,2:ny-1,2:nz-1)+phi(1:nx-2,2:ny-1,2:nz-1)+...
            phi(2:nx-1,3:ny,2:nz-1)+phi(2:nx-1,1:ny-2,2:nz-1)+phi(2:nx-1,2:ny-1,3:nz)+phi(2:nx-1,2:ny-1,1:nz-2));   %只求解2到nz-1的电势值，依据是泊松方程
        phi(2:nx-1,2:ny-1,2:nz-1) = phi(2:nx-1,2:ny-1,2:nz-1)+w*(g-phi(2:nx-1,2:ny-1,2:nz-1));                %做超松弛迭代
        phi(1,:,:) = phi(2,:,:);                                              %左边界的电势
        phi(nx,:,:) = phi(nx-1,:,:);                                          %右边界的电势
        phi(:,1,:) = phi(:,2,:);
        phi(:,ny,:) = phi(:,ny-1,:);
        phi(:,:,1) = phi(:,:,2);
        phi(:,:,nz) = phi(:,:,nz-1);
        if mod(i,25) == 0
            res = rho(2:nx-1,2:ny-1,2:nz-1)/EPS0+(phi(3:nx,2:ny-1,2:nz-1)+phi(1:nx-2,2:ny-1,2:nz-1)+...
            phi(2:nx-1,3:ny,2:nz-1)+phi(2:nx-1,1:ny-2,2:nz-1)+phi(2:nx-1,2:ny-1,3:nz)+phi(2:nx-1,2:ny-1,1:nz-2))/(dg*dg); %求出迭代的残差
            sum_res = sum(sum(sum(res.^2)));          %求出所有残差的和
            judge = sqrt(sum_res/(nx*ny*nz));       %求出平均残差
            if judge < tol                  %将平均残差与设定容差进行对比，若满足残差要求则停止循环
                break;
            end
        end
    end
    if i == solver_it
        fprintf('The interation is not converged!!!\n');
    end
end
        

