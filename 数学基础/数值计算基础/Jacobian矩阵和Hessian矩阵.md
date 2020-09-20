# Jacobian矩阵和Hessian矩阵





在向量分析中, 雅可比矩阵是一阶偏导数以一定方式排列成的矩阵, 其行列式称为雅可比行列式. 还有, 在代数几何中, 代数曲线的雅可比量表示雅可比簇：伴随该曲线的一个代数群, 曲线可以嵌入其中. 它们全部都以数学家卡尔·雅可比(Carl Jacob, 1804年10月4日－1851年2月18日)命名；

## 一、雅可比矩阵

雅可比矩阵的重要性在于它体现了一个可微方程与给出点的最优线性逼近. 因此, 雅可比矩阵类似于多元函数的导数.

假设 $${R_n} \to {R_m}$$ 是一个从欧式n维空间转换到欧式m维空间的函数. 
$$
\left\{\begin{array}{l}
y_{1}=f_{1}\left(x_{1}, \ldots, x_{n}\right) \\
y_{2}=f_{2}\left(x_{1}, \ldots, x_{n}\right) \\
\cdots \\
y_{m}=f_{n}\left(x_{1}, \ldots, x_{n}\right)
\end{array}\right.
$$


这些函数的偏导数(如果存在)可以组成一个m行n列的矩阵, 这就是所谓的雅可比矩阵：

$$
\mathbf{J}=\left[\begin{array}{ccc}
\frac{\partial \mathbf{f}}{\partial x_{1}} & \cdots & \frac{\partial \mathbf{f}}{\partial x_{n}}
\end{array}\right]=\left[\begin{array}{ccc}
\frac{\partial f_{1}}{\partial x_{1}} & \cdots & \frac{\partial f_{1}}{\partial x_{n}} \\
\vdots & \ddots & \vdots \\
\frac{\partial f_{m}}{\partial x_{1}} & \cdots & \frac{\partial f_{m}}{\partial x_{n}}
\end{array}\right]
$$


此矩阵表示为:$${J_F}({x_1}, \ldots ,{x_n})$$ 

如果  $$p$$ 是 $${R_n}$$ 中的一点, $$F$$ 在 $${p}$$点可微分, 那么在这一点的导数由 $${J_F}({p})$$ 给出(这是求该点导数最简便的方法). 在此情况下, 由 $$F(p)$$ 描述的线性算子即接近点 $$p$$ 的 $$F$$ 的最优线性逼近, $${{x}}$$ 逼近于$${{p}}$$:

<center><img src=".\img\image-20200717000649311.png" alt="image-20200717000649311"  /></center>

与泰勒一阶展开近似

## 二、雅可比行列式

如果m = n, 那么 $$F$$ 是从 n 维空间到 n 维空间的函数, 且它的雅可比矩阵是一个方块矩阵. 于是我们可以取它的行列式, 称为雅可比行列式.

在某个给定点的雅可比行列式提供了 在接近该点时的表现的重要信息.

- 如果连续可微函数 $$F$$ 在 $$p$$ 点的雅可比行列式不是零, 那么它在该点附近具有反函数. 这称为反函数定理. 更进一步,
- 如果 $$p$$ 点的雅可比行列式是正数, 则 $$F$$ 在 $$p$$ 点的取向不变；
- 如果是负数, 则 $$F$$ 的取向相反

而从雅可比行列式的绝对值, 就可以知道函数 $$F$$ 在 $$p$$ 点的缩放因子；这就是为什么它出现在换元积分法中.

对于取向问题可以这么理解, 例如一个物体在平面上匀速运动, 如果施加一个正方向的力 $$F$$ , 即取向相同, 则加速运动, 类比于速度的导数加速度为正；如果施加一个反方向的力 $$F$$ , 即取向相反, 则减速运动, 类比于速度的导数加速度为负.



## 三、海森Hessian矩阵

1. 二阶导数 $$f^{\prime \prime}(x)$$  刻画了曲率。假设有一个二次函数（实际任务中，很多函数不是二次的，但是在局部可以近似为二次函数）：

   - 如果函数的二阶导数为零，则它是一条直线。如果梯度为 1，则当沿着负梯度的步长为 $$\epsilon$$  时，函数值减少$$\epsilon$$  。

   - 如果函数的二阶导数为负，则函数向下弯曲。如果梯度为1，则当沿着负梯度的步长为 $$\epsilon$$  时，函数值减少的量大于 $$\epsilon$$  。

   - 如果函数的二阶导数为正，则函数向上弯曲。如果梯度为1，则当沿着负梯度的步长为 $$\epsilon$$  时，函数值减少的量少于 $$\epsilon$$  。

     <img src=".\img\20a66edb9b4c2272719593a5566e18d2.png" alt="20a66edb9b4c2272719593a5566e18d2" style="zoom:67%;" />

2. 在数学中, 海森矩阵(Hessian matrix或Hessian)是一个自变量为向量的实值函数的<u>二阶偏导数</u>组成的方块矩阵, 此函数如下：
   $$
   f({x_1},{x_2} \ldots ,{x_n})
   $$
   如果 $$f$$ 的所有二阶导数都存在, 那么 $$f$$ 的海森矩阵即：
   $$
   H{(f)_{ij}}(x) = {D_i}{D_j}f(x)
   $$
   其中 $$x = ({x_1},{x_2} \ldots ,{x_n})$$  即 $$H(f)$$ 为: 
   $$
   \begin{bmatrix} 
   \frac{\partial^2 f}{\partial x_1^2} & \frac{\partial^2 f}{\partial x_1\,\partial x_2} & \cdots & \frac{\partial^2 f}{\partial x_1\,\partial x_n} \\  \\ 
   \frac{\partial^2 f}{\partial x_2\,\partial x_1} & \frac{\partial^2 f}{\partial x_2^2} & \cdots & \frac{\partial^2 f}{\partial x_2\,\partial x_n} \\  \\ 
   \vdots & \vdots & \ddots & \vdots \\  \\ 
   \frac{\partial^2 f}{\partial x_n\,\partial x_1} & \frac{\partial^2 f}{\partial x_n\,\partial x_2} & \cdots & \frac{\partial^2 f}{\partial x_n^2} 
   \end{bmatrix}
   $$
   即海森矩阵的第 $$i$$ 行 $$j$$ 列元素为：$$\mathbf{H}_{i, j}=\frac{\partial^{2}}{\partial x_{i} \partial x_{j}} f(\vec{\mathbf{x}})$$ 

   海森矩阵被应用于牛顿法解决的大规模优化问题.

3. 当二阶偏导是连续时，海森矩阵是对称阵，即有： $$\mathbf{H}=\mathbf{H}^{T}$$

4. 对于特定方向 $$\vec{\mathrm{d}}$$ 上的二阶导数为：$$\vec{\mathrm{d}}^{T} \mathbf{H} \vec{\mathrm{d}}$$ 。

   - 如果  $$\vec{\mathrm{d}}$$  是海森矩阵的特征向量，则该方向的二阶导数就是对应的特征值。
   - 如果  $$\vec{\mathrm{d}}$$  不是海森矩阵的特征向量，则该方向的二阶导数就是所有特征值的加权平均，权重在 `(0,1)`之间。且与  $$\vec{\mathrm{d}}$$  夹角越小的特征向量对应的特征值具有更大的权重。
   - 最大特征值确定了最大二阶导数，最小特征值确定最小二阶导数。





## 四、海森矩阵与学习率

1. 将 $$f(\vec{\mathbf{x}})$$ 在 $$\vec{\mathbf{x_0}}$$ 处泰勒展开：$$f(\vec{\mathbf{x}}) \approx f\left(\vec{\mathbf{x}}_{0}\right)+\left(\vec{\mathbf{x}}-\vec{\mathbf{x}}_{0}\right)^{T} \vec{\mathbf{g}}+\frac{1}{2}\left(\vec{\mathbf{x}}-\vec{\mathbf{x}}_{0}\right)^{T} \mathbf{H}\left(\vec{\mathbf{x}}-\vec{\mathbf{x}}_{0}\right)$$  。其中：  $$\vec{\mathbf{g}}$$ 为  $$\vec{\mathbf{x_0}}$$ 处的梯度 ；$$H$$ 为  $$\vec{\mathbf{x_0}}$$ 处的海森矩阵。

   根据梯度下降法：$$\vec{\mathbf{x}}^{\prime}=\vec{\mathbf{x}}-\epsilon \nabla_{\vec{\mathbf{x}}} f(\vec{\mathbf{x}})$$ 

   应用在点  $$\vec{\mathbf{x_0}}$$  ，有：$$f\left(\vec{\mathbf{x}}_{0}-\epsilon \vec{\mathbf{g}}\right) \approx f\left(\vec{\mathbf{x}}_{0}\right)-\epsilon \vec{\mathbf{g}}^{T} \vec{\mathbf{g}}+\frac{1}{2} \epsilon^{2} \vec{\mathbf{g}}^{T} \mathbf{H} \vec{\mathbf{g}}$$ 

   - 第一项代表函数在点   $$\vec{\mathbf{x_0}}$$   处的值。
   - 第二项代表由于斜率的存在，导致函数值的变化。
   - 第三项代表由于曲率的存在，对于函数值变化的矫正。

2. 注意：如果 $$\frac{1}{2} \epsilon^{2} \vec{\mathbf{g}}^{T} \mathbf{H} \vec{\mathbf{g}}$$ 较大，则很有可能导致：沿着负梯度的方向，函数值反而增加！

   - 如果  $$\vec{\mathbf{g}}^{T} \mathbf{H} \vec{\mathbf{g}} \leq 0$$，则无论  $$\epsilon$$ 取多大的值， 可以保证函数值是减小的。

   - 如果  $$\vec{\mathbf{g}}^{T} \mathbf{H} \vec{\mathbf{g}} > 0$$ ，则学习率   $$\epsilon$$ 不能太大。若   $$\epsilon$$  太大则函数值增加。

     - 根据 $$f\left(\vec{\mathbf{x}}_{0}-\epsilon \vec{\mathbf{g}}\right)-f\left(\vec{\mathbf{x}}_{0}\right)<0$$ ，则需要满足：$$\epsilon<\frac{2 \vec{\mathrm{g}}^{T} \vec{\mathrm{g}}}{\vec{\mathrm{g}}^{T} \mathrm{H} \vec{\mathrm{g}}}$$  。若  $$\epsilon \geq \frac{2 \vec{\mathrm{g}}^{T} \vec{\mathrm{g}}}{\vec{\mathrm{g}}^{T} \mathrm{H} \vec{\mathrm{g}}}$$ ，则会导致沿着负梯度的方向函数值在增加。

     - 考虑最速下降法，选择使得 $$f$$ 下降最快的    $$\epsilon$$  ，则有：$$\epsilon^{*}=\arg \min _{\epsilon, \epsilon>0} f\left(\vec{\mathbf{x}}_{0}-\epsilon \vec{\mathbf{g}}\right)$$ 。求解 $$\frac{\partial}{\partial \epsilon} f\left(\vec{\mathbf{x}}_{0}-\epsilon \vec{\mathbf{g}}\right)=0$$  有 ：$$\epsilon^{*}=\frac{\vec{\mathbf{g}}^{T} \vec{\mathbf{g}}}{\vec{\mathbf{g}}^{T} \mathbf{H} \vec{\mathbf{g}}}$$ 。

       根据  $$\vec{\mathbf{g}}^{T} \mathbf{H} \vec{\mathbf{g}} > 0$$ ，很明显有： ：$$\epsilon<\frac{2 \vec{\mathrm{g}}^{T} \vec{\mathrm{g}}}{\vec{\mathrm{g}}^{T} \mathrm{H} \vec{\mathrm{g}}}$$ 。

3. 由于海森矩阵为实对称阵，因此它可以进行特征值分解。假设其特征值从大到小排列为：$$\lambda_{1} \geq \lambda_{2} \geq \cdots \geq \lambda_{n}$$ 。

   海森矩阵的瑞利商为：$$R(\vec{\mathbf{x}})=\frac{\vec{\mathbf{x}}^{T} \mathbf{H} \vec{\mathbf{x}}}{\vec{\mathbf{x}}^{T} \vec{\mathbf{x}}}, \vec{\mathbf{x}} \neq \vec{\mathbf{0}}$$ 。可以证明：

   $$
   \begin{aligned}
   \lambda_{n} & \leq R(\vec{\mathbf{x}}) \leq \lambda_{1} \\
   \lambda_{1} &=\max _{\vec{\mathbf{x}} \neq \vec{\mathbf{0}}} R(\vec{\mathbf{x}}) \\
   \lambda_{n} &=\min _{\vec{\mathbf{x}} \neq \vec{0}} R(\vec{\mathbf{x}})
   \end{aligned}
   $$


   根据 $$\epsilon^{*}=\frac{\vec{\mathbf{g}}^{T} \vec{\mathbf{g}}}{\vec{\mathbf{g}}^{T} \mathbf{H} \vec{\mathbf{g}}}=\frac{1}{R(\vec{\mathbf{g}})}$$ 可知：海森矩阵决定了学习率的取值范围。最坏的情况下，梯度  $$\vec{\mathbf{g}}$$与海森矩阵最大特征值 $$\lambda_1$$ 对应的特征向量平行，则此时最优学习率为  $$\frac{1}{\lambda_{1}}$$ 。



## 五、驻点与全局极小点

1. 满足导数为零的点（即 $$f^{\prime \prime}(x)=0$$）称作驻点。驻点可能为下面三种类型之一：

   - 局部极小点：在  $$x$$ 的一个邻域内，该点的值最小。
   - 局部极大点：在 $$x$$  的一个邻域内，该点的值最大。
   - 鞍点：既不是局部极小，也不是局部极大。

   <img src=".\img\critical_point.png" alt="critical_point" style="zoom:80%;" />

2. 全局极小点： $$x^{*}=\arg \min _{x} f(x)$$ 

   - 全局极小点可能有一个或者多个。

   - 在深度学习中，目标函数很可能具有非常多的局部极小点，以及许多位于平坦区域的鞍点。这使得优化非常不利。

     因此通常选取一个非常低的目标函数值，而不一定要是全局最小值。

     <img src=".\img\deeplearning_optimization.png" alt="deeplearning_optimization" style="zoom:80%;" />

3. 二阶导数可以配合一阶导数来决定驻点的类型：

   - 局部极小点： $$f^{\prime}(x)=0, f^{\prime \prime}(x)>0$$ 。
   - 局部极大点： $$f^{\prime}(x)=0, f^{\prime \prime}(x)<0$$。
   - $$f^{\prime}(x)=0, f^{\prime \prime}(x)=0$$ ：驻点的类型可能为任意三者之一。

4. 对于多维的情况类似有：

   - 局部极小点：$$\nabla_{\mathbf{x}}f(\mathbf{x})=0$$ ，且海森矩阵为**正定**的（即所有的特征值都是正的）,当海森矩阵为正定时，任意方向的二阶偏导数都是正的。

   - 局部极大点：$$\nabla_{\mathbf{x}} f(\mathbf{x})=0$$，且海森矩阵为**负定**的（即所有的特征值都是负的）。

     当海森矩阵为负定时，任意方向的二阶偏导数都是负的。

   - $$\nabla_{\mathbf{x}} f(\mathbf{x})=0$$ ，且海森矩阵的特征值中至少一个正值、至少一个负值时，为鞍点。

   - 当海森矩阵非上述情况时，驻点类型无法判断。

   下图为 $$f({\mathbf{x}})=x_{1}^{2}-x_{2}^{2}$$ 在原点附近的等值线。其海森矩阵为一正一负。

   - 沿着  $$x_1$$ 方向，曲线向上弯曲；沿着  $$x_2$$ 方向，曲线向下弯曲。
   - 鞍点就是在一个横截面内的局部极小值，另一个横截面内的局部极大值。

<img src=".\img\saddle.png" alt="saddle" style="zoom:80%;" />









## 参考

- <a href="http://www.huaxiaozhuan.com/数学基础/chapters/3_numerical_computaion.html" target="">数值计算:三、二阶导数与海森矩阵</a> 
- <a href="http://jacoxu.com/jacobian矩阵和hessian矩阵/" target="">Jacobian矩阵和Hessian矩阵</a> 