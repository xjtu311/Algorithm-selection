Solver : ( SATenstein-LS )
SLS随机局部搜索

Algorithm Configurator : (FocusedILS)  2010
	MSAC(Sequential Model-based Algorithm Configuration)  2015
	FocusedILS是产生于ParamILS框架
	对于一个给定的算法存在一个参数空间，搜索算法的参数空间中可能的参数配置，通过给定的基准实例集合运行算法，评价这个参数配置。
算法：SAT4J，SAPS，GLS+
创建算法的命令行，执行命令行，保存数据，从结果中选择性能最好的
算法文件为wrapper，参数文件为param.txt
	MSAC 与ParamILS 的区别
MSAC支持连续的参数, ParamILS只支持分类参数,MSAC不能运行ParamILS的所有的对象，所实例文件中的实例排序不重要，时间预算中比ParamILS更加精细，处理CPU时间还有其他的时钟信息等 。。。。。

Performance metric :  (PAR  score)

Portfolio builder :  (SATzilla)			   定义一个算法集合S，训练的数据集合I，D表示数据集I中每个实例被对应每个算法集合的算法运行的性能值，用F表示数据集I中每个实例对应的特征集合。
首先，根据局部搜索技术选择集合S中的至少两个求解器（算法）作为初始的求解器集合，并设置固定的时间预算为{0s，2s，5s，10s}，使用前向选择特征来减少特征集的规模，基于二次函数的前向选择，新特征集为F’。对于算法集合S中的每个求解器使用岭回归的变体去预测性能集合D与新特征集F’，构建一个F’中的特征与预测性能之间的映射关系。顺序的执行求解器集合，如果实例没有被解决，计算特征值，预测器评估每个算法的性能，执行预测性能最好的算法。
	
