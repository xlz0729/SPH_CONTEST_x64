/*
FLUIDS v.3 - SPH Fluid Simulator for CPU and GPU
Copyright (C) 2012. Rama Hoetzlein, http://fluids3.com

Fluids-ZLib license (* see part 1 below)
This software is provided 'as-is', without any express or implied
warranty.  In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
claim that you wrote the original software. Acknowledgement of the
original author is required if you publish this in a paper, or use it
in a product. (See fluids3.com for details)
2. Altered source versions must be plainly marked as such, and must not be
misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/


#ifndef DEF_FLUID_SYS
#define DEF_FLUID_SYS

// 标准库头文件
#include <iostream>
#include <string>

// 公共开源库头文件
#include "common\vector.h"
#include "common\xml_settings.h"
#include "common\gl_helper.h"
#include "common\camera3d.h"


// 宏定义
#define MAX_PARAM				60
#define GRID_UCHAR				0xFF
#define GRID_UNDEF				0xFFFFFFFF // 设定0xFFFFFFFF为未定义值

#define RUN_CPU_SPH				0
#define RUN_CUDA_INDEX_SPH		1	
#define RUN_CUDA_FULL_SPH		2
#define RUN_CPU_PCISPH			3
#define RUN_CUDA_INDEX_PCISPH	4 
#define RUN_CUDA_FULL_PCISPH	5

// 标量参数
#define PRUN_MODE				0	
#define PMAXNUM					1	// 最大粒子数量
#define PEXAMPLE				2	
#define PSIMSIZE				3	  
#define PSIMSCALE				4	// 缩放比例1:2
#define PGRID_DENSITY			5	// grid的密度
#define PGRIDSIZEREALSCALE		6	// grid的真实缩放比例
#define PVISC					7	// 体积
#define PRESTDENSITY			8	// 密度
#define PMASS					9	// 质量
#define PCOLLISIONRADIUS		10  // 碰撞半径
#define PSPACINGREALWORLD		11	// 真实世界的间隔
#define PSMOOTHRADIUS			12  // 光滑半径
#define PGASCONSTANT			13	// 气体含量
#define PBOUNDARYSTIFF			14	// 边界弹性
#define PBOUNDARYDAMP			15	// 边界阻力
#define PACCEL_LIMIT			16  // 加速度上限
#define PVEL_LIMIT				17	// 速度上限
#define PSPACINGGRAPHICSWORLD	18	// 图形世界的间隔
#define PGROUND_SLOPE			19	// 地面倾斜度
#define PFORCE_MIN				20	// 最小力
#define PFORCE_MAX				21	// 最大力
#define PMAX_FRAC				22	
#define PDRAWMODE				23	// 绘制方式
#define PDRAWSIZE				24	
#define PDRAWTEXT				26	// 绘制内容
#define PCLR_MODE				27
#define PPOINT_GRAV_AMT			28	
#define PSTAT_OCCUPANCY			29	// 有粒子的grid的数量
#define PSTAT_GRIDCOUNT			30	// grid中粒子的总数
#define PSTAT_NEIGHCNT			31	
#define PSTAT_NEIGHCNTMAX		32	// 领域内最大粒子数
#define PSTAT_SEARCHCNT			33	
#define PSTAT_SEARCHCNTMAX		34	// 搜索范围内最大粒子数
#define PSTAT_PMEM				35	
#define PSTAT_GMEM				36	
#define PTIME_INSERT			37	
#define PTIME_SORT				38
#define PTIME_COUNT				39
#define PTIME_PRESS				40
#define PTIME_FORCE				41
#define PTIME_ADVANCE			42
#define PTIME_RECORD			43	
#define PTIME_RENDER			44	
#define PTIME_TOGPU				45	
#define PTIME_FROMGPU			46	
#define PFORCE_FREQ				47	// 受力频率
#define PTIME_OTHER_FORCE		48
#define PTIME_PCI_STEP			49
#define	PDENSITYERRORFACTOR		50
#define PMINLOOPPCISPH			51  // PCI最小循环次数
#define PMAXLOOPPCISPH			52  // PCI最大循环次数
#define PMAXDENSITYERRORALLOWED 53  // 允许的最大密度误差
#define PKERNELSELF				54
#define PINITIALIZEDENSITY		55

// 向量参数
#define PGRIDVOLUMEMIN			0  // grid中编号最小的cell对应的最小的点
#define PGRIDVOLUMEMAX			1	
#define PBOUNDARYMIN			2	
#define PBOUNDARYMAX			3	
#define PINITPARTICLEMIN		4	
#define PINITPARTICLEMAX		5	
#define PEMIT_POS				6
#define PEMIT_ANG				7
#define PEMIT_DANG				8	
#define PEMIT_SPREAD			9
#define PEMIT_RATE				10
#define PPOINT_GRAV_POS			11	
#define PPLANE_GRAV_DIR			12	

// 布尔参数
#define PPAUSE					0
#define PDEBUG					1	
#define PUSE_CUDA				2	
#define	PUSE_GRID				3	
#define PWRAP_X					4	
#define PWALL_BARRIER			5	
#define PLEVY_BARRIER			6	
#define PDRAIN_BARRIER			7	
#define PPLANE_GRAV_ON			11	
#define PPROFILE				12
#define PCAPTURE				13	
#define PDRAWGRIDCELLS			14  
#define PPRINTDEBUGGINGINFO		15
#define PDRAWDOMAIN				16
#define	PDRAWGRIDBOUND			17
#define PUSELOADEDSCENE			18


// 全局变量
const int max_num_adj_grid_cells_cpu = 27;


// 在需要高频率（如渲染循环中）访问数据的时候，一般情况下SOA的效率高于AOS，
// 因为将需要频繁访问的数据连续存放会大大提高访问速度。
// 虽然AOS的结构可能更适合面向对象设计，但是在高度依赖效率的地方应该使用SOA。
// 
// 使用SOA 代替 AOS 结构体Particle只是用来计算粒子所占用内存大小的辅助数据结构
struct Particle {

	// offset - TOTAL: 120 (must be multiple of 12 = sizeof(Vector3DF) )
	Vector3DF		fpos;									// 0
	Vector3DF		fpredicted_pos;							// 12
	Vector3DF		fvel;									// 24
	Vector3DF		fvel_eval;								// 36
	Vector3DF		fcorrection_pressure_force;				// 48
	Vector3DF		fforce;									// 60
	Vector3DF		fsum_gradw;								// 72
	float			fsum_gradw_dot;							// 84
	float			fpressure;								// 88
	float			fcorrection_pressure;					// 92
	float			fdensity_reciprocal;					// 96
	float			fdensityError;							// 100
	int				fparticle_grid_cell_index;				// 104
	int				fnext_particle_index_in_the_same_cell;	// 108			
	DWORD			fclr;									// 112
	int				fpadding1;								// 116
	int				fpadding2;								// 120	填充字节 成员对齐，参见：《高质量程序设计指南C/C++语言》（第三版）P147， 8.1.4 成员对齐
};


// 函数命名规范:
// set和get函数开头第一个字母小写，函数名中每个单词开头字母大写
//     其余函数开头第一个字母大写，函数名中每个单词开头字母大写
class ParticleSystem {

public:
	ParticleSystem();
	~ParticleSystem();

	// set函数
	void        setup(bool bStart);
	void        setRender();
	inline void setToggle(int p);

	// get函数
	std::string      getModeStr();
	inline bool      getToggle(int p);
	inline int       getNumPoints();
	inline Vector3DF getGridRes();
	inline int       getSelected();
	inline int       getGridTotal();
	inline int       getGridAdjCnt();
	inline float     getParam(int p);
	inline Vector3DF getVec(int p);
	inline float     getDT();

	// 运行函数
	void Run();
	void RunCPUSPH();
	void RunCPUPCISPH();

	// 绘制函数
	inline void DrawParticleInfo();
	void        Draw(Camera3D& cam, float rad);

private:
	static int const lutSize = 100000;

	float lutKernelM4[lutSize];
	float lutKernelPressureGrad[lutSize];

	std::string	scene_name_;

	// 模拟参数
	float     param_[MAX_PARAM];
	Vector3DF vec_[MAX_PARAM];
	bool      toggle_[MAX_PARAM];

	// XML 设置文件
	XmlSettings	xml;

	// 绘制函数相关变量
	int    selected_;
	int	   vbo_[3];
	DWORD* clr_;
	Image  image_;

	// SPH光滑核函数系数
	float poly6_kern_;
	float lap_kern_;
	float spiky_kern_;

	// 时间步长
	float time_;
	float time_step_;
	float time_step_sph_;
	float time_step_wcsph_;
	float time_step_pcisph_;
	int	  frame_;

	// 粒子相关属性
	int	       num_points_;
	Vector3DF* pos_;
	Vector3DF* vel_;
	Vector3DF* vel_eval_;
	Vector3DF* force_;
	Vector3DF* correction_pressure_force_;
	float*     pressure_;
	float*     density_;
	uint*      particle_grid_cell_index_;             // 每个粒子所在的cell的索引
	uint*      next_particle_index_in_the_same_cell_; // 同一个cell中下个粒子的索引
	uint*      neighbor_index_;
	uint*      neighbor_particle_numbers_;

	// 加速数据结构---网格相关变量
	uint*     grid_head_cell_particle_index_array_;  // grid中每个cell中第一个粒子的索引
	uint*     grid_particles_number_;                // grid中每个cell中粒子的数量
	int       grid_total_;                           // grid的数量
	Vector3DI grid_res_;                             // grid的分辨率，即x，y，z方向上各有多少个cell
	Vector3DF grid_min_;
	Vector3DF grid_max_;
	//Vector3DF grid_size_;
	Vector3DF grid_delta_;
	uint*     cluster_cell_;                         // cell集群
	int	      grid_adj_cnt_;
	int       grid_neighbor_cell_index_offset_[max_num_adj_grid_cells_cpu];

	// 加速数据结构---邻居表相关变量
	int    neighbor_particles_num_;
	int    neighbor_particles_max_num_;
	int*   neighbor_table_;
	float* neighbor_dist_;

	// 边界碰撞处理相关变量
	bool  addBoundaryForce;
	float maxBoundaryForce;
	float boundaryForceFactor;
	float forceDistance;


	// 私有函数
	// SPH函数
	void InsertParticlesCPU(const uint& num_particle);
	void ComputePressureGrid();
	void ComputeForceGrid();
	void AdvanceStepSimpleCollision(float time_step);

	int       getGridCell(const Vector3DF& pos, Vector3DI& gc);
	int       getGridCell(int p, Vector3DI& gc);
	Vector3DI getCell(int c);

	// 光滑核函数
	float KernelM4(float dist, float sr);
	float KernelM4Lut(float dist, float sr);
	float KernelPressureGrad(float dist, float sr);
	float KernelPressureGradLut(float dist, float sr);

	void ComputeGasConstAndTimeStep(float densityVariation);
	void ClearNeighborTable();

	// 边界碰撞函数
	Vector3DF BoxBoundaryForce(const uint i);
	void      CollisionHandlingSimScale(Vector3DF& pos, Vector3DF& vel);

	// 绘制函数
	void DrawParticleInfo(int p);
	void DrawGrid();
	void DrawDomain(Vector3DF& domain_min, Vector3DF& domain_max);
	void DrawText();

	void setExampleParams(bool bStart);
	void setKernels();
	void setSpacing();

	void ParseXML(std::string name, int id, bool bStart);
};


inline bool ParticleSystem::getToggle(int p) {
	return this->toggle_[p];
}

inline int ParticleSystem::getNumPoints() {
	return this->num_points_;
}

inline Vector3DF ParticleSystem::getGridRes() {
	return this->grid_res_;
}

inline int ParticleSystem::getSelected() {
	return this->selected_;
}

inline int ParticleSystem::getGridTotal() {
	return this->grid_total_;
}

inline int ParticleSystem::getGridAdjCnt() {
	return this->grid_adj_cnt_;
}

inline float ParticleSystem::getParam(int p) {
	return (float)this->param_[p];
}

inline Vector3DF ParticleSystem::getVec(int p) {
	return this->vec_[p];
}

inline float ParticleSystem::getDT() {
	return this->time_step_;
}

inline void ParticleSystem::setToggle(int p) {
	this->toggle_[p] = !this->toggle_[p];
}

inline void ParticleSystem::DrawParticleInfo() {
	DrawParticleInfo(this->selected_);
}

#endif