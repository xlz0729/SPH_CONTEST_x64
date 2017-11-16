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
#include <vector>
#include <fstream>

#include <time.h>

// 公共开源库头文件
#include "common\vector.h"
#include "common\xml_settings.h"
#include "common\gl_helper.h"
#include "common\camera3d.h"


// 宏定义
#define MAX_PARAM				60         // 参数个数
#define GRID_UCHAR				0xFF
#define UNDEFINE				0xFFFFFFFF // 设定0xFFFFFFFF为未定义值

#define RUN_CPU_SPH				0
#define RUN_CPU_PBF				1

// 标量参数
#define PRUN_MODE				0	
#define PMAXNUM					1	// 最大粒子数量
#define PEXAMPLE				2	
#define PSIMSIZE				3	  
#define PSIMSCALE				4	// 缩放比例
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
#define PDRAWMODE				23	// 绘制方式，已弃用
#define PDRAWSIZE				24	
#define PDRAWTEXT				26	// 绘制内容
#define PCLR_MODE				27
#define PPOINT_GRAV_AMT			28	
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
#define PTIME_PBF_STEP			49
#define	PDENSITYERRORFACTOR		50
#define PMINLOOPPCISPH			51  // PCI最小循环次数
#define PMAXLOOPPCISPH			52  // PCI最大循环次数
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
const int grid_adj_cnt = 27;


// 粒子数据结构
struct Particle
{
	Vector3DF pos;                                  // 位置
	Vector3DF deta_pos;                             // PBF中的位置修正量
	Vector3DF vel;                                  // 瞬时速度
	Vector3DF vel_eval;                             // 平均速度
	Vector3DF force;                                // 受力
	float     pressure;	                            // 压力(无方向)
	float     density;                              // 密度
	float     mass;                                 // 质量
	int	      particle_grid_cell_index;             // 所在网格的索引
	int       next_particle_index_in_the_same_cell;	// 同一个网格中下一个粒子的索引
	DWORD	  clr;									// 用于渲染
};

// 网格数据结构
struct Cell
{
	uint first_particle_index;   // 每个cell中第一个粒子的索引
	uint particle_number;        // 每个cell中粒子的数量
	Cell() : first_particle_index(UNDEFINE), particle_number(0){}
	Cell(uint index, uint number) : first_particle_index(index), particle_number(number){}
};


// 函数命名规范:
// set和get函数开头第一个字母小写，函数名中每个单词开头字母大写
//     其余函数开头第一个字母大写，函数名中每个单词开头字母大写
//
// 变量命名规则:
// 全局类变量每个单词用_分割，所有单词小写
// 函数内局部变量每个单词开头字母大写
class ParticleSystem
{

public:
	ParticleSystem();
	~ParticleSystem();

	// set函数
	void         setUp(bool bStart);
	void         setRender();
	inline void  setToggle(int p);
	inline void  setParam(int p, int v);
	inline void  setParam(int p, float v);
	inline float IncParam(int p, float v, float mn, float mx);
	inline void  setVec(int p, Vector3DF v);

	int SelectParticle(int x, int y, int wx, int wy, Camera3D& cam);

	// get函数
	std::string      getModeStr();
	inline bool      getToggle(int p);
	inline int       getNumPoints();
	inline Vector3DF getGridRes();
	inline int       getGridTotal();
	inline int       getGridAdjCnt();
	inline float     getParam(int p);
	inline Vector3DF getVec(int p);
	inline float     getDT();

	// 运行函数
	void Run();
	void RunCPUSPH();
	void RunCPUPBF();

	// 绘制函数
	void        Draw(Camera3D& cam, float rad);

private:
	static int const lut_size = 100000;

	float lut_kernel_m4_[lut_size];
	float lut_kernel_pressure_grad_[lut_size];

	std::string	scene_name;
	int	        texture_[1];
	int	        sphere_points;

	// 模拟参数
	float     param_[MAX_PARAM];
	Vector3DF vec_[MAX_PARAM];
	bool      toggle_[MAX_PARAM];

	// XML 设置文件
	XmlSettings	xml;

	// 绘制函数相关变量
	int	   vbo_[3];
	Image  image;

	// SPH光滑核函数系数
	float poly6_kern;
	float lap_kern;
	float spiky_kern;

	// 时间步长
	float time_;
	float time_step;
	float time_step_sph;
	int	  frame;

	// 粒子相关属性
	std::vector<Particle> points;
	int	                  points_number;

	// 加速数据结构---网格相关变量
	std::vector<Cell>	grid;
	int					grid_total;                           // grid中cell的数量
	Vector3DI			grid_res;                             // x，y，z方向上各有多少个cell
	Vector3DF			grid_min;
	Vector3DF			grid_max;
	Vector3DF			grid_size;
	Vector3DF			grid_delta;
	int					grid_neighbor_cell_index_offset_[grid_adj_cnt];
	int					grid_nonempty_cell_number;              // grid中有粒子的cell的数量
	int					grid_particle_number;					// grid中粒子的数量

	// 边界碰撞相关变量
	bool  addBoundaryForce;
	float maxBoundaryForce;
	float boundaryForceFactor;
	float forceDistance;


	// 私有函数
	// SPH函数
	void InsertParticlesCPU();
	void ComputePressureGrid();
	void ComputeForceGrid();
	void AdvanceStepSimpleCollision();
	void ComputeDensity() {}
	void PositionBasedFluid();

	// 加速数据结构---网格相关函数
	int       getGridCell(const Vector3DF& pos, Vector3DI& gc);
	int       getGridCell(int p, Vector3DI& gc);
	Vector3DI getCell(int c);

	// 光滑核函数
	float KernelM4(float dist, float sr);
	float KernelM4Lut(float dist, float sr);
	float KernelPressureGrad(float dist, float sr);
	float KernelPressureGradLut(float dist, float sr);

	void ComputeGasConstAndTimeStep();

	// 边界碰撞函数
	Vector3DF BoxBoundaryForce(const uint i);
	void      CollisionHandlingSimScale(Vector3DF& pos, Vector3DF& vel);

	// 绘制函数
	void DrawGrid();
	void DrawDomain(Vector3DF& domain_min, Vector3DF& domain_max);
	void DrawText();

	void setDefaultParams();
	void setExampleParams(bool bStart);
	void setKernels();
	void setSpacing();
	void setInitParticleVolumeFromFile(const Vector3DF& minVec, const Vector3DF& maxVec);
	void setInitParticleVolume(const Vector3DI& numParticlesXYZ, const Vector3DF& lengthXYZ, const float jitter);
	void setGridAllocate(const float border);

	void AllocateParticlesMemory(int cnt);

	// 读取XML文件
	void ParseXML(std::string name, int id, bool bStart);

	void Record(int param, std::string name, mint::Time& start);

	// 从文件中读取粒子模型
	int ReadInFluidParticles(const char* filename, Vector3DF& minVec, Vector3DF& maxVec);

	inline float frand();
};


inline bool ParticleSystem::getToggle(int p) {
	return this->toggle_[p];
}

inline int ParticleSystem::getNumPoints() {
	return this->points_number;
}

inline Vector3DF ParticleSystem::getGridRes() {
	return this->grid_res;
}

inline int ParticleSystem::getGridTotal() {
	return this->grid_total;
}

inline int ParticleSystem::getGridAdjCnt() {
	return grid_adj_cnt;
}

inline float ParticleSystem::getParam(int p) {
	return (float)this->param_[p];
}

inline Vector3DF ParticleSystem::getVec(int p) {
	return this->vec_[p];
}

inline float ParticleSystem::getDT() {
	return this->time_step;
}

inline void ParticleSystem::setToggle(int p) {
	this->toggle_[p] = !this->toggle_[p];
}

inline float ParticleSystem::frand() {
	return rand() / (float)RAND_MAX;
}

inline float ParticleSystem::IncParam(int p, float v, float mn, float mx) {
	param_[p] += v;
	if (param_[p] < mn) param_[p] = mn;
	if (param_[p] > mx) param_[p] = mn;
	return param_[p];
}

inline void ParticleSystem::setParam(int p, float v) {
	param_[p] = v;
}

inline void ParticleSystem::setParam(int p, int v) {
	param_[p] = (float)v;
}

inline void ParticleSystem::setVec(int p, Vector3DF v) {
	vec_[p] = v;
}

#endif