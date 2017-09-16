#include "fluid_system.h"
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

#include "fluid_system.h"

ParticleSystem::ParticleSystem() {
	frame_ = 0;
	time_step_ = 0.003;
	time_ = 0;

	// 程序运行方法---6种方法可选
	param_[PRUN_MODE] = RUN_CPU_SPH;

	// 是否加载3D模型数据
	toggle_[PUSELOADEDSCENE] = false;

	num_points_ = 0;
	max_points_ = 0;
	good_points_ = 0;

	param_[PEXAMPLE] = 0;

	pos_ = 0x0;
	predictedPosition_ = 0x0;
	vel_ = 0x0;
	vel_eval_ = 0x0;
	correction_pressure_force_ = 0x0;
	force_ = 0x0;
	sumGradW_ = 0x0;
	labetalv = 0x0;
	sumGradWDot_ = 0x0;
	pressure_ = 0x0;
	correction_pressure_ = 0x0;
	density_ = 0x0;
	predicted_density_ = 0x0;
	densityError_ = 0x0;
	max_predicted_density_array_ = 0x0;
	labata_ = 0x0;
	beiyong1 = 0x0;
	beiyong2 = 0x0;
	particle_grid_cell_index_ = 0x0;
	next_particle_index_in_the_same_cell_ = 0x0;
	index_ = 0x0;
	clr_ = 0x0;

	cluster_cell_ = 0x0;
	age_ = 0x0;
	neighbor_index_ = 0x0;
	neighbor_particle_numbers_ = 0x0;

	grid_head_cell_particle_index_array_ = 0x0;
	grid_particles_number_ = 0x0;
	grid_total_ = 0;
	grid_search_ = 0;
	grid_adj_cnt_ = 0;

	neighbor_particles_num_ = 0;
	neighbor_particles_max_num_ = 0;
	neighbor_table_ = 0x0;
	neighbor_dist_ = 0x0;

	pack_fluid_particle_buf_ = 0x0;
	pack_grid_buf_ = 0x0;

	selected_ = -1;

	if (!xml.Load("scene.xml")) {
		printf("ERROR: Problem loading scene.xml. Check formatting.\n");
		system("pause");
		exit(-1);
	}
}

ParticleSystem::~ParticleSystem() {}

// 程序运行入口
void ParticleSystem::setup(bool bStart) {
	frame_ = 0;
	time_ = 0;

	setRender();

	setExampleParams(bStart);

	param_[PGRIDSIZEREALSCALE] = param_[PSMOOTHRADIUS] / param_[PGRID_DENSITY];

	setKernels();

	setSpacing();

	ComputeGasConstAndTimeStep(param_[PMAXDENSITYERRORALLOWED]);

	/*if (param_[PRUN_MODE] == RUN_CPU_PCISPH || param_[PRUN_MODE] == RUN_CUDA_INDEX_PCISPH || param_[PRUN_MODE] == RUN_CUDA_FULL_PCISPH)
	{
		time_step_ = time_step_pcisph_;
		ClearNeighborTable();
		AllocateTemporalParticlesMemory(param_[PMAXNUM]);
		Vector3DF init_sample_particle_volume_min = Vector3DF(-20.0, 0.0, -20);
		Vector3DF init_sample_particle_volume_max = Vector3DF(20.0, 40.0, 20);
		SetupInitParticleVolume(init_sample_particle_volume_min, init_sample_particle_volume_max, param_[PSPACINGGRAPHICSWORLD], 0.1);
		Vector3DF init_grid_volume_min = init_sample_particle_volume_min - (param_[PSMOOTHRADIUS] / param_[PSIMSCALE]);
		Vector3DF init_grid_volume_max = init_sample_particle_volume_max + (param_[PSMOOTHRADIUS] / param_[PSIMSCALE]);
		SetupSampleGridAllocatePCISPH(init_grid_volume_min, init_grid_volume_max, param_[PSIMSCALE], param_[PGRIDSIZEREALSCALE], 1.0);
		ComputeDensityErrorFactor(num_points_);
		DeallocateTemporalParticleMemory();
	}*/

	ClearNeighborTable();

	// 从文件中加载粒子位置信息,并设置其速度，颜色等属性值
	if (toggle_[PUSELOADEDSCENE] == true)
	{
		// 从文件中读取bunny模型数据
		particles.clear();
		const char* file_name = "models/bunny.txt";
		num_points_ = readInFluidParticleNum(file_name);

		Vector3DF minCorner = vec_[PINITPARTICLEMIN];
		Vector3DF maxCorner = vec_[PINITPARTICLEMAX];
		const float lengthX = maxCorner.x - minCorner.x;
		const float lengthY = maxCorner.y - minCorner.y;
		const float lengthZ = maxCorner.z - minCorner.z;

		const int numParticlesX = ceil(lengthX / param_[PSPACINGGRAPHICSWORLD]);
		const int numParticlesY = ceil(lengthY / param_[PSPACINGGRAPHICSWORLD]);
		const int numParticlesZ = ceil(lengthZ / param_[PSPACINGGRAPHICSWORLD]);
		const int numParticles = numParticlesX * numParticlesY * numParticlesZ;

		AllocateParticlesMemory(num_points_ + numParticles);

		Vector3DF minVec;
		Vector3DF maxVec;
		readInFluidParticles(file_name, num_points_, minVec, maxVec);
		SetupInitParticleVolumeFromFile(minVec, maxVec);
		SetupAdditonalParticleVolume(vec_[PINITPARTICLEMIN], vec_[PINITPARTICLEMAX], param_[PSPACINGGRAPHICSWORLD], 0.1, numParticles);

		num_points_ = num_points_ + numParticles;
	}
	else
	{
		num_points_ = 0;

		AllocateParticlesMemory(param_[PMAXNUM]);

		SetupInitParticleVolume(vec_[PINITPARTICLEMIN], vec_[PINITPARTICLEMAX], param_[PSPACINGGRAPHICSWORLD], 0.1);
	}

	AllocatePackBuf();

	SetupGridAllocate(vec_[PGRIDVOLUMEMIN], vec_[PGRIDVOLUMEMAX], param_[PSIMSCALE], param_[PGRIDSIZEREALSCALE], 1.0);

/*
#ifdef USE_CUDA

	ParticleClearCUDA();

	int3 grid_res = make_int3(grid_res_.x, grid_res_.y, grid_res_.z);
	float3 grid_size = make_float3(grid_size_.x, grid_size_.y, grid_size_.z);
	float3 grid_delta = make_float3(grid_delta_.x, grid_delta_.y, grid_delta_.z);
	float3 grid_min = make_float3(grid_min_.x, grid_min_.y, grid_min_.z);
	float3 grid_max = make_float3(grid_max_.x, grid_max_.y, grid_max_.z);
	ParticleSetupCUDA(num_points(), grid_search_, grid_res, grid_size, grid_delta, grid_min, grid_max, grid_total_, (int)vec_[PEMIT_RATE].x, param_[PGRIDSIZEREALSCALE], param_[PKERNELSELF]);

	Vector3DF grav = vec_[PPLANE_GRAV_DIR];
	float3 boundMin = make_float3(vec_[PBOUNDARYMIN].x, vec_[PBOUNDARYMIN].y, vec_[PBOUNDARYMIN].z);
	float3 boundMax = make_float3(vec_[PBOUNDARYMAX].x, vec_[PBOUNDARYMAX].y, vec_[PBOUNDARYMAX].z);
	FluidParamCUDA(param_[PSIMSCALE], param_[PSMOOTHRADIUS], param_[PCOLLISIONRADIUS], param_[PMASS], param_[PRESTDENSITY], boundMin, boundMax,
		param_[PBOUNDARYSTIFF], param_[PGASCONSTANT], param_[PVISC], param_[PBOUNDARYDAMP], param_[PFORCE_MIN], param_[PFORCE_MAX],
		param_[PFORCE_FREQ], param_[PGROUND_SLOPE], grav.x, grav.y, grav.z, param_[PACCEL_LIMIT], param_[PVEL_LIMIT], param_[PDENSITYERRORFACTOR]);

	TransferToCUDA();		// 数据拷贝CPU -> GPU

#endif 
*/
}

void ParticleSystem::setRender() {

	param_[PMAXNUM] = 1048576;             // 最大粒子数量
	param_[PSIMSCALE] = 0.005;             // 缩放比例1:2
	param_[PGRID_DENSITY] = 1.0;           // grid的密度
	param_[PVISC] = 1.002;                 // 体积
	param_[PRESTDENSITY] = 1000.0;         // 密度
	param_[PMASS] = 0.001953125;           // 质量
	param_[PCOLLISIONRADIUS] = 0.00775438; // 碰撞半径
	param_[PSPACINGREALWORLD] = 0.0125;    // 真实世界的间隔
	param_[PSMOOTHRADIUS] = 0.025025;      // 光滑半径
	param_[PGRIDSIZEREALSCALE] = param_[PSMOOTHRADIUS] / param_[PGRID_DENSITY]; // grid的真实缩放比例
	param_[PGASCONSTANT] = 1000.0;         // 气体含量
	param_[PBOUNDARYSTIFF] = 2.5;          // 边界弹性
	param_[PBOUNDARYDAMP] = 1.0;           // 边界阻力
	param_[PACCEL_LIMIT] = 150.0;          // 加速度上限
	param_[PVEL_LIMIT] = 3.0;              // 速度上限
	param_[PSPACINGGRAPHICSWORLD] = 0.0;   // 图形世界的间隔
	param_[PGROUND_SLOPE] = 0.0;           // 地面倾斜度
	param_[PFORCE_MIN] = 0.0;              // 最小力
	param_[PFORCE_MAX] = 0.0;              // 最大力
	param_[PDRAWMODE] = 2;                 // 绘制方式
	param_[PDRAWTEXT] = 0;                 // 绘制内容
	param_[PPOINT_GRAV_AMT] = 0.0;         // 
	param_[PSTAT_NEIGHCNTMAX] = 0;         // 领域内最大粒子数
	param_[PSTAT_SEARCHCNTMAX] = 0;        // 搜索范围内最大粒子数
	param_[PFORCE_FREQ] = 8.0;             // 受力频率
	param_[PMINLOOPPCISPH] = 3;            // PCI最小循环次数
	param_[PMAXLOOPPCISPH] = MAX_PCISPH_LOOPS; // PCI最大循环次数
	param_[PMAXDENSITYERRORALLOWED] = 5.0; // 允许的最大密度误差

	vec_[PEMIT_POS].Set(0, 0, 0);
	vec_[PEMIT_ANG].Set(0, 90, 1.0);
	vec_[PEMIT_RATE].Set(0, 0, 0);
	vec_[PPOINT_GRAV_POS].Set(0, 50, 0);   // 点受力的位置
	vec_[PPLANE_GRAV_DIR].Set(0, -9.8, 0); // 重力方向

	toggle_[PPAUSE] = false;
	toggle_[PWRAP_X] = false;
	toggle_[PWALL_BARRIER] = false;
	toggle_[PLEVY_BARRIER] = false;
	toggle_[PDRAIN_BARRIER] = false;
	toggle_[PPROFILE] = false;
	toggle_[PCAPTURE] = false;
	toggle_[PPRINTDEBUGGINGINFO] = false;
	toggle_[PDRAWDOMAIN] = false;
	toggle_[PDRAWGRIDBOUND] = false;
	toggle_[PDRAWGRIDCELLS] = false;
}

std::string ParticleSystem::getModeStr() {

	std::string buf = "SIMULATE CPU SPH";

	switch ((int)param_[PRUN_MODE]) {
	case RUN_CPU_SPH:
		buf = "SIMULATE CPU SPH";
		break;
	case RUN_CPU_PCISPH:
		buf = "SIMULATE CPU PCISPH";
		break;
	};
	return buf;
};

void ParticleSystem::Run() {

	// 计时器清零
	param_[PTIME_INSERT] = 0.0;
	param_[PTIME_SORT] = 0.0;
	param_[PTIME_COUNT] = 0.0;
	param_[PTIME_PRESS] = 0.0;
	param_[PTIME_FORCE] = 0.0;
	param_[PTIME_ADVANCE] = 0.0;
	param_[PTIME_OTHER_FORCE] = 0.0;
	param_[PTIME_PCI_STEP] = 0.0;

	// 运行程序
	switch ((int)param_[PRUN_MODE]) {
	case RUN_CPU_SPH:
		RunCPUSPH();
		break;
	case RUN_CPU_PCISPH:
		RunCPUPCISPH();
		break;
	}

	time_ += time_step_;
	frame_++;
}

void ParticleSystem::RunCPUSPH() {

	// 插入粒子
	InsertParticlesCPU(this->num_points_);

	// 计算压力
	ComputePressureGrid();

	// 计算外力
	ComputeForceGrid();

	// 移动粒子
	AdvanceStepSimpleCollision(time_step_);
}

void ParticleSystem::InsertParticlesCPU(const uint& num_particle) {

	memset(next_particle_index_in_the_same_cell_, GRID_UNDEF, num_particle * sizeof(uint));
	memset(particle_grid_cell_index_,             GRID_UNDEF, num_particle * sizeof(uint));
	memset(cluster_cell_,                         GRID_UNDEF, num_particle * sizeof(uint));

	memset(grid_head_cell_particle_index_array_, GRID_UNDEF, grid_total_ * sizeof(uint));
	memset(grid_particles_number_,               0,          grid_total_ * sizeof(uint));

	const int xns = grid_res_.x;
	const int yns = grid_res_.y;
	const int zns = grid_res_.z;

	param_[PSTAT_OCCUPANCY] = 0.0; // 有粒子的grid的数量
	param_[PSTAT_GRIDCOUNT] = 0.0; // grid中粒子的总数

	for (int idx = 0; idx < num_particle; ++idx) {
		Vector3DI gridCell;
		const int gridCellIndex = getGridCell(pos_[idx], gridCell);

		if (gridCell.x >= 0 && gridCell.x < xns && gridCell.y >= 0 && gridCell.y < yns && gridCell.z >= 0 && gridCell.z < zns) {
			particle_grid_cell_index_[idx] = gridCellIndex;
			next_particle_index_in_the_same_cell_[idx] = grid_head_cell_particle_index_array_[gridCellIndex];
			if (next_particle_index_in_the_same_cell_[idx] == GRID_UNDEF) {
				param_[PSTAT_OCCUPANCY] += 1.0;
			}
			grid_head_cell_particle_index_array_[gridCellIndex] = idx;
			grid_particles_number_[gridCellIndex] += 1;
			param_[PSTAT_GRIDCOUNT] += 1.0;
		}
		else {
			Vector3DF vel = *(vel_ + idx);
			Vector3DF ve = *(vel_eval_ + idx);
			float pr = *(pressure_ + idx);
			float dn = *(density_ + idx);
			printf("WARNING: ParticleSystem::InsertParticlesCPU(): Out of Bounds: %d, Position<%f %f %f>, Velocity<%f %f %f>, Pressure:%f, Density:%f\n", idx, pos_[idx].x, pos_[idx].y, pos_[idx].z, vel.x, vel.y, vel.z, pr, dn);
			pos_[idx].x = -1; pos_[idx].y = -1; pos_[idx].z = -1;
		}
	}
}

void ParticleSystem::ComputePressureGrid() {

	const float	sim_scale_square = param_[PSIMSCALE] * param_[PSIMSCALE];
	const float smooth_radius = param_[PSMOOTHRADIUS];
	const float smooth_radius_square = smooth_radius * smooth_radius;
	const float mass = param_[PMASS];

	const float own_density_contribution = param_[PKERNELSELF] * mass;

	float minDens = 10e10;
	float maxDens = 0.0;
	for (int i = 0; i < this->num_points_; ++i) {

		density_[i] = own_density_contribution;

		int neighbor_nums = 0;
		int search_nums = 0;
		float sum = 0.0;

		const uint i_cell_index = particle_grid_cell_index_[i];
		if (i_cell_index != GRID_UNDEF) {
			for (int cell = 0; cell < max_num_adj_grid_cells_cpu; cell++) {
				const int neighbor_cell_index = i_cell_index + grid_neighbor_cell_index_offset_[cell];
				if (neighbor_cell_index == GRID_UNDEF || neighbor_cell_index < 0 || neighbor_cell_index > grid_total_ - 1) {
					continue;
				}

				int j = grid_head_cell_particle_index_array_[neighbor_cell_index];
				while (j != GRID_UNDEF) {
					if (i == j) {
						j = next_particle_index_in_the_same_cell_[j];
						continue;
					}
					Vector3DF dst_graphics_scale = pos_[j];
					dst_graphics_scale -= pos_[i];
					const float dist_square_sim_scale = sim_scale_square*(dst_graphics_scale.x*dst_graphics_scale.x + dst_graphics_scale.y*dst_graphics_scale.y + dst_graphics_scale.z*dst_graphics_scale.z);
					if (dist_square_sim_scale <= smooth_radius_square) {
						const float dist = sqrt(dist_square_sim_scale);
						float kernelValue = KernelM4Lut(dist, smooth_radius);
						density_[i] += kernelValue * mass;

						neighbor_nums++;
					}
					search_nums++;
					j = next_particle_index_in_the_same_cell_[j];
				}
			}
		}

		if (density_[i] < minDens)
			minDens = density_[i];
		if (density_[i] > maxDens)
			maxDens = density_[i];

		pressure_[i] = max(0.0f, (density_[i] - param_[PRESTDENSITY]) * param_[PGASCONSTANT]);

		param_[PSTAT_NEIGHCNT] = float(neighbor_nums);
		param_[PSTAT_SEARCHCNT] = float(search_nums);
		if (param_[PSTAT_NEIGHCNT] > param_[PSTAT_NEIGHCNTMAX])
			param_[PSTAT_NEIGHCNTMAX] = param_[PSTAT_NEIGHCNT];
		if (param_[PSTAT_SEARCHCNT] > param_[PSTAT_SEARCHCNTMAX])
			param_[PSTAT_SEARCHCNTMAX] = param_[PSTAT_SEARCHCNT];
	}
}

void ParticleSystem::ComputeForceGrid() {

	const float mass = param_[PMASS];
	const float sim_scale = param_[PSIMSCALE];
	const float sim_scale_square = sim_scale * sim_scale;
	const float smooth_radius = param_[PSMOOTHRADIUS];
	const float smooth_radius_square = smooth_radius * smooth_radius;
	const float visc = param_[PVISC];
	Vector3DF   vec_gravity = vec_[PPLANE_GRAV_DIR];
	const float vterm = lap_kern_ * visc;

	for (int i = 0; i < this->num_points_; i++) {
		force_[i].Set(0, 0, 0);
		Vector3DF force(0, 0, 0);
		const uint i_cell_index = particle_grid_cell_index_[i];
		Vector3DF ipos = pos_[i];
		Vector3DF iveleval = vel_eval_[i];
		float	  ipress = pressure_[i];
		float	  idensity = density_[i];

		if (i_cell_index != GRID_UNDEF) {
			for (int cell = 0; cell < max_num_adj_grid_cells_cpu; cell++) {
				const int neighbor_cell_index = i_cell_index + grid_neighbor_cell_index_offset_[cell];
				if (neighbor_cell_index == GRID_UNDEF || (neighbor_cell_index < 0 || neighbor_cell_index > grid_total_ - 1)) {
					continue;
				}

				int j = grid_head_cell_particle_index_array_[neighbor_cell_index];
				while (j != GRID_UNDEF) {
					if (i == j) {
						j = next_particle_index_in_the_same_cell_[j];
						continue;
					}

					Vector3DF vector_i_minus_j = ipos - pos_[j];
					const float dx = vector_i_minus_j.x;
					const float dy = vector_i_minus_j.y;
					const float dz = vector_i_minus_j.z;

					const float dist_square_sim_scale = sim_scale_square*(dx*dx + dy*dy + dz*dz);
					if (dist_square_sim_scale <= smooth_radius_square && dist_square_sim_scale > 0) {
						const float jdist = sqrt(dist_square_sim_scale);
						const float jpress = pressure_[j];
						const float h_minus_r = smooth_radius - jdist;
						const float pterm = -0.5f * h_minus_r * spiky_kern_ * (ipress + jpress) / jdist;
						const float dterm = h_minus_r / (idensity * density_[j]);

						Vector3DF vel_j_minus_i = vel_eval_[j];
						vel_j_minus_i -= iveleval;

						force += vector_i_minus_j * sim_scale * pterm * dterm;

						force += vel_j_minus_i * vterm * dterm;
					}
					j = next_particle_index_in_the_same_cell_[j];
				}
			}
		}

		force *= mass * mass;
		force += vec_gravity * mass;

		if (addBoundaryForce) {
			force += BoxBoundaryForce(i);
		}

		force_[i] = force;

	}
}

void ParticleSystem::AdvanceStepSimpleCollision(float time_step) {

	const float acc_limit = param_[PACCEL_LIMIT];
	const float acc_limit_square = acc_limit*acc_limit;
	const float speed_limit = param_[PVEL_LIMIT];
	const float speed_limit_square = speed_limit*speed_limit;

	const float sim_scale = param_[PSIMSCALE];

	Vector3DF norm;
	Vector4DF clr;
	float adj;
	float speed;
	float diff;

	for (int i = 0; i < this->num_points_; i++) {
		if (particle_grid_cell_index_[i] == GRID_UNDEF)
			continue;

		Vector3DF acceleration = force_[i] + correction_pressure_force_[i];
		acceleration /= param_[PMASS];

		vel_eval_[i] += acceleration * time_step;
		pos_[i] *= sim_scale;
		pos_[i] += vel_eval_[i] * time_step;

		CollisionHandlingSimScale(pos_[i], vel_eval_[i]);

		pos_[i] /= sim_scale;
	}
}

// 得到当前pos所在的gridCell的编号
int ParticleSystem::getGridCell(const Vector3DF& pos, Vector3DI& gc) {

	float px = pos.x - vec_[PGRIDVOLUMEMIN].x;
	float py = pos.y - vec_[PGRIDVOLUMEMIN].y;
	float pz = pos.z - vec_[PGRIDVOLUMEMIN].z;

	if (px < 0.0)
		px = 0.0;
	if (py < 0.0)
		py = 0.0;
	if (pz < 0.0)
		pz = 0.0;

	const float cellSize = param_[PGRIDSIZEREALSCALE] / param_[PSIMSCALE];
	gc.x = (int)(px / cellSize);
	gc.y = (int)(py / cellSize);
	gc.z = (int)(pz / cellSize);

	if (gc.x > grid_res_.x - 1)
		gc.x = grid_res_.x - 1;
	if (gc.y > grid_res_.y - 1)
		gc.y = grid_res_.y - 1;
	if (gc.z > grid_res_.z - 1)
		gc.z = grid_res_.z - 1;

	// 将3维的cell索引转化为1维的整数索引
	return (int)(gc.y * grid_res_.x  * grid_res_.z + gc.z * grid_res_.x + gc.x); // x*y*z + x*y + x
}

int ParticleSystem::getGridCell(int p, Vector3DI& gc){
	return getGridCell(*(pos_ + p), gc);
}

Vector3DI ParticleSystem::getCell(int c) {

	Vector3DI gc;
	int xz = grid_res_.x*grid_res_.z;
	gc.y = c / xz;
	c -= gc.y*xz;
	gc.z = c / grid_res_.x;
	c -= gc.z*grid_res_.x;
	gc.x = c;
	return gc;
}

float ParticleSystem::KernelM4(float dist, float sr) {

	float s = dist / sr;
	float result;
	float factor = 2.546479089470325472f / (sr * sr * sr);
	if (dist < 0.0f || dist >= sr)
		return 0.0f;
	else
	{
		if (s < 0.5f)
		{
			result = 1.0f - 6.0 * s * s + 6.0f * s * s * s;
		}
		else
		{
			float tmp = 1.0f - s;
			result = 2.0 * tmp * tmp * tmp;
		}
	}
	return factor * result;
}

float ParticleSystem::KernelM4Lut(float dist, float sr) {

	int index = dist / sr * lutSize;

	if (index >= lutSize)
		return 0.0f;
	else
		return lutKernelM4[index];
}

float ParticleSystem::KernelPressureGrad(float dist, float sr) {

	if (dist == 0)
		return 0.0f;
	if (dist > sr)
		return 0.0f;

	float kernelPressureConst = -45.f / ((float(MY_PI)*sr*sr*sr*sr*sr*sr));
	return kernelPressureConst / dist * (sr - dist)*(sr - dist);
}

float ParticleSystem::KernelPressureGradLut(float dist, float sr) {

	int index = dist / sr * lutSize;
	if (index >= lutSize) return 0.0f;
	else return lutKernelPressureGrad[index];
}

void ParticleSystem::ComputeGasConstAndTimeStep(float densityVariation) {

	float maxParticleSpeed = 4.0f;
	float courantFactor = 0.4f;

	if (densityVariation >= 1.0) {
		time_step_pcisph_ = 0.001;

		param_[PGASCONSTANT] = 70000;
		float speedOfSound = sqrt(param_[PGASCONSTANT]);
		float relevantSpeed = max(speedOfSound, maxParticleSpeed);
		time_step_wcsph_ = courantFactor * param_[PSMOOTHRADIUS] / relevantSpeed;

		param_[PGASCONSTANT] = 1000.0f;
		speedOfSound = sqrt(param_[PGASCONSTANT]);
		relevantSpeed = max(speedOfSound, maxParticleSpeed);
		time_step_sph_ = courantFactor * param_[PSMOOTHRADIUS] / relevantSpeed;
	}
	else {
		time_step_pcisph_ = 0.0005;

		param_[PGASCONSTANT] = 6000000;
		float speedOfSound = sqrt(param_[PGASCONSTANT]);
		float relevantSpeed = max(speedOfSound, maxParticleSpeed);
		time_step_wcsph_ = courantFactor * param_[PSMOOTHRADIUS] / relevantSpeed;

		param_[PGASCONSTANT] = 1000.0f;
		speedOfSound = sqrt(param_[PGASCONSTANT]);
		relevantSpeed = max(speedOfSound, maxParticleSpeed);
		time_step_sph_ = courantFactor * param_[PSMOOTHRADIUS] / relevantSpeed;
	}

	time_step_ = time_step_sph_;

}

void ParticleSystem::ClearNeighborTable() {

	if (neighbor_table_ != 0x0)
		free(neighbor_table_);
	if (neighbor_dist_ != 0x0)
		free(neighbor_dist_);
	neighbor_table_ = 0x0;
	neighbor_dist_ = 0x0;
	neighbor_particles_num_ = 0;
	neighbor_particles_max_num_ = 0;
}

// 边界碰撞
Vector3DF ParticleSystem::BoxBoundaryForce(const uint i) {

	static const float invForceDist = 1.0 / forceDistance;
	static const float forceStrength = maxBoundaryForce;
	static const Vector3DF bound_min = vec_[PBOUNDARYMIN];
	static const Vector3DF bound_max = vec_[PBOUNDARYMAX];
	float distToWall, factor;
	Vector3DF force(0.0, 0.0, 0.0);

	if (pos_[i].y < bound_min.y + forceDistance) {
		distToWall = bound_min.y + forceDistance - pos_[i].y;
		factor = distToWall * invForceDist;
		force += Vector3DF(0, 1, 0) * factor * 2.0f * forceStrength;
	}
	else if (pos_[i].y > bound_max.y - forceDistance) {
		distToWall = pos_[i].y - (bound_max.y - forceDistance);
		factor = distToWall * invForceDist;
		force += Vector3DF(0, -1, 0) * factor * forceStrength;
	}

	if (pos_[i].x < bound_min.x + forceDistance) {
		distToWall = bound_min.x + forceDistance - pos_[i].x;
		factor = distToWall * invForceDist;
		force += Vector3DF(1, 0, 0) * factor * forceStrength;
	}
	if (pos_[i].x > bound_max.x - forceDistance) {
		distToWall = pos_[i].x - (bound_max.x - forceDistance);
		factor = distToWall * invForceDist;
		force += Vector3DF(-1, 0, 0) * 1 * factor * forceStrength;
	}

	if (pos_[i].z < bound_min.z + forceDistance) {
		distToWall = bound_min.z + forceDistance - pos_[i].z;
		factor = distToWall * invForceDist;
		force += Vector3DF(0, 0, 1) * factor * forceStrength;
	}
	if (pos_[i].z > bound_max.z - forceDistance) {
		distToWall = pos_[i].z - (bound_max.z - forceDistance);
		factor = distToWall * invForceDist;
		force += Vector3DF(0, 0, -1) * factor * forceStrength;
	}

	return force;
}

void ParticleSystem::CollisionHandlingSimScale(Vector3DF& pos, Vector3DF& vel) {

	const float     sim_scale = param_[PSIMSCALE];
	const Vector3DF vec_bound_min = vec_[PBOUNDARYMIN] * sim_scale;
	const Vector3DF vec_bound_max = vec_[PBOUNDARYMAX] * sim_scale;

	float damping = 0.1;

	float reflect = 1.1;

	// 碰撞处理
	if (pos.x < vec_bound_min.x) {
		pos.x = vec_bound_min.x;
		Vector3DF axis(1, 0, 0);
		vel -= axis * (float)axis.Dot(vel) * reflect;
		vel.x *= damping;
	}

	if (pos.x > vec_bound_max.x) {
		pos.x = vec_bound_max.x;
		Vector3DF axis(-1, 0, 0);
		vel -= axis * (float)axis.Dot(vel) * reflect;
		vel.x *= damping;
	}

	if (pos.y < vec_bound_min.y) {
		pos.y = vec_bound_min.y;
		Vector3DF axis(0, 1, 0);
		vel -= axis * (float)axis.Dot(vel) * reflect;
		vel.y *= damping;
	}

	if (pos.y > vec_bound_max.y) {
		pos.y = vec_bound_max.y;
		Vector3DF axis(0, -1, 0);
		vel -= axis * (float)axis.Dot(vel) * reflect;
		vel.y *= damping;
	}

	if (pos.z < vec_bound_min.z) {
		pos.z = vec_bound_min.z;
		Vector3DF axis(0, 0, 1);
		vel -= axis * (float)axis.Dot(vel) * reflect;
		vel.z *= damping;
	}

	if (pos.z > vec_bound_max.z) {
		pos.z = vec_bound_max.z;
		Vector3DF axis(0, 0, -1);
		vel -= axis * (float)axis.Dot(vel) * reflect;
		vel.z *= damping;
	}
}

void ParticleSystem::DrawParticleInfo(int p) {

	char disp[256];

	glColor4f(1.0, 1.0, 1.0, 1.0);
	sprintf(disp, "Particle: %d", p);
	drawText2D(10, 20, disp);

	Vector3DI gc;
	int gs = getGridCell(p, gc);
	sprintf(disp, "Grid Cell:    <%d, %d, %d> id: %d", gc.x, gc.y, gc.z, gs);		drawText2D(10, 40, disp);

	int cc = *(cluster_cell_ + p);
	gc = getCell(cc);
	sprintf(disp, "Cluster Cell: <%d, %d, %d> id: %d", gc.x, gc.y, gc.z, cc);		drawText2D(10, 50, disp);

	sprintf(disp, "Neighbors:    ");

	int cnt = *(neighbor_particle_numbers_ + p);
	int ndx = *(neighbor_index_ + p);
	for (int n = 0; n < cnt; n++) {
		sprintf(disp, "%s%d, ", disp, neighbor_table_[ndx]);
		ndx++;
	}
	drawText2D(10, 70, disp);

	if (cc != -1) {
		sprintf(disp, "Cluster Group: ");		drawText2D(10, 90, disp);
		int cadj;
		int stotal = 0;
		for (int n = 0; n < max_num_adj_grid_cells_cpu; n++) {
			cadj = cc + grid_neighbor_cell_index_offset_[n];

			if (cadj == GRID_UNDEF || (cadj < 0 || cadj > grid_total_ - 1))
			{
				continue;
			}
			gc = getCell(cadj);
			sprintf(disp, "<%d, %d, %d> id: %d, cnt: %d ", gc.x, gc.y, gc.z, cc + grid_neighbor_cell_index_offset_[n], grid_particles_number_[cadj]);
			drawText2D(20, 100 + n * 10, disp);
			stotal += grid_particles_number_[cadj];
		}

		sprintf(disp, "Search Overhead: %f (%d of %d), %.2f%% occupancy", float(stotal) / cnt, cnt, stotal, float(cnt)*100.0 / stotal);
		drawText2D(10, 380, disp);
	}
}

// 绘制函数
void ParticleSystem::Draw(Camera3D& cam, float rad) {

	float* pdens;

	glDisable(GL_LIGHTING);

	if (toggle_[PDRAWGRIDCELLS]) {
		glColor4f(0.0, 0.0, 1.0, 0.1);
		DrawGrid();
	}

	if (toggle_[PDRAWDOMAIN]) {
		DrawDomain(vec_[PBOUNDARYMIN], vec_[PBOUNDARYMAX]);
	}

	if (toggle_[PDRAWGRIDBOUND]) {
		DrawDomain(vec_[PGRIDVOLUMEMIN], vec_[PGRIDVOLUMEMAX]);
	}

	if (param_[PDRAWTEXT] == 1.0) {
		DrawText();
	};

	switch ((int)param_[PDRAWMODE]) {
	case 0:
		{
			glPointSize(6);
			glEnable(GL_POINT_SIZE);
			glEnable(GL_BLEND);
			glBindBuffer(GL_ARRAY_BUFFER, vbo_[0]);
			glBufferData(GL_ARRAY_BUFFER, this->num_points_ * sizeof(Vector3DF), pos_, GL_DYNAMIC_DRAW);
			glVertexPointer(3, GL_FLOAT, 0, 0x0);
			glBindBuffer(GL_ARRAY_BUFFER, vbo_[1]);
			glBufferData(GL_ARRAY_BUFFER, this->num_points_ * sizeof(uint), clr_, GL_DYNAMIC_DRAW);
			glColorPointer(4, GL_UNSIGNED_BYTE, 0, 0x0);
			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_COLOR_ARRAY);
			glNormal3f(0, 0.001, 1);
			glColor3f(1, 1, 1);
			glDrawArrays(GL_POINTS, 0, this->num_points_);
			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_COLOR_ARRAY);
		}
		break;

	case 1:
		{
			glEnable(GL_LIGHTING);
			glEnable(GL_BLEND);
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GREATER, 0.5);
			glEnable(GL_COLOR_MATERIAL);
			// 根据glColor所设置的值来指定材料参数
			glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

			// 点精灵大小
			glPointSize(32);
			glEnable(GL_POINT_SIZE);
			glEnable(GL_POINT_SPRITE);
			float quadratic[] = { 0.0f, 0.3f, 0.00f };
			glEnable(GL_POINT_DISTANCE_ATTENUATION);
			glPointParameterfv(GL_POINT_DISTANCE_ATTENUATION, quadratic);
			float maxSize = 64.0f;
			glGetFloatv(GL_POINT_SIZE_MAX, &maxSize);
			glPointSize(maxSize);
			glPointParameterf(GL_POINT_SIZE_MAX, maxSize);
			glPointParameterf(GL_POINT_SIZE_MIN, 1.0f);

			// 纹理&混合模式
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, image_.getID());
			glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

			glBindBuffer(GL_ARRAY_BUFFER, vbo_[0]);
			glBufferData(GL_ARRAY_BUFFER, this->num_points_ * sizeof(Vector3DF), pos_, GL_DYNAMIC_DRAW);
			glVertexPointer(3, GL_FLOAT, 0, 0x0);
			glBindBuffer(GL_ARRAY_BUFFER, vbo_[1]);
			glBufferData(GL_ARRAY_BUFFER, this->num_points_ * sizeof(uint), clr_, GL_DYNAMIC_DRAW);
			glColorPointer(4, GL_UNSIGNED_BYTE, 0, 0x0);
			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_COLOR_ARRAY);

			glNormal3f(0, 1, 0.001);
			glColor3f(1, 1, 1);
			glDrawArrays(GL_POINTS, 0, this->num_points_);

			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_COLOR_ARRAY);
			glDisable(GL_POINT_SPRITE);
			glDisable(GL_ALPHA_TEST);
			glDisable(GL_TEXTURE_2D);
		}
		break;

	case 2:
		{
			glEnable(GL_LIGHTING);
			pdens = density_;

			for (int n = 0; n < this->num_points_; n++) {
				glPushMatrix();
				glTranslatef(pos_[n].x, pos_[n].y, pos_[n].z);
				glScalef(rad, rad, rad);
				glColor4f(RED(clr_[n]), GRN(clr_[n]), BLUE(clr_[n]), ALPH(clr_[n]));
				drawSphere();
				glPopMatrix();
			}
		}
		break;
	}
}

// 绘制grid
void ParticleSystem::DrawGrid() {

	Vector3DF gd(1, 1, 1);
	Vector3DF gc;
	gd /= grid_delta_;

	glBegin(GL_LINES);
	for (int z = 0; z <= grid_res_.z; z++) {
		for (int y = 0; y <= grid_res_.y; y++) {
			gc.Set(1, y, z);
			gc /= grid_delta_;
			gc += grid_min_;
			glVertex3f(grid_min_.x, gc.y, gc.z);
			glVertex3f(grid_max_.x, gc.y, gc.z);
		}
	}
	for (int z = 0; z <= grid_res_.z; z++) {
		for (int x = 0; x <= grid_res_.x; x++) {
			gc.Set(x, 1, z);
			gc /= grid_delta_;
			gc += grid_min_;
			glVertex3f(gc.x, grid_min_.y, gc.z);
			glVertex3f(gc.x, grid_max_.y, gc.z);
		}
	}
	for (int y = 0; y <= grid_res_.y; y++) {
		for (int x = 0; x <= grid_res_.x; x++) {
			gc.Set(x, y, 1);
			gc /= grid_delta_;
			gc += grid_min_;
			glVertex3f(gc.x, gc.y, grid_min_.z);
			glVertex3f(gc.x, gc.y, grid_max_.z);
		}
	}
	glEnd();
}

// 绘制边界
void ParticleSystem::DrawDomain(Vector3DF& domain_min, Vector3DF& domain_max) {
	/*****************************************************
	   /8------------- /7
	  / |             / |
	 /	|			 /  |
	/5--|-----------6	|
	|   |			|   |
	|	|			|	|
	|	|			|	|
	|	4-----------|---3
	|  /			|  /
	| /			    | /
	|1 -------------2

	*******************************************************/

	glColor3f(1.0, 0.0, 0.0);

	// ground (1234)
	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_min.y, domain_max.z);//1
	glVertex3f(domain_max.x, domain_min.y, domain_max.z);//2
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_min.y, domain_min.z);//4
	glVertex3f(domain_max.x, domain_min.y, domain_min.z);//3
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_min.y, domain_min.z);//4
	glVertex3f(domain_min.x, domain_min.y, domain_max.z);//1
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_max.x, domain_min.y, domain_max.z);//2
	glVertex3f(domain_max.x, domain_min.y, domain_min.z);//3
	glEnd();

	//ceil(5678)
	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_max.y, domain_max.z);//5
	glVertex3f(domain_max.x, domain_max.y, domain_max.z);//6
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_max.y, domain_min.z);//8
	glVertex3f(domain_max.x, domain_max.y, domain_min.z);//7
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_max.y, domain_min.z);//8
	glVertex3f(domain_min.x, domain_max.y, domain_max.z);//5
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_max.x, domain_max.y, domain_max.z);//6
	glVertex3f(domain_max.x, domain_max.y, domain_min.z);//7
	glEnd();

	//left face (14,58,15,48)
	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_min.y, domain_max.z);//1
	glVertex3f(domain_min.x, domain_min.y, domain_min.z);//4
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_max.y, domain_min.z);//8
	glVertex3f(domain_min.x, domain_max.y, domain_max.z);//5
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_min.y, domain_max.z);//1
	glVertex3f(domain_min.x, domain_max.y, domain_max.z);//5
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_max.y, domain_min.z);//8
	glVertex3f(domain_min.x, domain_min.y, domain_min.z);//4
	glEnd();

	//right face (23,67,26,37)
	glBegin(GL_LINES);
	glVertex3f(domain_max.x, domain_min.y, domain_max.z);//2
	glVertex3f(domain_max.x, domain_min.y, domain_min.z);//3
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_max.x, domain_max.y, domain_max.z);//6
	glVertex3f(domain_max.x, domain_max.y, domain_min.z);//7
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_max.x, domain_min.y, domain_max.z);//2
	glVertex3f(domain_max.x, domain_max.y, domain_max.z);//6
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_max.x, domain_min.y, domain_min.z);//3
	glVertex3f(domain_max.x, domain_max.y, domain_min.z);//7
	glEnd();

	//back face(43,78,37,48)
	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_min.y, domain_min.z);//4
	glVertex3f(domain_max.x, domain_min.y, domain_min.z);//3
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_max.y, domain_min.z);//8
	glVertex3f(domain_max.x, domain_max.y, domain_min.z);//7
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_max.x, domain_min.y, domain_min.z);//3
	glVertex3f(domain_max.x, domain_max.y, domain_min.z);//7
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_min.y, domain_min.z);//4
	glVertex3f(domain_min.x, domain_max.y, domain_min.z);//8
	glEnd();

	//front face(12,56,15,26)
	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_min.y, domain_max.z);//1
	glVertex3f(domain_max.x, domain_min.y, domain_max.z);//2
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_max.y, domain_max.z);//5
	glVertex3f(domain_max.x, domain_max.y, domain_max.z);//6
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_min.x, domain_min.y, domain_max.z);//1
	glVertex3f(domain_min.x, domain_max.y, domain_max.z);//5
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(domain_max.x, domain_min.y, domain_max.z);//2
	glVertex3f(domain_max.x, domain_max.y, domain_max.z);//6
	glEnd();
}

void ParticleSystem::DrawText() {
	char msg[100];
	for (int n = 0; n < this->num_points_; n++) {

		sprintf(msg, "%d", n);
		glColor4f((RED(clr_[n]) + 1.0)*0.5, (GRN(clr_[n]) + 1.0)*0.5, (BLUE(clr_[n]) + 1.0)*0.5, ALPH(clr_[n]));
		drawText3D(pos_[n].x, pos_[n].y, pos_[n].z, msg);
	}
}

void ParticleSystem::setExampleParams(bool bStart) {

	switch ((int)param_[PEXAMPLE]) {
	case 0:
		if (toggle_[PUSELOADEDSCENE] == true) {
			vec_[PBOUNDARYMIN].Set(-80, 0, -80);
			vec_[PBOUNDARYMAX].Set(80, 160, 80);
			vec_[PINITPARTICLEMIN].Set(-80, 0, -80);
			vec_[PINITPARTICLEMAX].Set(80, 10, 80);
		}
		else {
			vec_[PBOUNDARYMIN].Set(-30, 0, -50);
			vec_[PBOUNDARYMAX].Set(30, 90, 0);
			vec_[PINITPARTICLEMIN].Set(0, 0, -50);
			vec_[PINITPARTICLEMAX].Set(30, 40, 0);
		}
		param_[PFORCE_MIN] = 0.0;
		param_[PGROUND_SLOPE] = 0.0;
		break;

	case 1:
		if (toggle_[PUSELOADEDSCENE] == true) {
			vec_[PBOUNDARYMIN].Set(-80, 0, -80);
			vec_[PBOUNDARYMAX].Set(80, 160, 80);
			vec_[PINITPARTICLEMIN].Set(-80, 0, -80);
			vec_[PINITPARTICLEMAX].Set(80, 20, 80);
		}
		else {
			vec_[PBOUNDARYMIN].Set(-80, 0, -80);
			vec_[PBOUNDARYMAX].Set(80, 160, 80);
			vec_[PINITPARTICLEMIN].Set(-80, 0, -80);
			vec_[PINITPARTICLEMAX].Set(80, 20, 80);
		}
		param_[PFORCE_MIN] = 20.0;
		param_[PGROUND_SLOPE] = 0.10;
		break;

	case 2:
		if (toggle_[PUSELOADEDSCENE] == true) {
			vec_[PBOUNDARYMIN].Set(-80, 0, -80);
			vec_[PBOUNDARYMAX].Set(80, 160, 80);
			vec_[PINITPARTICLEMIN].Set(-80, 0, -80);
			vec_[PINITPARTICLEMAX].Set(80, 30, 80);
		}
		else {
			vec_[PBOUNDARYMIN].Set(-80, 0, -80);
			vec_[PBOUNDARYMAX].Set(80, 160, 80);
			vec_[PINITPARTICLEMIN].Set(-80, 0, -80);
			vec_[PINITPARTICLEMAX].Set(80, 60, 80);
		}
		param_[PFORCE_MIN] = 20.0;
		param_[PFORCE_MAX] = 20.0;
		vec_[PPLANE_GRAV_DIR].Set(0.0, -9.8, 0);
		break;
	}

	// 从xml文件中加载场景
	ParseXML("Scene", (int)param_[PEXAMPLE], bStart);
}

void ParticleSystem::setKernels() {

	if (param_[PSPACINGREALWORLD] == 0)
	{
		param_[PSPACINGREALWORLD] = pow((float)param_[PMASS] / (float)param_[PRESTDENSITY], 1 / 3.0f);
	}
	param_[PKERNELSELF] = KernelM4(0.0f, param_[PSMOOTHRADIUS]);
	poly6_kern_ = 315.0f / (64.0f * MY_PI * pow(param_[PSMOOTHRADIUS], 9));
	spiky_kern_ = -45.0f / (MY_PI * pow(param_[PSMOOTHRADIUS], 6));
	lap_kern_ = 45.0f / (MY_PI * pow(param_[PSMOOTHRADIUS], 6));

	float sr = param_[PSMOOTHRADIUS];
	for (int i = 0; i<lutSize; i++)
	{
		float dist = sr * i / lutSize;
		lutKernelM4[i] = KernelM4(dist, sr);
		lutKernelPressureGrad[i] = KernelPressureGrad(dist, sr);
	}

}

void ParticleSystem::setSpacing() {

	if (param_[PSPACINGGRAPHICSWORLD] == 0) {

		if (param_[PSPACINGREALWORLD] == 0)
		{
			param_[PSPACINGREALWORLD] = pow((float)param_[PMASS] / param_[PRESTDENSITY], 1 / 3.0f);
		}

		param_[PSPACINGGRAPHICSWORLD] = param_[PSPACINGREALWORLD] / param_[PSIMSCALE];
	}
	else {
		param_[PSPACINGREALWORLD] = param_[PSPACINGGRAPHICSWORLD] * param_[PSIMSCALE];
		param_[PRESTDENSITY] = param_[PMASS] / pow((float)param_[PSPACINGREALWORLD], 3.0f);
	}

	vec_[PGRIDVOLUMEMIN] = vec_[PBOUNDARYMIN] - (param_[PSMOOTHRADIUS] / param_[PSIMSCALE]);
	vec_[PGRIDVOLUMEMAX] = vec_[PBOUNDARYMAX] + (param_[PSMOOTHRADIUS] / param_[PSIMSCALE]);
}

void ParticleSystem::ParseXML(std::string name, int id, bool bStart) {

	xml.setBase(name, id);

	xml.assignValueF(&time_step_, "DT");
	xml.assignValueStr(scene_name_, "Name");
	if (bStart)	xml.assignValueF(&param_[PMAXNUM], "Num");
	xml.assignValueF(&param_[PGRID_DENSITY], "GridDensity");
	xml.assignValueF(&param_[PSIMSCALE], "SimScale");
	xml.assignValueF(&param_[PVISC], "Viscosity");
	xml.assignValueF(&param_[PRESTDENSITY], "RestDensity");
	xml.assignValueF(&param_[PSPACINGGRAPHICSWORLD], "SpacingGraphicsWorld");
	xml.assignValueF(&param_[PMASS], "Mass");
	xml.assignValueF(&param_[PCOLLISIONRADIUS], "Radius");
	xml.assignValueF(&param_[PSPACINGREALWORLD], "SearchDist");
	xml.assignValueF(&param_[PGASCONSTANT], "IntStiff");
	xml.assignValueF(&param_[PBOUNDARYSTIFF], "BoundStiff");
	xml.assignValueF(&param_[PBOUNDARYDAMP], "BoundDamp");
	xml.assignValueF(&param_[PACCEL_LIMIT], "AccelLimit");
	xml.assignValueF(&param_[PVEL_LIMIT], "VelLimit");
	xml.assignValueF(&param_[PPOINT_GRAV_AMT], "PointGravAmt");
	xml.assignValueF(&param_[PGROUND_SLOPE], "GroundSlope");
	xml.assignValueF(&param_[PFORCE_MIN], "WaveForceMin");
	xml.assignValueF(&param_[PFORCE_MAX], "WaveForceMax");
	xml.assignValueF(&param_[PFORCE_FREQ], "WaveForceFreq");
	xml.assignValueF(&param_[PDRAWMODE], "DrawMode");
	xml.assignValueF(&param_[PDRAWTEXT], "drawText2D");

	xml.assignValueV3(&vec_[PBOUNDARYMIN], "BoundaryMin");
	xml.assignValueV3(&vec_[PBOUNDARYMAX], "BoundaryMax");
	xml.assignValueV3(&vec_[PINITPARTICLEMIN], "InitParticleVolumeMin");
	xml.assignValueV3(&vec_[PINITPARTICLEMAX], "InitParticleVolumMax");
	xml.assignValueV3(&vec_[PPOINT_GRAV_POS], "PointGravPos");
	xml.assignValueV3(&vec_[PPLANE_GRAV_DIR], "PlaneGravDir");

}