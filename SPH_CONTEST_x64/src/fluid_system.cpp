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

	param_[PEXAMPLE] = 0;

	points.clear();

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

ParticleSystem::~ParticleSystem() {
	points.clear();
}

// 程序运行入口, 内存分配
void ParticleSystem::setup(bool bStart) {
	frame_ = 0;
	time_ = 0;

	setDefaultParams();

	setExampleParams(bStart);

	param_[PGRIDSIZEREALSCALE] = param_[PSMOOTHRADIUS] / param_[PGRID_DENSITY];

	setKernels();

	setSpacing();

	ComputeGasConstAndTimeStep(param_[PMAXDENSITYERRORALLOWED]);

	ClearNeighborTable();

	float lengthX = vec_[PINITPARTICLEMAX].x - vec_[PINITPARTICLEMIN].x;
	float lengthY = vec_[PINITPARTICLEMAX].y - vec_[PINITPARTICLEMIN].y;
	float lengthZ = vec_[PINITPARTICLEMAX].z - vec_[PINITPARTICLEMIN].z;
	Vector3DF lengthXYZ(lengthX, lengthY, lengthZ);

	int numParticlesX = ceil(lengthX / param_[PSPACINGGRAPHICSWORLD]);
	int numParticlesY = ceil(lengthY / param_[PSPACINGGRAPHICSWORLD]);
	int numParticlesZ = ceil(lengthZ / param_[PSPACINGGRAPHICSWORLD]);
	int numParticles = numParticlesX * numParticlesY * numParticlesZ;
	Vector3DI numParticlesXYZ(numParticlesX, numParticlesY, numParticlesZ);

	AllocateParticlesMemory(param_[PMAXNUM]);

	// 从文件中加载粒子位置信息,并设置其速度，颜色等属性值
	if (toggle_[PUSELOADEDSCENE] == true) {

		// 从文件中读取bunny模型数据
		const char* file_name = "models/bunny.txt";

		Vector3DF minVec;
		Vector3DF maxVec;
		num_points_ = ReadInFluidParticles(file_name, minVec, maxVec);
		
		setInitParticleVolumeFromFile(minVec, maxVec);
		setInitParticleVolume(numParticlesXYZ, lengthXYZ, 0.1);
	}
	else {
		num_points_ = 0;
		setInitParticleVolume(numParticlesXYZ, lengthXYZ, 0.1);
	}

	AllocatePackBuf();

	setGridAllocate(1.0);
}

void ParticleSystem::setRender() {

	glEnable(GL_TEXTURE_2D);

	glGenTextures(1, (GLuint*)texture_);
	glBindTexture(GL_TEXTURE_2D, texture_[0]);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, 8, 8, 0, GL_RGB, GL_FLOAT, 0);

	glGenBuffers(3, (GLuint*)vbo_);

	const int udiv = 6;
	const int vdiv = 6;
	const float du = 180.0 / udiv;
	const float dv = 360.0 / vdiv;
	const float r = 1.0;

	Vector3DF* buf = (Vector3DF*)malloc(sizeof(Vector3DF) * (udiv + 2)*(vdiv + 2) * 2);
	Vector3DF* dat = buf;

	sphere_points_ = 0;
	for (float tilt = -90; tilt <= 90.0; tilt += du) {
		for (float ang = 0; ang <= 360; ang += dv) {
			float x = sin(ang*DEGtoRAD) * cos(tilt*DEGtoRAD);
			float y = cos(ang*DEGtoRAD) * cos(tilt*DEGtoRAD);
			float z = sin(tilt*DEGtoRAD);
			float x1 = sin(ang*DEGtoRAD) * cos((tilt + du)*DEGtoRAD);
			float y1 = cos(ang*DEGtoRAD) * cos((tilt + du)*DEGtoRAD);
			float z1 = sin((tilt + du)*DEGtoRAD);

			dat->x = x*r;
			dat->y = y*r;
			dat->z = z*r;
			dat++;
			dat->x = x1*r;
			dat->y = y1*r;
			dat->z = z1*r;
			dat++;
			sphere_points_ += 2;
		}
	}
	glBindBuffer(GL_ARRAY_BUFFER, vbo_[2]);
	glBufferData(GL_ARRAY_BUFFER, sphere_points_ * sizeof(Vector3DF), buf, GL_STATIC_DRAW);
	glVertexPointer(3, GL_FLOAT, 0, 0x0);

	free(buf);

	image_.read("ball32.bmp", "ball32a.bmp");
}

int ParticleSystem::SelectParticle(int x, int y, int wx, int wy, Camera3D& cam) {

	Vector4DF pnt;

	for (int n = 0; n < num_points_; n++) {
		pnt = cam.project(points[n].pos);
		pnt.x = (pnt.x + 1.0)*0.5 * wx;
		pnt.y = (pnt.y + 1.0)*0.5 * wy;

		if (x > pnt.x - 8 && x < pnt.x + 8 && y > pnt.y - 8 && y < pnt.y + 8) {
			selected_ = n;
			return n;
		}
	}
	selected_ = -1;
	return -1;
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

	mint::Time start;
	
	// 插入粒子
	start.SetSystemTime(ACC_NSEC);
	InsertParticlesCPU();
	Record(PTIME_INSERT, "Insert CPU SPH", start);

	// 计算压力
	start.SetSystemTime(ACC_NSEC);
	ComputePressureGrid();
	Record(PTIME_PRESS, "Press CPU SPH", start);

	// 计算外力
	start.SetSystemTime(ACC_NSEC);
	ComputeForceGrid();
	Record(PTIME_FORCE, "Force CPU SPH", start);

	// 移动粒子
	start.SetSystemTime(ACC_NSEC);
	AdvanceStepSimpleCollision();
	Record(PTIME_ADVANCE, "Advance CPU SPH", start);
}

void ParticleSystem::InsertParticlesCPU() {

	memset(grid_head_cell_particle_index_array_, GRID_UNDEF, grid_total_ * sizeof(uint));
	memset(grid_particles_number_,               0,          grid_total_ * sizeof(uint));

	const int xns = grid_res_.x;
	const int yns = grid_res_.y;
	const int zns = grid_res_.z;

	param_[PSTAT_OCCUPANCY] = 0.0; // 有粒子的grid的数量
	param_[PSTAT_GRIDCOUNT] = 0.0; // grid中粒子的总数

	for (int idx = 0; idx < num_points_; ++idx) {
		Vector3DI gridCell;
		const int gridCellIndex = getGridCell(points[idx].pos, gridCell);

		if (gridCell.x >= 0 && gridCell.x < xns && gridCell.y >= 0 && gridCell.y < yns && gridCell.z >= 0 && gridCell.z < zns) {
			points[idx].particle_grid_cell_index = gridCellIndex;
			points[idx].next_particle_index_in_the_same_cell = grid_head_cell_particle_index_array_[gridCellIndex];
			if (points[idx].next_particle_index_in_the_same_cell == GRID_UNDEF) {
				param_[PSTAT_OCCUPANCY] += 1.0;
			}
			grid_head_cell_particle_index_array_[gridCellIndex] = idx;
			grid_particles_number_[gridCellIndex] += 1;
			param_[PSTAT_GRIDCOUNT] += 1.0;
		}
		else {
			Vector3DF vel = points[idx].vel;
			Vector3DF ve = points[idx].vel_eval;
			float pr = points[idx].pressure;
			float dn = points[idx].density;
			printf("WARNING: ParticleSystem::InsertParticlesCPU(): Out of Bounds: %d, Position<%f %f %f>, Velocity<%f %f %f>, Pressure:%f, Density:%f\n", idx, points[idx].pos.x, points[idx].pos.y, points[idx].pos.z, vel.x, vel.y, vel.z, pr, dn);
			points[idx].pos.x = -1; points[idx].pos.y = -1; points[idx].pos.z = -1;
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
	for (int i = 0; i < num_points_; ++i) {

		points[i].density = own_density_contribution;

		int neighbor_nums = 0;
		int search_nums = 0;
		float sum = 0.0;

		const uint i_cell_index = points[i].particle_grid_cell_index;
		if (i_cell_index != GRID_UNDEF) {
			for (int cell = 0; cell < max_num_adj_grid_cells_cpu; cell++) {
				const int neighbor_cell_index = i_cell_index + grid_neighbor_cell_index_offset_[cell];
				if (neighbor_cell_index == GRID_UNDEF || neighbor_cell_index < 0 || neighbor_cell_index > grid_total_ - 1) {
					continue;
				}

				int j = grid_head_cell_particle_index_array_[neighbor_cell_index];
				while (j != GRID_UNDEF) {
					if (i == j) {
						j = points[j].next_particle_index_in_the_same_cell;
						continue;
					}
					Vector3DF dst_graphics_scale = points[i].pos;
					dst_graphics_scale -= points[j].pos;
					const float dist_square_sim_scale = sim_scale_square*(dst_graphics_scale.x*dst_graphics_scale.x + dst_graphics_scale.y*dst_graphics_scale.y + dst_graphics_scale.z*dst_graphics_scale.z);
					if (dist_square_sim_scale <= smooth_radius_square) {
						const float dist = sqrt(dist_square_sim_scale);
						float kernelValue = KernelM4Lut(dist, smooth_radius);
						points[i].density += kernelValue * mass;

						neighbor_nums++;
					}
					search_nums++;
					j = points[j].next_particle_index_in_the_same_cell;
				}
			}
		}

		if (points[i].density < minDens)
			minDens = points[i].density;
		if (points[i].density > maxDens)
			maxDens = points[i].density;

		points[i].pressure = max(0.0f, (points[i].density - param_[PRESTDENSITY]) * param_[PGASCONSTANT]);

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
		points[i].force.Set(0, 0, 0);
		Vector3DF force(0, 0, 0);
		const uint i_cell_index = points[i].particle_grid_cell_index;
		Vector3DF ipos = points[i].pos;
		Vector3DF iveleval = points[i].vel_eval;
		float	  ipress = points[i].pressure;
		float	  idensity = points[i].density;

		if (i_cell_index != GRID_UNDEF) {
			for (int cell = 0; cell < max_num_adj_grid_cells_cpu; cell++) {
				const int neighbor_cell_index = i_cell_index + grid_neighbor_cell_index_offset_[cell];
				if (neighbor_cell_index == GRID_UNDEF || (neighbor_cell_index < 0 || neighbor_cell_index > grid_total_ - 1)) {
					continue;
				}

				int j = grid_head_cell_particle_index_array_[neighbor_cell_index];
				while (j != GRID_UNDEF) {
					if (i == j) {
						j = points[j].next_particle_index_in_the_same_cell;
						continue;
					}

					Vector3DF vector_i_minus_j = ipos - points[j].pos;
					const float dx = vector_i_minus_j.x;
					const float dy = vector_i_minus_j.y;
					const float dz = vector_i_minus_j.z;

					const float dist_square_sim_scale = sim_scale_square*(dx*dx + dy*dy + dz*dz);
					if (dist_square_sim_scale <= smooth_radius_square && dist_square_sim_scale > 0) {
						const float jdist = sqrt(dist_square_sim_scale);
						const float jpress = points[j].pressure;
						const float h_minus_r = smooth_radius - jdist;
						const float pterm = -0.5f * h_minus_r * spiky_kern_ * (ipress + jpress) / jdist;
						const float dterm = h_minus_r / (idensity * points[j].density);

						Vector3DF vel_j_minus_i = points[j].vel_eval;
						vel_j_minus_i -= iveleval;

						force += vector_i_minus_j * sim_scale * pterm * dterm;

						force += vel_j_minus_i * vterm * dterm;
					}
					j = points[j].next_particle_index_in_the_same_cell;
				}
			}
		}

		force *= mass * mass;
		force += vec_gravity * mass;

		if (addBoundaryForce) {
			force += BoxBoundaryForce(i);
		}

		points[i].force = force;

	}
}

void ParticleSystem::AdvanceStepSimpleCollision() {

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

	for (int i = 0; i < num_points_; i++) {
		if (points[i].particle_grid_cell_index == GRID_UNDEF)
			continue;

		Vector3DF acceleration = points[i].force;// + points[i].correction_pressure_force;
		acceleration /= param_[PMASS];

		points[i].vel_eval += acceleration * time_step_;
		points[i].pos *= sim_scale;
		points[i].pos += points[i].vel_eval * time_step_;

		CollisionHandlingSimScale(points[i].pos, points[i].vel_eval);

		points[i].pos /= sim_scale;
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
	return getGridCell(points[p].pos, gc);
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

	if (points[i].pos.y < bound_min.y + forceDistance) {
		distToWall = bound_min.y + forceDistance - points[i].pos.y;
		factor = distToWall * invForceDist;
		force += Vector3DF(0, 1, 0) * factor * 2.0f * forceStrength;
	}
	else if (points[i].pos.y > bound_max.y - forceDistance) {
		distToWall = points[i].pos.y - (bound_max.y - forceDistance);
		factor = distToWall * invForceDist;
		force += Vector3DF(0, -1, 0) * factor * forceStrength;
	}

	if (points[i].pos.x < bound_min.x + forceDistance) {
		distToWall = bound_min.x + forceDistance - points[i].pos.x;
		factor = distToWall * invForceDist;
		force += Vector3DF(1, 0, 0) * factor * forceStrength;
	}
	if (points[i].pos.x > bound_max.x - forceDistance) {
		distToWall = points[i].pos.x - (bound_max.x - forceDistance);
		factor = distToWall * invForceDist;
		force += Vector3DF(-1, 0, 0) * 1 * factor * forceStrength;
	}

	if (points[i].pos.z < bound_min.z + forceDistance) {
		distToWall = bound_min.z + forceDistance - points[i].pos.z;
		factor = distToWall * invForceDist;
		force += Vector3DF(0, 0, 1) * factor * forceStrength;
	}
	if (points[i].pos.z > bound_max.z - forceDistance) {
		distToWall = points[i].pos.z - (bound_max.z - forceDistance);
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
	sprintf(disp, "Grid Cell:    <%d, %d, %d> id: %d", gc.x, gc.y, gc.z, gs);
	drawText2D(10, 40, disp);

	sprintf(disp, "Neighbors:    ");

	int cnt = *(neighbor_particle_numbers_ + p);
	int ndx = *(neighbor_index_ + p);
	for (int n = 0; n < cnt; n++) {
		sprintf(disp, "%s%d, ", disp, neighbor_table_[ndx]);
		ndx++;
	}
	drawText2D(10, 70, disp);
}

// 绘制函数
void ParticleSystem::Draw(Camera3D& cam, float rad) {

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

	glEnable(GL_LIGHTING);

	for (int n = 0; n < this->num_points_; n++) {
		glPushMatrix();
		glTranslatef(points[n].pos.x, points[n].pos.y, points[n].pos.z);
		glScalef(rad, rad, rad);
		glColor4f(RED(points[n].clr), GRN(points[n].clr), BLUE(points[n].clr), ALPH(points[n].clr));
		drawSphere();
		glPopMatrix();
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
		glColor4f((RED(points[n].clr) + 1.0)*0.5, (GRN(points[n].clr) + 1.0)*0.5, (BLUE(points[n].clr) + 1.0)*0.5, ALPH(points[n].clr));
		drawText3D(points[n].pos.x, points[n].pos.y, points[n].pos.z, msg);
	}
}

void ParticleSystem::setDefaultParams() {

	param_[PMAXNUM] = 1048576;
	param_[PSIMSCALE] = 0.005;
	param_[PGRID_DENSITY] = 1.0;
	param_[PVISC] = 1.002; //4.0;			
	param_[PRESTDENSITY] = 1000.0;
	param_[PMASS] = 0.001953125;
	param_[PCOLLISIONRADIUS] = 0.00775438;
	param_[PSPACINGREALWORLD] = 0.0125;
	param_[PSMOOTHRADIUS] = 0.025025;
	param_[PGRIDSIZEREALSCALE] = param_[PSMOOTHRADIUS] / param_[PGRID_DENSITY];
	param_[PGASCONSTANT] = 1000.0;
	param_[PBOUNDARYSTIFF] = 2.5;
	param_[PBOUNDARYDAMP] = 1.0;
	param_[PACCEL_LIMIT] = 150.0;
	param_[PVEL_LIMIT] = 3.0;
	param_[PSPACINGGRAPHICSWORLD] = 0.0;
	param_[PGROUND_SLOPE] = 0.0;
	param_[PFORCE_MIN] = 0.0;
	param_[PFORCE_MAX] = 0.0;
	param_[PDRAWTEXT] = 0;
	param_[PPOINT_GRAV_AMT] = 0.0;
	param_[PSTAT_NEIGHCNTMAX] = 0;
	param_[PSTAT_SEARCHCNTMAX] = 0;
	param_[PFORCE_FREQ] = 8.0;
	param_[PMINLOOPPCISPH] = 3;
	param_[PMAXLOOPPCISPH] = MAX_PCISPH_LOOPS;
	param_[PMAXDENSITYERRORALLOWED] = 5.0;

	vec_[PEMIT_POS].Set(0, 0, 0);
	vec_[PEMIT_ANG].Set(0, 90, 1.0);
	vec_[PEMIT_RATE].Set(0, 0, 0);
	vec_[PPOINT_GRAV_POS].Set(0, 50, 0);
	vec_[PPLANE_GRAV_DIR].Set(0, -9.8, 0);

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

	if (param_[PSPACINGREALWORLD] == 0) {
		param_[PSPACINGREALWORLD] = pow((float)param_[PMASS] / (float)param_[PRESTDENSITY], 1 / 3.0f);
	}
	param_[PKERNELSELF] = KernelM4(0.0f, param_[PSMOOTHRADIUS]);
	poly6_kern_ = 315.0f / (64.0f * MY_PI * pow(param_[PSMOOTHRADIUS], 9));
	spiky_kern_ = -45.0f / (MY_PI * pow(param_[PSMOOTHRADIUS], 6));
	lap_kern_ = 45.0f / (MY_PI * pow(param_[PSMOOTHRADIUS], 6));

	float sr = param_[PSMOOTHRADIUS];
	for (int i = 0; i<lutSize; i++) {
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

void ParticleSystem::setInitParticleVolumeFromFile(const Vector3DF& minVec, const Vector3DF& maxVec) {

	const float lengthX = maxVec.x - minVec.x;
	const float lengthY = maxVec.y - minVec.y;
	const float lengthZ = maxVec.z - minVec.z;
	const float inv_sim_scale = 1.0f / param_[PSIMSCALE];

	for (int i = 0; i < num_points_; ++i) {
		points[i].pos *= inv_sim_scale;
		points[i].vel_eval.Set(0.0, 0.0, 0.0);
		points[i].clr = COLORA((points[i].pos.x - minVec.x) / lengthX, (points[i].pos.y - minVec.y) / lengthY, (points[i].pos.z - minVec.z) / lengthZ, 1);
	}
}

void ParticleSystem::setInitParticleVolume(const Vector3DI& numParticlesXYZ, const Vector3DF& lengthXYZ, const float jitter) {

	srand(time(0x0));

	float spacingRealWorldSize = param_[PSPACINGREALWORLD];
	float particleVolumeRealWorldSize = spacingRealWorldSize * spacingRealWorldSize *spacingRealWorldSize;
	float mass = param_[PRESTDENSITY] * particleVolumeRealWorldSize;

	float tmpX, tmpY, tmpZ;
	if (numParticlesXYZ.x % 2 == 0) tmpX = 0.0;
	else                            tmpX = 0.5;
	if (numParticlesXYZ.y % 2 == 0) tmpY = 0.0;
	else                            tmpY = 0.5;
	if (numParticlesXYZ.z % 2 == 0) tmpZ = 0.0;
	else                            tmpZ = 0.5;

	int i = num_points_;
	for (int iy = 0; iy < numParticlesXYZ.y; iy++) {

		float y = 0.0 + (iy + tmpY) * param_[PSPACINGGRAPHICSWORLD];

		for (int ix = 0; ix < numParticlesXYZ.x; ix++) {

			float x = vec_[PINITPARTICLEMIN].x + (ix + tmpX) * param_[PSPACINGGRAPHICSWORLD];

			for (int iz = 0; iz < numParticlesXYZ.z; iz++) {

				float z = vec_[PINITPARTICLEMIN].z + (iz + tmpZ) * param_[PSPACINGGRAPHICSWORLD];

				if (num_points_ < param_[PMAXNUM]) {

					points[i].pos.Set(x + (frand() - 0.5) * jitter, y + (frand() - 0.5) * jitter, z + (frand() - 0.5) * jitter);
					points[i].vel_eval.Set(0.0, 0.0, 0.0);
					points[i].clr = COLORA((x - vec_[PINITPARTICLEMIN].x) / lengthXYZ.x, (y - vec_[PINITPARTICLEMIN].y) / lengthXYZ.y, (z - vec_[PINITPARTICLEMIN].z) / lengthXYZ.z, 1);
					num_points_++;
				}
				++i;
			}
		}
	}
}

void ParticleSystem::setGridAllocate(const float border) {

	float world_cellsize = param_[PGRIDSIZEREALSCALE] / param_[PSIMSCALE];

	grid_min_ = vec_[PGRIDVOLUMEMIN];
	grid_max_ = vec_[PGRIDVOLUMEMAX];

	grid_size_ = vec_[PGRIDVOLUMEMAX] - vec_[PGRIDVOLUMEMIN];
	grid_res_.x = (int)ceil(grid_size_.x / world_cellsize);
	grid_res_.y = (int)ceil(grid_size_.y / world_cellsize);
	grid_res_.z = (int)ceil(grid_size_.z / world_cellsize);

	if (grid_res_.x == 0)
		grid_res_.x = 1;
	if (grid_res_.y == 0)
		grid_res_.y = 1;
	if (grid_res_.z == 0)
		grid_res_.z = 1;

	grid_total_ = (int)(grid_res_.x * grid_res_.y * grid_res_.z);
	if (grid_total_ > 10000000) {
		printf("ERROR: ParticleSystem::setGridAllocate(...): number: %d, too many cells in initGrid().\n", grid_total_);
	}

	grid_size_.x = grid_res_.x * world_cellsize;
	grid_size_.y = grid_res_.y * world_cellsize;
	grid_size_.z = grid_res_.z * world_cellsize;
	grid_delta_ = grid_res_;
	grid_delta_ /= grid_size_;

	if (grid_head_cell_particle_index_array_ != 0x0)
		free(grid_head_cell_particle_index_array_);
	if (grid_particles_number_ != 0x0)
		free(grid_particles_number_);
	grid_head_cell_particle_index_array_ = (uint*)malloc(sizeof(uint*) * grid_total_);
	grid_particles_number_ = (uint*)malloc(sizeof(uint*) * grid_total_);
	memset(grid_head_cell_particle_index_array_, GRID_UNDEF, grid_total_ * sizeof(uint));
	memset(grid_particles_number_, GRID_UNDEF, grid_total_ * sizeof(uint));

	param_[PSTAT_GMEM] = 12 * grid_total_;

	grid_search_ = 3;
	int cell = 0;
	for (int y = -1; y < 2; y++)
		for (int z = -1; z < 2; z++)
			for (int x = -1; x < 2; x++)
				grid_neighbor_cell_index_offset_[cell++] = y * grid_res_.x *grid_res_.z + z * grid_res_.x + x;

	grid_adj_cnt_ = grid_search_ * grid_search_ * grid_search_;

	if (pack_grid_buf_ != 0x0)
		free(pack_grid_buf_);
	pack_grid_buf_ = (int*)malloc(sizeof(int) * grid_total_);
}

void ParticleSystem::AllocatePackBuf() {

	if (pack_fluid_particle_buf_ != 0x0)
		free(pack_fluid_particle_buf_);
	pack_fluid_particle_buf_ = (char*)malloc(sizeof(Particle) * param_[PMAXNUM]);
}

// 分配内存
void ParticleSystem::AllocateParticlesMemory(int cnt) {
	if (!points.empty())
		points.clear();
	points.resize(cnt);

	if (neighbor_index_ != 0x0)
		free(neighbor_index_);
	neighbor_index_ = (uint*)malloc(sizeof(uint) * param_[PMAXNUM]);

	if (neighbor_particle_numbers_ != 0x0)
		free(neighbor_particle_numbers_);
	neighbor_particle_numbers_ = (uint*)malloc(sizeof(uint) * param_[PMAXNUM]);

	param_[PSTAT_PMEM] = sizeof(Particle) * cnt;
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
	xml.assignValueF(&param_[PDRAWTEXT], "drawText2D");

	xml.assignValueV3(&vec_[PBOUNDARYMIN], "BoundaryMin");
	xml.assignValueV3(&vec_[PBOUNDARYMAX], "BoundaryMax");
	xml.assignValueV3(&vec_[PINITPARTICLEMIN], "InitParticleVolumeMin");
	xml.assignValueV3(&vec_[PINITPARTICLEMAX], "InitParticleVolumMax");
	xml.assignValueV3(&vec_[PPOINT_GRAV_POS], "PointGravPos");
	xml.assignValueV3(&vec_[PPLANE_GRAV_DIR], "PlaneGravDir");

}

void ParticleSystem::Record(int param, std::string name, mint::Time& start) {

	mint::Time stop;
	stop.SetSystemTime(ACC_NSEC);
	stop = stop - start;
	param_[param] = stop.GetMSec();
	if (toggle_[PPROFILE]) printf("%s:  %s\n", name.c_str(), stop.GetReadableTime().c_str());
}

int ParticleSystem::ReadInFluidParticles(const char * filename, Vector3DF & minVec, Vector3DF & maxVec) {

	float px, py, pz;
	float min_x, min_y, min_z;
	float max_x, max_y, max_z;
	int   cnt, readCnt = 0;
	try {
		infileParticles.open(filename);

		// 文件第一行为粒子数量
		infileParticles >> cnt;

		infileParticles >> min_x >> min_y >> min_z;
		infileParticles >> max_x >> max_y >> max_z;

		minVec.Set(min_x, min_y, min_z);
		maxVec.Set(max_x, max_y, max_z);

		while (infileParticles && readCnt < cnt) {

			infileParticles >> px >> py >> pz;
			points[readCnt].pos.Set(px, py, pz);
			readCnt++;
		}
	}
	catch (...) {
		printf("ERROR: ParticleSystem::ReadInFluidParticles(...): File %s open failed.\n", filename);
	}

	infileParticles.close();
	return cnt;
}
