#include "fluid_system.h"

ParticleSystem::ParticleSystem() {
	frame = 0;
	time_ = 0;

	// 程序运行方法---3种方法可选
	param_[PRUN_MODE] = RUN_CPU_SPH;

	points_number = 0;

	param_[PEXAMPLE] = 0;

	param_[PDRAWMODE] = 0;

	points.clear();

	grid_total = 0;

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
void ParticleSystem::setUp() {
	frame = 0;
	time_ = 0;

	setDefaultParams();

	setExampleParams();

	param_[PGRIDSIZEREALSCALE] = param_[PSMOOTHRADIUS] / param_[PGRID_DENSITY];

	setKernels();

	setSpacing();

	ComputeGasConstAndTimeStep();

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

	points_number = 0;
	switch ((int)param_[PEXAMPLE]) {
		case 0:
			setInitParticleVolume(numParticlesXYZ, lengthXYZ, 0.1);
			break;
		case 1:
			setInitParticleVolumeNew(numParticlesXYZ, lengthXYZ, 0.1);
			break;
	}
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

	sphere_points = 0;
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
			sphere_points += 2;
		}
	}
	glBindBuffer(GL_ARRAY_BUFFER, vbo_[2]);
	glBufferData(GL_ARRAY_BUFFER, sphere_points * sizeof(Vector3DF), buf, GL_STATIC_DRAW);
	glVertexPointer(3, GL_FLOAT, 0, 0x0);

	free(buf);

	image.read("ball32.bmp", "ball32a.bmp");
}

int ParticleSystem::SelectParticle(int x, int y, int wx, int wy, Camera3D& cam) {

	Vector4DF pnt;

	for (int n = 0; n < points_number; n++) {
		pnt = cam.project(points[n].pos);
		pnt.x = (pnt.x + 1.0)*0.5 * wx;
		pnt.y = (pnt.y + 1.0)*0.5 * wy;

		if (x > pnt.x - 8 && x < pnt.x + 8 && y > pnt.y - 8 && y < pnt.y + 8) {
			return n;
		}
	}
	return -1;
}

std::string ParticleSystem::getModeStr() {

	std::string buf = "SIMULATE CPU SPH";

	switch ((int)param_[PRUN_MODE]) {
	case RUN_CPU_SPH:
		buf = "SIMULATE CPU SPH";
		break;
	case RUN_CPU_PBF:
		buf = "SIMULATE CPU PBF";
		break;
	case RUN_CPU_XLZ:
		buf = "SIMULATE CPU XLZ";
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
	param_[PTIME_PBF_STEP] = 0.0;

	// 运行程序
	switch ((int)param_[PRUN_MODE]) {
	case RUN_CPU_SPH:
		RunCPUSPH();
		break;
	case RUN_CPU_PBF:
		RunCPUPBF();
		break;
	case RUN_CPU_XLZ:
		RunCPUXLZ();
		break;
	}

	time_ += time_step;
	frame++;

	if (frame == 1001) {
		system("pause");
	}
	else if (frame == 2001) {
		system("pause");
	}
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

void ParticleSystem::RunCPUPBF() {
	mint::Time start;

	// 插入粒子
	start.SetSystemTime(ACC_NSEC);
	InsertParticlesCPU();
	Record(PTIME_INSERT, "Insert CPU SPH", start);

	// 计算压力和密度
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

	// PBF
	start.SetSystemTime(ACC_NSEC);
	PositionBasedFluid();
	Record(PTIME_PBF_STEP, "PBF Step Time", start);
}

void ParticleSystem::RunCPUXLZ()
{
	mint::Time start;

	// 插入粒子
	start.SetSystemTime(ACC_NSEC);
	InsertParticlesCPU();
	Record(PTIME_INSERT, "Insert CPU SPH", start);

	// 计算压力和密度
	start.SetSystemTime(ACC_NSEC);
	ComputePressureGrid();
	Record(PTIME_PRESS, "Press CPU SPH", start);

	// 计算外力
	start.SetSystemTime(ACC_NSEC);
	ComputeNewForceGrid();
	Record(PTIME_FORCE, "Force CPU SPH", start);

	// 移动粒子
	start.SetSystemTime(ACC_NSEC);
	AdvanceStepSimpleCollision();
	Record(PTIME_ADVANCE, "Advance CPU SPH", start);
}

void ParticleSystem::InsertParticlesCPU() {

	grid.clear();
	grid.resize(grid_total, Cell());

	const int xns = grid_res.x;
	const int yns = grid_res.y;
	const int zns = grid_res.z;

	grid_nonempty_cell_number = 0;
	grid_particle_number = 0;

	for (int idx = 0; idx < points_number; ++idx) {
		Vector3DI gridCell;
		const int gridCellIndex = getGridCell(points[idx].pos, gridCell);

		if (gridCell.x >= 0 && gridCell.x < xns && gridCell.y >= 0 && gridCell.y < yns && gridCell.z >= 0 && gridCell.z < zns) {
			points[idx].particle_grid_cell_index = gridCellIndex;
			points[idx].next_particle_index_in_the_same_cell = grid[gridCellIndex].first_particle_index;
			if (points[idx].next_particle_index_in_the_same_cell == UNDEFINE) {
				grid_nonempty_cell_number += 1;
			}
			grid[gridCellIndex].first_particle_index = idx;
			grid[gridCellIndex].particle_number += 1;
			grid_particle_number += 1;
		}
		else {
			Vector3DF vel = points[idx].vel;
			Vector3DF ve = points[idx].vel_eval;
			float pr = points[idx].pressure;
			float dn = points[idx].density;
			printf("WARNING: ParticleSystem::InsertParticlesCPU(): Out of Bounds: %d, Position<%f %f %f>, Velocity<%f %f %f>, Pressure:%f, Density:%f\t", idx, points[idx].pos.x, points[idx].pos.y, points[idx].pos.z, vel.x, vel.y, vel.z, pr, dn);
			printf("Grid cell x: %d, y: %d, z: %d\n", gridCell.x, gridCell.y, gridCell.z);
			points[idx].pos.x = -1; points[idx].pos.y = -1; points[idx].pos.z = -1;
		}
	}
}

void ParticleSystem::ComputePressureGrid() {

	const float	simScaleSquare = param_[PSIMSCALE] * param_[PSIMSCALE];
	const float smoothRadiusSquare = param_[PSMOOTHRADIUS] * param_[PSMOOTHRADIUS];

	for (int i = 0; i < points_number; ++i) {

		points[i].density = 0.0;

		float sum = 0.0;

		const uint iCellIndex = points[i].particle_grid_cell_index;
		if (iCellIndex != UNDEFINE) {
			for (int cell = 0; cell < grid_adj_cnt; cell++) {
				const int neighborCellIndex = iCellIndex + grid_neighbor_cell_index_offset_[cell];
				if (neighborCellIndex == UNDEFINE || neighborCellIndex < 0 || neighborCellIndex > grid_total - 1) {
					continue;
				}

				int j = grid[neighborCellIndex].first_particle_index;
				while (j != UNDEFINE) {
					if (i == j) {
						j = points[j].next_particle_index_in_the_same_cell;
						continue;
					}
					Vector3DF ri = points[i].pos;
					Vector3DF ri_rj = ri - points[j].pos;
					const float ri_rjSquare = simScaleSquare * (ri_rj.x*ri_rj.x + ri_rj.y*ri_rj.y + ri_rj.z*ri_rj.z);
					if (ri_rjSquare <= smoothRadiusSquare) {
						points[i].density += points[j].mass * pow(smoothRadiusSquare - ri_rjSquare, 3.0f);
					}
					j = points[j].next_particle_index_in_the_same_cell;
				}
			}
		}
		points[i].density *= poly6_kern;
		points[i].pressure = max(0.0f, (points[i].density - param_[PRESTDENSITY]) * param_[PGASCONSTANT]);
	}
}

void ParticleSystem::ComputeForceGrid() {

	const float simScale = param_[PSIMSCALE];
	const float simScaleSquare = simScale * simScale;
	const float smoothRadius = param_[PSMOOTHRADIUS];
	const float smoothRadiusSquare = smoothRadius * smoothRadius;
	const float visc = param_[PVISC];
	Vector3DF   gravity = vec_[PPLANE_GRAV_DIR];
	const float vterm = lap_kern * visc;
	const float eps = 1.0e-8;

	for (int i = 0; i < points_number; i++) {
		points[i].force.Set(0, 0, 0);
		Vector3DF force(0, 0, 0);
		Vector3DF ipos = points[i].pos;
		Vector3DF iveleval = points[i].vel_eval;
		float	  ipress = points[i].pressure;
		float	  idensity = points[i].density;

		const uint iCellIndex = points[i].particle_grid_cell_index;
		if (iCellIndex != UNDEFINE) {
			for (int cell = 0; cell < grid_adj_cnt; cell++) {
				const int neighbor_cell_index = iCellIndex + grid_neighbor_cell_index_offset_[cell];
				if (neighbor_cell_index == UNDEFINE || (neighbor_cell_index < 0 || neighbor_cell_index > grid_total - 1)) {
					continue;
				}

				int j = grid[neighbor_cell_index].first_particle_index;
				while (j != UNDEFINE) {
					if (i == j) {
						j = points[j].next_particle_index_in_the_same_cell;
						continue;
					}

					Vector3DF ri_rj = ipos - points[j].pos;
					const float ri_rjSquare = simScaleSquare*(ri_rj.x*ri_rj.x + ri_rj.y*ri_rj.y + ri_rj.z*ri_rj.z);
					if (ri_rjSquare <= smoothRadiusSquare) {
						//
						const float jdist = sqrt(ri_rjSquare);
						const float jpress = points[j].pressure;
						const float h_r = smoothRadius - jdist;
						const float pterm = -0.5f * h_r * spiky_kern * (ipress + jpress) / (jdist + eps);
						const float dterm = h_r / (idensity * points[j].density);

						Vector3DF uj_ui = points[j].vel_eval;
						uj_ui -= iveleval;

						force += ri_rj * simScale * pterm * dterm  * points[j].mass;

						force += uj_ui * vterm * dterm * points[j].mass;
					}
					j = points[j].next_particle_index_in_the_same_cell;
				}
			}
		}

		force *= points[i].mass;
		force += gravity * points[i].mass;

		points[i].force = force;

	}
}

void ParticleSystem::ComputeNewForceGrid() {

	const float simScale = param_[PSIMSCALE];
	const float simScaleSquare = simScale * simScale;
	const float smoothRadius = param_[PSMOOTHRADIUS];
	const float smoothRadiusSquare = smoothRadius * smoothRadius;
	Vector3DF   gravity = vec_[PPLANE_GRAV_DIR];
	const float eps = 1.0e-8;

	for (int i = 0; i < points_number; i++) {
		points[i].force.Set(0.0f, 0.0f, 0.0f);
		Vector3DF pforce(0.0f, 0.0f, 0.0f);
		Vector3DF vforce(0.0f, 0.0f, 0.0f);
		Vector3DF force(0.0f, 0.0f, 0.0f);
		Vector3DF ipos = points[i].pos;
		Vector3DF iveleval = points[i].vel_eval;
		float	  ipress = points[i].pressure;
		float	  idensity = points[i].density;
		Vector3DF D(0.0f, 0.0f, 0.0f);

		const uint iCellIndex = points[i].particle_grid_cell_index;
		if (iCellIndex != UNDEFINE) {
			for (int cell = 0; cell < grid_adj_cnt; cell++) {
				const int neighbor_cell_index = iCellIndex + grid_neighbor_cell_index_offset_[cell];
				if (neighbor_cell_index == UNDEFINE || (neighbor_cell_index < 0 || neighbor_cell_index > grid_total - 1)) {
					continue;
				}

				int j = grid[neighbor_cell_index].first_particle_index;
				while (j != UNDEFINE) {
					if (i == j) {
						j = points[j].next_particle_index_in_the_same_cell;
						continue;
					}

					Vector3DF ri_rj = ipos - points[j].pos;
					const float ri_rjSquare = simScaleSquare*(ri_rj.x*ri_rj.x + ri_rj.y*ri_rj.y + ri_rj.z*ri_rj.z);
					if (ri_rjSquare <= smoothRadiusSquare) {
						//
						const float jdist = sqrt(ri_rjSquare);
						const float jpress = points[j].pressure;
						const float h_r = smoothRadius - jdist;
						const float pterm = -0.5f * h_r * spiky_kern * (ipress + jpress) / (jdist + eps);
						const float dterm = h_r / (idensity * points[j].density);

						Vector3DF uj_ui = points[j].vel_eval;
						uj_ui -= iveleval;

						pforce += ri_rj * simScale * pterm * dterm  * points[j].mass;

						vforce += uj_ui * lap_kern * dterm * points[j].mass;

						Vector3DF vr = uj_ui * spiky_kern * points[i].mass;

						vr /= points[i].density + eps;

						D += vr;
					}
					j = points[j].next_particle_index_in_the_same_cell;
				}
			}
		}

		float d = 1.414 * fabs(D.x + D.y + D.z);

		float K = 1.0f;
		float visc = 64.0f - 48.0f / (1.0f + d);
		float tau = visc * sqrtf(D.x * D.x + D.y * D.y + D.z * D.z);
		//printf("%.6f\t%.6f\t%.6f\t", d, visc, tau);
		visc = 10.0f * (tau - 0.07f) / pow(sqrtf(tau) - sqrtf(0.07f), 2);
		if (tau < 7.0f) visc = 0.0f;
		//printf("%.6f\n", visc);

		vforce *= visc;
		force += (gravity + pforce + vforce) * points[i].mass;

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
	for (int i = 0; i < points_number; i++) {
		if (points[i].particle_grid_cell_index == UNDEFINE)
			continue;

		Vector3DF acceleration = points[i].force;// + points[i].correction_pressure_force;
		acceleration /= points[i].mass;

		points[i].vel_eval += acceleration * time_step;
		points[i].pos *= sim_scale;
		points[i].pos += points[i].vel_eval * time_step;

		CollisionHandlingSimScale(points[i].pos, points[i].vel_eval);

		points[i].pos /= sim_scale;
	}
}

void ParticleSystem::ComputeDensity() {

	const float	simScaleSquare = param_[PSIMSCALE] * param_[PSIMSCALE];
	const float smoothRadiusSquare = param_[PSMOOTHRADIUS] * param_[PSMOOTHRADIUS];
	const float mass = param_[PMASS];

	for (int i = 0; i < points_number; i++) {

		points[i].density = 0.0;
		float sum = 0.0;

		const uint i_cell_index = points[i].particle_grid_cell_index;
		if (i_cell_index != UNDEFINE) {

			for (int cell = 0; cell < grid_adj_cnt; cell++) {

				const int neighbor_cell_index = i_cell_index + grid_neighbor_cell_index_offset_[cell];
				if (neighbor_cell_index == UNDEFINE || (neighbor_cell_index < 0 || neighbor_cell_index > grid_total - 1)){
					continue;
				}

				int j = grid[neighbor_cell_index].first_particle_index;
				while (j != UNDEFINE) {
					if (i == j) {
						j = points[j].next_particle_index_in_the_same_cell;
						continue;
					}
					Vector3DF ri = points[i].pos;
					Vector3DF ri_rj = ri - points[j].pos;
					const float ri_rjSquare = simScaleSquare*(ri_rj.x*ri_rj.x + ri_rj.y*ri_rj.y + ri_rj.z*ri_rj.z);
					if (ri_rjSquare <= smoothRadiusSquare) {
						points[i].density += pow(smoothRadiusSquare - ri_rjSquare, 3.0f);
					}
					j = points[j].next_particle_index_in_the_same_cell;
				}
			}
		}
	}


}

void ParticleSystem::PositionBasedFluid() {

	const float mass = param_[PMASS];
	const float simScale = param_[PSIMSCALE];
	const float simScaleSquare = simScale * simScale;
	const float smoothRadius = param_[PSMOOTHRADIUS];
	const float smoothRadiusSquare = smoothRadius * smoothRadius;
	const float denstity0 = param_[PRESTDENSITY];
	const float eps = 1.0e-6;
	float C = 0.0;
	int m_iterations = 0;
	float avg_density_devergence_err = 0;
	float simScale5 = pow(simScale, 4);

	while (((avg_density_devergence_err > 0.132) || (m_iterations < 2)) && (m_iterations < 5)) {
		ComputeDensity();
		avg_density_devergence_err = 0;
		for (int i = 0; i < points_number; i++) {
			const uint	i_cell_index = points[i].particle_grid_cell_index;
			Vector3DF	ipos = points[i].pos;

			float		sumdeta2 = 0.0;
			float		sumdevergencei = 0.0;
			Vector3DF	sumdetai(0, 0, 0);
			if (i_cell_index != UNDEFINE) {
				for (int cell = 0; cell < grid_adj_cnt; cell++) {

					const int neighbor_cell_index = i_cell_index + grid_neighbor_cell_index_offset_[cell];
					if (neighbor_cell_index == UNDEFINE || (neighbor_cell_index < 0 || neighbor_cell_index > grid_total - 1)) {
						continue;
					}

					int j = grid[neighbor_cell_index].first_particle_index;
					while (j != UNDEFINE) {
						if (i == j) {
							j = points[j].next_particle_index_in_the_same_cell;
							continue;
						}

						Vector3DF ri_rj = ipos - points[j].pos;

						const float ri_rjSquare = ri_rj.x*ri_rj.x + ri_rj.y*ri_rj.y + ri_rj.z*ri_rj.z;
						if (ri_rjSquare * simScaleSquare <= smoothRadiusSquare) {
							const float jdist = sqrt(ri_rjSquare);
							const float h_minus_r = smoothRadius/simScale - jdist;//h-|pi-pj|
							const float h_minus_r2 = h_minus_r * h_minus_r;//(h-|pi-pj|)2
							Vector3DF pi_pjnorm = ri_rj * (1.0f / (jdist + eps));//pi-pj|pi-pj|
							Vector3DF labataj = (pi_pjnorm * h_minus_r2) * lap_kern * simScale5 * mass * (1 / denstity0);
							const float ladx = labataj.x;
							const float lady = labataj.y;
							const float ladz = labataj.z;
							const float lababta2_i = (ladx*ladx + lady*lady + ladz*ladz);

							sumdeta2 += lababta2_i;

							sumdetai += labataj;
						}
						j = points[j].next_particle_index_in_the_same_cell;
					}
				}
			}

			C = (points[i].density) / denstity0 - 1;

			float sumdetai2 = (sumdetai.x*sumdetai.x + sumdetai.y*sumdetai.y + sumdetai.z*sumdetai.z);
			Vector3DF detaposi = sumdetai * ((-1)*C / (sumdeta2 + sumdetai2 + eps));
			points[i].pos += detaposi * 0.00000001f;
			avg_density_devergence_err += abs(C);
		}
		avg_density_devergence_err /= points_number;
		m_iterations++;
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

	if (gc.x > grid_res.x - 1)
		gc.x = grid_res.x - 1;
	if (gc.y > grid_res.y - 1)
		gc.y = grid_res.y - 1;
	if (gc.z > grid_res.z - 1)
		gc.z = grid_res.z - 1;

	// 将3维的cell索引转化为1维的整数索引
	return (int)(gc.y * grid_res.x  * grid_res.z + gc.z * grid_res.x + gc.x); // x*y*z + x*y + x
}

int ParticleSystem::getGridCell(int p, Vector3DI& gc){
	return getGridCell(points[p].pos, gc);
}

Vector3DI ParticleSystem::getCell(int c) {

	Vector3DI gc;
	int xz = grid_res.x*grid_res.z;
	gc.y = c / xz;
	c -= gc.y*xz;
	gc.z = c / grid_res.x;
	c -= gc.z*grid_res.x;
	gc.x = c;
	return gc;
}

void ParticleSystem::ComputeGasConstAndTimeStep() {

	float maxParticleSpeed = 4.0f;
	float courantFactor = 0.4f;

	param_[PGASCONSTANT] = 1000.0f;
	float speedOfSound = sqrt(param_[PGASCONSTANT]);
	float relevantSpeed = max(speedOfSound, maxParticleSpeed);
	time_step_sph = courantFactor * param_[PSMOOTHRADIUS] / relevantSpeed;
	time_step_pbf = 0.001;

	if (param_[PRUN_MODE] == RUN_CPU_PBF)
		time_step = time_step_pbf;
	else
		time_step = time_step_sph;
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

// 绘制函数
void ParticleSystem::Draw(Camera3D& cam, float rad) {

	glDisable(GL_LIGHTING);

	if (toggle_[PDRAWGRIDCELLS]) {
		glColor4f(1.0, 0.0, 0.0, 1.0);
		DrawGrid();
	}

	if (toggle_[PDRAWDOMAIN]) {
		DrawDomain(vec_[PBOUNDARYMIN], vec_[PBOUNDARYMAX]);
	}

	if (toggle_[PDRAWGRIDBOUND]) {
		DrawDomain(vec_[PGRIDVOLUMEMIN], vec_[PGRIDVOLUMEMAX]);
	}

	glEnable(GL_LIGHTING);

	for (int n = 0; n < this->points_number; n++) {
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
	gd /= grid_delta;

	glBegin(GL_LINES);
	for (int z = 0; z <= grid_res.z; z++) {
		for (int y = 0; y <= grid_res.y; y++) {
			gc.Set(1, y, z);
			gc /= grid_delta;
			gc += grid_min;
			glVertex3f(grid_min.x, gc.y, gc.z);
			glVertex3f(grid_max.x, gc.y, gc.z);
		}
	}
	for (int z = 0; z <= grid_res.z; z++) {
		for (int x = 0; x <= grid_res.x; x++) {
			gc.Set(x, 1, z);
			gc /= grid_delta;
			gc += grid_min;
			glVertex3f(gc.x, grid_min.y, gc.z);
			glVertex3f(gc.x, grid_max.y, gc.z);
		}
	}
	for (int y = 0; y <= grid_res.y; y++) {
		for (int x = 0; x <= grid_res.x; x++) {
			gc.Set(x, y, 1);
			gc /= grid_delta;
			gc += grid_min;
			glVertex3f(gc.x, gc.y, grid_min.z);
			glVertex3f(gc.x, gc.y, grid_max.z);
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

void ParticleSystem::setDefaultParams() {

	param_[PMAXNUM] = 8388608;
	param_[PSIMSCALE] = 0.005;
	param_[PGRID_DENSITY] = 1.0;
	param_[PVISC] = 4.0;			
	param_[PRESTDENSITY] = 1051.0; //1000.0;
	param_[PMASS] = 0.002052734; //0.001953125;
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
	param_[PFORCE_FREQ] = 8.0;
	param_[PMINLOOPPCISPH] = 3;
	param_[PMAXLOOPPCISPH] = MAX_PCISPH_LOOPS;

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

void ParticleSystem::setExampleParams() {

	switch ((int)param_[PEXAMPLE]) {
	case 0:  
		vec_[PBOUNDARYMIN].Set(-60, 0, -100);
		vec_[PBOUNDARYMAX].Set(60, 90, 0);
		vec_[PINITPARTICLEMIN].Set(0, 0, -100);
		vec_[PINITPARTICLEMAX].Set(60, 60, 0);
		param_[PFORCE_MIN] = 0.0;
		param_[PGROUND_SLOPE] = 0.0;
		break;

	case 1:
		vec_[PBOUNDARYMIN].Set(-80, 0, -80);
		vec_[PBOUNDARYMAX].Set(80, 160, 80);
		vec_[PINITPARTICLEMIN].Set(-80, 0, -80);
		vec_[PINITPARTICLEMAX].Set(80, 20, 80);
		param_[PFORCE_MIN] = 0.0;
		param_[PGROUND_SLOPE] = 0.0;
		break;
	}
}

void ParticleSystem::setKernels() {

	if (param_[PSPACINGREALWORLD] == 0) {
		param_[PSPACINGREALWORLD] = pow((float)param_[PMASS] / (float)param_[PRESTDENSITY], 1 / 3.0f);
	}
	poly6_kern = 315.0f / (64.0f * MY_PI * pow(param_[PSMOOTHRADIUS], 9));
	spiky_kern = -45.0f / (MY_PI * pow(param_[PSMOOTHRADIUS], 6));
	lap_kern = 45.0f / (MY_PI * pow(param_[PSMOOTHRADIUS], 6));
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

void ParticleSystem::setInitParticleVolume(const Vector3DI& numParticlesXYZ, const Vector3DF& lengthXYZ, const float jitter) {

	std::random_device rd;
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

	int i = points_number;
	for (int iy = 0; iy < numParticlesXYZ.y; iy++) {

		float y = 0.0 + (iy + tmpY) * param_[PSPACINGGRAPHICSWORLD];

		for (int ix = 0; ix < numParticlesXYZ.x; ix++) {

			float x = vec_[PINITPARTICLEMIN].x + (ix + tmpX) * param_[PSPACINGGRAPHICSWORLD];

			for (int iz = 0; iz < numParticlesXYZ.z; iz++) {

				float z = vec_[PINITPARTICLEMIN].z + (iz + tmpZ) * param_[PSPACINGGRAPHICSWORLD];

				if (points_number < param_[PMAXNUM]) {

					points[i].pos.Set(x + (frand() - 0.5) * jitter, y + (frand() - 0.5) * jitter, z + (frand() - 0.5) * jitter);
					points[i].vel.Set(0.0, 0.0, 0.0);
					points[i].vel_eval.Set(0.0, 0.0, 0.0);
					if (rd() % 100 < 60) {
						points[i].mass = 0.002001953;
						points[i].clr = COLORA(0.45, 0.45, 0.45, 1);
						points[i].flag = PLASMA;
					}
					else {
						points[i].mass = 0.002128906;
						points[i].clr = COLORA(0.98, 0.05, 0.05, 1);
						points[i].flag = HAEMOCYTES;
					}
					
					points_number++;
				}
				++i;
			}
		}
	}
}

void ParticleSystem::setInitParticleVolumeNew(const Vector3DI& numParticlesXYZ, const Vector3DF& lengthXYZ, const float jitter) {

	std::random_device rd;
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

	int i = points_number;
	for (int iy = 0; iy < numParticlesXYZ.y; iy++) {

		float y = vec_[PINITPARTICLEMIN].y + (iy + tmpY) * param_[PSPACINGGRAPHICSWORLD];

		for (int ix = 0; ix < numParticlesXYZ.x; ix++) {

			float x = vec_[PINITPARTICLEMIN].x + (ix + tmpX) * param_[PSPACINGGRAPHICSWORLD];

			for (int iz = 0; iz < numParticlesXYZ.z; iz++) {

				float z = vec_[PINITPARTICLEMIN].z + (iz + tmpZ) * param_[PSPACINGGRAPHICSWORLD];

				if (points_number < param_[PMAXNUM]) {

					points[i].pos.Set(x + (frand() - 0.5) * jitter, y + (frand() - 0.5) * jitter, z + (frand() - 0.5) * jitter);
					points[i].vel.Set(0.0, 0.0, 0.0);
					points[i].vel_eval.Set(0.0, 0.0, 0.0);
					if (iy > numParticlesXYZ.y * 0.3) {
						points[i].mass = 0.002001953;
						points[i].clr = COLORA(0.45, 0.45, 0.45, 1);
						points[i].flag = PLASMA;
					}
					else {
						points[i].mass = 0.002128906;
						points[i].clr = COLORA(0.98, 0.05, 0.05, 1);
						points[i].flag = HAEMOCYTES;
					}

					points_number++;
				}
				++i;
			}
		}
	}
}

void ParticleSystem::setGridAllocate(const float border) {

	float world_cellsize = param_[PGRIDSIZEREALSCALE] / param_[PSIMSCALE];

	grid_min = vec_[PGRIDVOLUMEMIN];
	grid_max = vec_[PGRIDVOLUMEMAX];

	grid_size = vec_[PGRIDVOLUMEMAX] - vec_[PGRIDVOLUMEMIN];
	grid_res.x = (int)ceil(grid_size.x / world_cellsize);
	grid_res.y = (int)ceil(grid_size.y / world_cellsize);
	grid_res.z = (int)ceil(grid_size.z / world_cellsize);

	if (grid_res.x == 0)
		grid_res.x = 1;
	if (grid_res.y == 0)
		grid_res.y = 1;
	if (grid_res.z == 0)
		grid_res.z = 1;

	grid_total = (int)(grid_res.x * grid_res.y * grid_res.z);
	if (grid_total > 10000000) {
		printf("ERROR: ParticleSystem::setGridAllocate(...): number: %d, too many cells in initGrid().\n", grid_total);
	}

	grid_size.x = grid_res.x * world_cellsize;
	grid_size.y = grid_res.y * world_cellsize;
	grid_size.z = grid_res.z * world_cellsize;
	grid_delta = grid_res;
	grid_delta /= grid_size;

	grid.reserve(grid_total + 1);

	param_[PSTAT_GMEM] = 12 * grid_total;

	int cell = 0;
	for (int y = -1; y < 2; y++)
		for (int z = -1; z < 2; z++)
			for (int x = -1; x < 2; x++)
				grid_neighbor_cell_index_offset_[cell++] = y * grid_res.x *grid_res.z + z * grid_res.x + x;
}

// 分配内存
void ParticleSystem::AllocateParticlesMemory(int cnt) {
	if (!points.empty())
		points.clear();
	points.resize(cnt);

	param_[PSTAT_PMEM] = sizeof(Particle) * cnt;
}

void ParticleSystem::Record(int param, std::string name, mint::Time& start) {

	mint::Time stop;
	stop.SetSystemTime(ACC_NSEC);
	stop = stop - start;
	param_[param] = stop.GetMSec();
	if (toggle_[PPROFILE]) printf("%s:  %s\n", name.c_str(), stop.GetReadableTime().c_str());
}

