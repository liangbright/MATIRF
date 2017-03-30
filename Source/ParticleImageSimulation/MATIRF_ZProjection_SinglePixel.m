function Value=MATIRF_ZProjection_SinglePixel(x, y, z_min, z_max, ParticleFeature, MATIRF_Param, AngleIndex)

Value=integral(@(z)MATIRF_Particle_ZProfile(z, x, y, MATIRF_Param, AngleIndex, ParticleFeature), z_min, z_max);        