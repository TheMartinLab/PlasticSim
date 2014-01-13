// HOST FUNCTIONS
float calc_force(lj_params params, float r)
{
	float sigma6 = pow(params.r_min, 6);
		
	float sigma12 = pow(sigma6, 2);
		
	float A = 4 * params.depth * sigma12;
		
	float B = 4 * params.depth * sigma6;
		
	float r6 = pow(r, 6);
		
	float r7 = r * r6;
		
	float r12 = pow(r6, 2);
		
	float r13 = r12 * r;
		
	float potential = 12 * A / r13 - 6 * B / r7;
		
	return potential;
}
float calc_potential_1(lj_params params, float r)
{
	float sigma6 = pow(params.r_min, 6);
		
	float sigma12 = pow(sigma6, 2);
		
	float A = 4 * params.depth * sigma12;
		
	float B = 4 * params.depth * sigma6;
		
	float r6 = pow(r, 6);
		
	float r12 = pow(r6, 2);
		
	float potential = A / r12 - B / r6;
		
	return potential;

}
void test(cTetrahedron *lattice, int num_tetra, lj_params params, float delta_t, int num_steps, int walk_steps)
{
	// vars
	clock_t start, finish;
	int i;
	float time_taken, time_for_simul;
	// init vars
	// do stuff
	start = clock();

	for(i = 0; i < num_steps; i++)
	{

		//step_start = clock();
		step(lattice, num_tetra, params, delta_t, walk_steps);
		printf("\n%d", i);
		//print_lattice_xyz("steps.xyz", lattice, num_tetra);
		/*step_finish = clock();
		time_taken = (float) (step_finish - step_start);
		time_taken /= CLOCKS_PER_SEC;
		time_taken *= 1000;
		time_for_simul = num_tetra * (12*4*4) * time_taken / 1e9;
		printf("\n(%d) Time to calc the forces on %d cTetrahedron: %g", i, num_tetra, time_for_simul);*/
	}
	finish = clock();
	time_taken = (float) (finish - start);
	time_taken /= CLOCKS_PER_SEC;
	time_taken *= 1000;
	time_for_simul = num_tetra * (12*4*4) * time_taken / 1e9;
	printf("\nTotal time to calc the forces on %d cTetrahedron: %g", num_tetra, time_for_simul);
}
void *get_worst(cTetrahedron *curr, lj_params params)
{
	// vars
	int i, j, k;
	cTetrahedron *worst, *test;
	float test_ener, worst_ener, r;
	cVector pos;
	// init vars
	worst_ener = 1e6;
	// allocate memory
	// init ptrs
	// do stuff
	for(i = 0; i < curr->num_surr; i++)
	{
		test_ener = 0;
		test = curr->first_shell[i];
		// calculate the interaction potential of each molecule-molecule interaction
		for(j = 1; j < 5;  j++)
		{
			for(k = 1; k < 5; k++)
			{
				v_subtract(&curr->cAtoms[j].v, &test->cAtoms[k].v, &pos);
				r = v_length(&pos);
				test_ener += calc_potential_1(params, r);
			}
		}
		if(test_ener < worst_ener)
		{
			worst_ener = test_ener;
			worst = test;
		}
	}
	// free memory
	// return val
	return worst;
}
void resolve(cTetrahedron *curr, cTetrahedron *worst, lj_params params, float t, quaternion *q_space)
{
	// vars
	int i, j;
	cVector f_trans, f_torque, force, br_br, br_curr, br_worst, c_br,
		total_trans, total_torque, axis;
	float r, f, phi;
	// init vars
	r = 0;
	f = 0;
	// allocate memory
	// init ptrs

	// zero out certain vectors
	total_trans.x = 0;
	total_trans.y = 0;
	total_trans.z = 0;

	total_torque.x = 0;
	total_torque.y = 0;
	total_torque.z = 0;

	// do stuff
	// loop through the bromines on the first cTetrahedron
	for(i = 1; i < 5; i++)
	{
		br_curr = curr->cAtoms[i].v;
		v_subtract(&curr->cAtoms[0].v, &curr->cAtoms[i].v, &c_br);
		v_unit(&c_br, &c_br);

		for(j = 1; j < 5; j++)
		{
			br_worst = worst->cAtoms[j].v;

			v_subtract(&br_curr, &br_worst, &br_br);
			
			// calc the distance between the bromines
			r = v_length(&br_br);
			// calculate the magnitude of the force
			f = calc_force(params, r);
			// transform the magnitude into a vector
			v_unit(&br_br, &br_br);
			v_scale(&br_br, f, &force);
			
			// calc the distance moved
			v_scale(&force, 1/(4*M_BR + M_C), &f_trans);

			// calc the torque
			v_cross(&c_br, &force, &f_torque);

			v_add(&f_trans, &total_trans, &total_trans);
			v_add(&f_torque, &total_torque, &total_torque);
		}
	}
	v_scale(&total_trans, .5, &total_trans);
	// determine the amount to be moved
	v_scale(&total_trans, 1/(4*M_BR + M_C), &total_trans);
	v_scale(&total_trans, 2*t*t, &total_trans);
	v_scale(&total_torque, .5, &total_torque);
	v_unit(&total_torque, &axis);
	phi = v_length(&total_torque) / I * t * t;
	// move the first
	mol_translate(&total_trans, 5, curr);
	// rotate the first
	mol_rotate(&axis, &curr->cAtoms[0].v, phi, q_space, curr);
	// move and rotate any linked molecules
	for(i = 0; i < curr->num_linked; i++)
	{
		mol_translate(&total_trans, 5, curr->link[i]);
		mol_rotate(&axis, &curr->link[i]->cAtoms[0].v, phi, q_space, curr);
	}

	v_scale(&total_trans, -1, &total_trans);
	// move the second
	mol_translate(&total_trans, 5, curr);
	// rotate the second
	mol_rotate(&axis, &worst->cAtoms[0].v, -phi, q_space, worst);

	// move and rotate any linked molecules
	for(i = 0; i < worst->num_linked; i++)
	{
		mol_translate(&total_trans, 5, worst->link[i]);
		mol_rotate(&axis, &worst->link[i]->cAtoms[0].v, phi, q_space, curr);
	}
}
void step(cTetrahedron *lattice, int num_tetra, lj_params params, float t, int walk_steps)
{
	// vars
	int i, j;
	cTetrahedron *worst, *curr;
	quaternion *q_space;
	// init vars
	// allocate memory
	q_space = (quaternion *) malloc(sizeof(quaternion) * 7);
	// init ptrs
	// do stuff
		// loop through the lattice
	for(i = 0; i < num_tetra; i++)
	{
		curr = &lattice[i];
		for(j = 0; j < walk_steps; j++)
		{
			worst = (cTetrahedron *) get_worst(curr, params);
			resolve(curr, worst, params, t, q_space);
			if(worst->num_surr != 12)
				curr = worst->link[0];
			else
				curr = worst;
		}


	}
	// free memory
	free(q_space);
	// return val
}
void build_surroundings(cTetrahedron *lattice, int num_tetra, cTetrahedron *surface, int num_surface, float edge_length)
{
	// vars
	int i, j, too_many_surrounding, cur_num, dist_offset, max_idx, shell_idx;
	float *dist, max_val, curr_dist;
	cTetrahedron *(surr[13*26*26*20]);
	float time_taken;
	clock_t start, finish;

	// init vars

	// allocate memory
	dist = (float *) malloc(sizeof(float) * num_tetra * 20);
	memset(dist, 0.0, sizeof(float) * num_tetra * 20);
	// do stuff
	copy_surface(lattice, num_tetra, surface, num_surface, edge_length);
	find_surroundings(lattice, num_tetra, surface, num_surface, surr, dist);
	
	start = clock();
	for(i = 0; i < num_tetra; i++)
	{
		dist_offset = 20 * i;
		cur_num = 0;
		// count how many were never set
		for(j = 0; j < 20; j++)
		{
			if(dist[dist_offset + j] > 0.0)
				cur_num++;
		}

		if(cur_num > 12)
			too_many_surrounding = 1;
		else
			too_many_surrounding = 0;
		max_val = 0;
		while(too_many_surrounding)
		{
			max_idx = 0;
			max_val=0;
			for(j = 0; j < cur_num; j++)
			{
				curr_dist = dist[dist_offset+j];
				if(curr_dist > max_val)
				{
					max_val = dist[dist_offset + j];
					max_idx = j;
				}
			}

			dist[dist_offset + max_idx] = 0.0;
			cur_num--;

			if(cur_num < 13)
				too_many_surrounding = 0;
		}
		shell_idx = 0;
		lattice[i].num_surr=0;
		for(j = 0; j < 20; j++)
		{
			curr_dist = dist[dist_offset+j];
			if(curr_dist > 0.0)
			{
				lattice[i].first_shell[shell_idx] = surr[dist_offset + j];
				shell_idx++;
				lattice[i].num_surr++;
			}
		}

	}
	finish = clock();
	time_taken = (float) (finish - start);
	time_taken /= CLOCKS_PER_SEC;
	time_taken *= 1000;

	printf("\nTotal time to find the 12 closest molecules: %g ms", time_taken);
	for(i = 0; i < num_tetra; i++)
	{
		if(lattice[i].num_surr != 12)
			printf("\ncTetrahedron %d has %d surrounding", i, lattice[i].num_surr);
	}
	// free memory
	free(dist);
}
void find_surroundings(cTetrahedron *lattice, int num_tetra, 
						cTetrahedron *surface, int num_surface, 
						cTetrahedron *surr[13*26*26*20], float *dist)
{
	// vars
	int num_surr[13*26*26];
	int i, j, total, min, max, dist_idx;
	cVector pos;
	float r, avg;
	float cutoff;
	FILE *fp;
	float time_taken;
	clock_t start, finish;

	// init vars
	cutoff = 8.82;
	min = 1e6;
	total = 0;
	max = 0;
	
	fp = fopen("num_surr_output.txt", "w");
	if(fp == NULL)
		printf("\nError opening file at line %d", __LINE__);

	start = clock();
	for(j = 0; j < num_tetra; j++)
	{
		num_surr[j] = 0;
		r = 0;
		dist_idx = j * 20;
		for(i = 0; i < num_tetra; i++)
		{
			if(i == j)
				continue;
			v_subtract(&lattice[j].cAtoms[0].v, &lattice[i].cAtoms[0].v, &pos);
			r = v_length(&pos);
			if(r < cutoff)
			{
				surr[dist_idx] = &lattice[i];
				dist[dist_idx] = r;
				dist_idx++;
				num_surr[j]++;
			}
		}
		for(i = 0; i < num_surface; i++)
		{
			v_subtract(&lattice[j].cAtoms[0].v, &surface[i].cAtoms[0].v, &pos);
			r = v_length(&pos);
			if(r < cutoff)
			{
				surr[dist_idx] = &surface[i];
				dist[dist_idx] = r;
				dist_idx++;
				num_surr[j]++;
			}
		}
	}
	for(i = 0; i < num_tetra; i++)
	{
		if(num_surr[i] < min)
			min = num_surr[i];
		if(num_surr[i] > max)
			max = num_surr[i];
		total += num_surr[i];
		fprintf(fp, "\n%d\t%d", i, num_surr[i]);
	}
	finish = clock();
	time_taken = (float) (finish - start);
	time_taken /= CLOCKS_PER_SEC;
	time_taken *= 1000;

	printf("\nTotal time to find all molecules closer than %g angstroms for %d molecules.\n%g ms", cutoff, num_tetra, time_taken);
	avg = ((float) total) / num_tetra;
	printf("\nAverage number of molecules and (min/max) within %g angstroms of each molecule is: \n%g (%d/%d)", 
			cutoff, avg, min, max);

	fclose(fp);
}
void copy_surface(cTetrahedron *lattice, int num_tetra, cTetrahedron *surface, int num_surface, float edge_length)
{
	int i, j, k, idx, mol_idx=num_tetra+1;
	int surfIdx = 0, link_idx;
	cVector trans[7];
	cTetrahedron curr, *(link[8]), *center, copied;
	for(i = 0; i < num_tetra; i++)
	{
		link_idx = 0;
		if(i%10000 == 0)
			printf("\ni: %d", i);
		switch(lattice[i].loc)
		{
		//000001 0yz face
		case 1:
			v_set(&trans[0], edge_length, 0, 0);
			curr = lattice[i];
			memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
			mol_translate(&trans[0], 5, &surface[surfIdx]);
			copied = surface[surfIdx];
			surface[surfIdx].mol_idx = mol_idx++;
			link[0] = &lattice[i];
			link[1] = &surface[surfIdx];
			link_idx = 2;
			surfIdx++;
			break;
		//000010 1yz face
		case 2:
			v_set(&trans[0], -edge_length, 0, 0);
			curr = lattice[i];
			memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
			mol_translate(&trans[0], 5, &surface[surfIdx]);
			copied = surface[surfIdx];
			surface[surfIdx].mol_idx = mol_idx++;
			link[0] = &lattice[i];
			link[1] = &surface[surfIdx];
			link_idx = 2;
			surfIdx++;
			break;
		//000100 x0z face
		case 4:
			v_set(&trans[0], 0, edge_length, 0);
			curr = lattice[i];
			memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
			mol_translate(&trans[0], 5, &surface[surfIdx]);
			copied = surface[surfIdx];
			surface[surfIdx].mol_idx = mol_idx++;
			link[0] = &lattice[i];
			link[1] = &surface[surfIdx];
			link_idx = 2;
			surfIdx++;
			break;
		//000101 00z edge
		case 5:
			v_set(&trans[0], 0, edge_length, 0);
			v_set(&trans[1], edge_length, 0, 0);
			v_set(&trans[2], edge_length, edge_length, 0);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 3; j++)
			{
				curr = lattice[i];				
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		// 000110 10z edge
		case 6:
			v_set(&trans[0], 0, edge_length, 0);
			v_set(&trans[1], -edge_length, 0, 0);
			v_set(&trans[2], -edge_length, edge_length, 0);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 3; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//001000 x1z face
		case 8:
			v_set(&trans[0], 0, -edge_length, 0);
			curr = lattice[i];
			memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
			mol_translate(&trans[0], 5, &surface[surfIdx]);
			copied = surface[surfIdx];
			surface[surfIdx].mol_idx = mol_idx++;
			link[0] = &lattice[i];
			link[1] = &surface[surfIdx];
			link_idx = 2;
			surfIdx++;
			break;
		//001001 01z edge
		case 9:
			v_set(&trans[0], 0, -edge_length, 0);
			v_set(&trans[1], edge_length, 0, 0);
			v_set(&trans[2], edge_length, -edge_length, 0);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 3; j++)
			{		
				curr = lattice[i];			
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//001010 11z edge
		case 10:
			v_set(&trans[0], 0, -edge_length, 0);
			v_set(&trans[1], -edge_length, 0, 0);
			v_set(&trans[2], -edge_length, -edge_length, 0);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 3; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//010000 xy0 face
		case 16:
			v_set(&trans[0], 0, 0, edge_length);
			curr = lattice[i];
			memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
			mol_translate(&trans[0], 5, &surface[surfIdx]);
			copied = surface[surfIdx];
			surface[surfIdx].mol_idx = mol_idx++;
			lattice[i].link[0] = &surface[surfIdx];
			link[0] = &lattice[i];
			link[1] = &surface[surfIdx];
			link_idx = 2;
			surfIdx++;
			break;
		//010001 0y0 edge
		case 17:
			v_set(&trans[0], edge_length, 0, 0);
			v_set(&trans[1], 0, 0, edge_length);
			v_set(&trans[2], edge_length, 0, edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 3; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//010010 1y0 edge
		case 18:
			v_set(&trans[0], -edge_length, 0, 0);
			v_set(&trans[1], 0, 0, edge_length);
			v_set(&trans[2], -edge_length, 0, edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 3; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//010100 x00 edge
		case 20:
			v_set(&trans[0], 0, edge_length, 0);
			v_set(&trans[1], 0, 0, edge_length);
			v_set(&trans[2], 0, edge_length, edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 3; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//010101 000 corner
		case 21:
			v_set(&trans[0], edge_length, 0, 0);
			v_set(&trans[1], 0, edge_length, 0);
			v_set(&trans[2], 0, 0, edge_length);
			v_set(&trans[3], edge_length, edge_length, 0);
			v_set(&trans[4], edge_length, 0, edge_length);
			v_set(&trans[5], 0, edge_length, edge_length);
			v_set(&trans[6], edge_length, edge_length, edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 7; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//010110 100 corner
		case 22:
			v_set(&trans[0], -edge_length, 0, 0);
			v_set(&trans[1], 0, edge_length, 0);
			v_set(&trans[2], 0, 0, edge_length);
			v_set(&trans[3], -edge_length, edge_length, 0);
			v_set(&trans[4], -edge_length, 0, edge_length);
			v_set(&trans[5], 0, edge_length, edge_length);
			v_set(&trans[6], -edge_length, edge_length, edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 7; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//011000 x10 edge
		case 24:
			v_set(&trans[0], 0, -edge_length, 0);
			v_set(&trans[1], 0, 0, edge_length);
			v_set(&trans[2], 0, -edge_length, edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 3; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//011001 010 corner
		case 25:
			v_set(&trans[0], edge_length, 0, 0);
			v_set(&trans[1], 0, -edge_length, 0);
			v_set(&trans[2], 0, 0, edge_length);
			v_set(&trans[3], edge_length, -edge_length, 0);
			v_set(&trans[4], edge_length, 0, edge_length);
			v_set(&trans[5], 0, -edge_length, edge_length);
			v_set(&trans[6], edge_length, -edge_length, edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 7; j++)
			{			
				curr = lattice[i];		
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//011010 110 corner
		case 26:
			v_set(&trans[0], -edge_length, 0, 0);
			v_set(&trans[1], 0, -edge_length, 0);
			v_set(&trans[2], 0, 0, edge_length);
			v_set(&trans[3], -edge_length, -edge_length, 0);
			v_set(&trans[4], -edge_length, 0, edge_length);
			v_set(&trans[5], 0, -edge_length, edge_length);
			v_set(&trans[6], -edge_length, -edge_length, edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 7; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//100000 xy1 face
		case 32:
			v_set(&trans[0], 0, 0, -edge_length);
			curr = lattice[i];
			memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
			mol_translate(&trans[0], 5, &surface[surfIdx]);
			copied = surface[surfIdx];
			surface[surfIdx].mol_idx = mol_idx++;
			link[0] = &lattice[i];
			link[1] = &surface[surfIdx];
			link_idx = 2;
			surfIdx++;
			break;
		//100001 0y1 edge
		case 33:
			v_set(&trans[0], edge_length, 0, 0);
			v_set(&trans[1], 0, 0, -edge_length);
			v_set(&trans[2], edge_length, 0, -edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 3; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//100010 1y1 edge
		case 34:
			v_set(&trans[0], -edge_length, 0, 0);
			v_set(&trans[1], 0, 0, -edge_length);
			v_set(&trans[2], -edge_length, 0, -edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 3; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//100100 x01 edge
		case 36:
			v_set(&trans[0], 0,edge_length, 0);
			v_set(&trans[1], 0, 0, -edge_length);
			v_set(&trans[2], 0, edge_length, -edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 3; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//100101 001 corner
		case 37:
			v_set(&trans[0], edge_length, 0, 0);
			v_set(&trans[1], 0, edge_length, 0);
			v_set(&trans[2], 0, 0, -edge_length);
			v_set(&trans[3], edge_length, edge_length, 0);
			v_set(&trans[4], edge_length, 0, -edge_length);
			v_set(&trans[5], 0, edge_length, -edge_length);
			v_set(&trans[6], edge_length, edge_length, -edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 7; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//100110 101 corner
		case 38:
			v_set(&trans[0], -edge_length, 0, 0);
			v_set(&trans[1], 0, edge_length, 0);
			v_set(&trans[2], 0, 0, -edge_length);
			v_set(&trans[3], -edge_length, edge_length, 0);
			v_set(&trans[4], -edge_length, 0, -edge_length);
			v_set(&trans[5], 0, edge_length, -edge_length);
			v_set(&trans[6], -edge_length, edge_length, -edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 7; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//101000 x11 edge
		case 40:
			v_set(&trans[0], 0, -edge_length, 0);
			v_set(&trans[1], 0, 0, -edge_length);
			v_set(&trans[2], 0, -edge_length, -edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 3; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//101001 011 corner
		case 41:
			v_set(&trans[0], edge_length, 0, 0);
			v_set(&trans[1], 0, -edge_length, 0);
			v_set(&trans[2], 0, 0, -edge_length);
			v_set(&trans[3], edge_length, -edge_length, 0);
			v_set(&trans[4], edge_length, 0, -edge_length);
			v_set(&trans[5], 0, -edge_length, -edge_length);
			v_set(&trans[6], edge_length, -edge_length, -edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 7; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		//101010 111 corner
		case 42:
			v_set(&trans[0], -edge_length, 0, 0);
			v_set(&trans[1], 0, -edge_length, 0);
			v_set(&trans[2], 0, 0, -edge_length);
			v_set(&trans[3], -edge_length, -edge_length, 0);
			v_set(&trans[4], -edge_length, 0, -edge_length);
			v_set(&trans[5], 0, -edge_length, -edge_length);
			v_set(&trans[6], -edge_length, -edge_length, -edge_length);
			link[0] = &lattice[i];
			link_idx = 1;
			for(j = 0; j < 7; j++)
			{				
				curr = lattice[i];	
				memcpy(&surface[surfIdx], &curr, sizeof(cTetrahedron));
				mol_translate(&trans[j], 5, &surface[surfIdx]);
				copied = surface[surfIdx];
				surface[surfIdx].mol_idx = mol_idx++;
				link[j+1] = &surface[surfIdx];
				link_idx++;
				surfIdx++;
			}
			break;
		
		}
		if(link_idx > 8)
			printf("\nSomehow there are more than 8 molecules linked");
		// fill out the linked molecules
		if(link_idx == 0)
			lattice[i].num_linked = 0;
		for(j = 0; j < link_idx; j++)
		{
			idx = 0;
			center = link[j];
			center->num_linked = link_idx-1;
			for(k = 0; k < link_idx; k++)
			{
				if(j == k)
					continue;
				center->link[idx++] = link[k];
			}
		}
	}
	
	printf("\nFinal surfIdx: %d", surfIdx);
	
}
int set_face_flags(cTetrahedron *lattice, int num_tetra, float a_constant, float edge_length)
{
	// vars
	int i;
	cVector pos;
	cTetrahedron *ptr_lattice;
	int idx = 0;
	int num_at_faces = 0;
	int num_at_edges = 0;
	int num_at_corners = 0;
	//a_constant /= 2;
	// init vars

	// allocate memory

	// init ptrs

	ptr_lattice = lattice;

	for(i = 0; i < num_tetra; i++)
	{
		pos = ptr_lattice[i].cAtoms[0].v;
		CLEAR(ptr_lattice[i].loc);
		idx = 0;
		if(pos.x > edge_length - a_constant)
		{
			SET(ptr_lattice[i].loc, X_POS);
			idx++;
		}
		else if(pos.x < a_constant)
		{		
			SET(ptr_lattice[i].loc, X_NEG);
			idx++;
		}
		
		if(pos.y > edge_length - a_constant)
		{
			SET(ptr_lattice[i].loc, Y_POS);
			idx++;
		}
		else if(pos.y < a_constant)
		{
			SET(ptr_lattice[i].loc, Y_NEG);
			idx++;
		}
		
		if(pos.z > edge_length - a_constant)
		{
			SET(ptr_lattice[i].loc, Z_POS);
			idx++;
		}
		else if(pos.z < a_constant)
		{
			SET(ptr_lattice[i].loc, Z_NEG);
			idx++;
		}
		switch(idx)
		{
			case 1:
				num_at_faces++;
				break;
			case 2:
				num_at_edges++;
				break;
			case 3:
				num_at_corners++;
				break;
		}
	}
	printf("\nNumber of molecules on the faces, edges, and corners: %d %d %d",
			num_at_faces, num_at_edges, num_at_corners);
	return num_at_faces + num_at_edges*3 + num_at_corners*7;
}

// DEVICE FUNCTIONS
// need to make find and build surroundings 

__global__ void build_surroundings(cTetrahedron *lattice, int num_mols, cTetrahedron *surface, int num_surface); 