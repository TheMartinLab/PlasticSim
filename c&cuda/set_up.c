void mol_translate(cVector *translate, int num_cAtom_types, cTetrahedron *target)
{
	int i;
	for(i = 0; i < num_cAtom_types; i++)
		v_add(&(target->cAtoms[i].v), translate, &(target->cAtoms[i].v));
}
/* The q_space parameter must have enough storage capacity for 7 quaternion structures */
void mol_rotate(cVector *axis, cVector *origin, float phi, quaternion *q_space, cTetrahedron *target)
{
	int i;
	quaternion *to_rotate, *q_target;
	cVector *ptr_axis;
	to_rotate = q_space++;
	q_target = q_space++;
	
	ptr_axis = (cVector *) malloc(sizeof(cVector));
	for(i = 0; i < 5; i++)
	{
		v_set(ptr_axis, axis->x, axis->y, axis->z);
		q_set(q_target, 0, 0, 0, 0);
		to_rotate->s = 0;
		to_rotate->v.x = target->cAtoms[i].v.x;
		to_rotate->v.y = target->cAtoms[i].v.y;
		to_rotate->v.z = target->cAtoms[i].v.z;

		q_rotate(to_rotate, ptr_axis, origin, phi, q_space, q_target);

		target->cAtoms[i].v.x = q_target->v.x;
		target->cAtoms[i].v.y = q_target->v.y;
		target->cAtoms[i].v.z = q_target->v.z;
	}

	free(ptr_axis);

}
void print_lattice_xyz(char *filename, cTetrahedron *lattice, int num_mols)
{
	FILE *fp;
	int i, j;
	char *abbrev;
	int print = 1;
	cTetrahedron *ptr_lattice;
	ptr_lattice = lattice;
	fp = fopen(filename, "w");
	fprintf(fp, "%d\n", num_mols*5);
	for(i = 0; i < num_mols; i++, ptr_lattice++)
	{
		for(j = 0; j < 5; j++)
		{
			
			if(ptr_lattice->cAtoms[j].Z == 6)
			{
				abbrev = "C";
				print = 1;
			}
			else if(ptr_lattice->cAtoms[j].Z == 35)
			{
				abbrev = "Br";
				print = 1;
			}
			else
			{
				printf("\nError reading the cAtom type at line %d in set_up.c", __LINE__);
				print = 0;
			}
			if(print)
			{
				fprintf(fp, "\n%s\t%lf\t%lf\t%lf", abbrev,
						ptr_lattice->cAtoms[j].v.x, ptr_lattice->cAtoms[j].v.y,
						ptr_lattice->cAtoms[j].v.z);
				fflush(fp);
			}
		}
		//fprintf(fp, "\nnum_surr: %d", ptr_lattice->index);
		//if(ptr_lattice->num_surr != 12)
		//	printf("\nMore or less than 12 mols surrounding %d index (%d)", i, ptr_lattice->num_surr);
		//fprintf(fp, "\n");
	}
	fclose(fp);
}
 
void print_lattice_xyz_2(char *filename, cTetrahedron *lattice, int num_mols)
{
	FILE *fp;
	int i, j;
	char *abbrev;
	int print = 1;
	cTetrahedron *ptr_lattice;
	ptr_lattice = lattice;
	fp = fopen(filename, "w");
	fprintf(fp, "%d\n", num_mols*5);
	for(i = 0; i < num_mols; i++, ptr_lattice++)
	{
		for(j = 0; j < 5; j++)
		{
			if(print)
			{
				fprintf(fp, "\n%d\t%lf\t%lf\t%lf", ptr_lattice->cAtoms[j].Z,
						ptr_lattice->cAtoms[j].v.x, ptr_lattice->cAtoms[j].v.y,
						ptr_lattice->cAtoms[j].v.z);
				fflush(fp);
			}
		}
		//fprintf(fp, "\nnum_surr: %d", ptr_lattice->index);
		//if(ptr_lattice->num_surr != 12)
		//	printf("\nMore or less than 12 mols surrounding %d index (%d)", i, ptr_lattice->num_surr);
		//fprintf(fp, "\n");
	}
	fclose(fp);
}
 
void read_shells(int num_shells, int num_per_shell, cTetrahedron *shells, FILE *fp)
{
	// host vars
	cTetrahedron curr;
	cAtom c, br;
	char buf[100], *ptr_buf, *ptr_tok;
	int i, j, k;

	// init host vars
	
	// allocate memory

	// init pointers

	// do stuff
	for(i = 0; i < num_shells; i++)
	{
		//fgets(buf, sizeof(buf), fp);
		for(j = 0; j < num_per_shell; j++)
		{
			fgets(buf, sizeof(buf), fp);
			ptr_buf = buf;
			ptr_tok = strtok(ptr_buf, ",");
			c.Z = 6;
		
			ptr_tok = strtok(NULL, ",");
			c.v.x = atof(ptr_tok);
		
			ptr_tok = strtok(NULL, ",");
			c.v.y = atof(ptr_tok);
		
			ptr_tok = strtok(NULL, ",");
			c.v.z = atof(ptr_tok);		
			memcpy(&curr.cAtoms[0], &c, sizeof(cAtom));

			for(k = 1; k < 5; k++)
			{
				fgets(buf, sizeof(buf), fp);
				ptr_buf = buf;
				ptr_tok = strtok(ptr_buf, ",");
				br.Z = 35;
		
				ptr_tok = strtok(NULL, ",");
				br.v.x = atof(ptr_tok);
		
				ptr_tok = strtok(NULL, ",");
				br.v.y = atof(ptr_tok);
		
				ptr_tok = strtok(NULL, ",");
				br.v.z = atof(ptr_tok);		
				memcpy(&curr.cAtoms[k], &br, sizeof(cAtom));
			}
			memcpy(shells + num_per_shell * i + j, &curr, sizeof(cTetrahedron));
		}
		printf("%d\n", i);
	}
}
void make_first_shell_centers(cVector *centers, int num_layers, float layer_thickness)
{
	// host vars
	cVector *ptr_centers;
	int x, x_zero, y, j, k, index=0;

	// init host vars
	x_zero = 15, y = 0;	
	
	// allocate memory
	
	// init pointers
	ptr_centers = centers;
	
	// do stuff
	for(k = 0; k < num_layers; k++)
	{
		x = (x_zero + 9*k) % num_layers;
		y = 0;
		for(j = 0; j < num_layers; j++, ptr_centers++, y++)
		{
			ptr_centers->x = x;
			ptr_centers->y = y;
			ptr_centers->z = k;
			//v_print_t(ptr_centers);
			v_scale(ptr_centers, layer_thickness, ptr_centers);
			x = (x + 3) % num_layers;
			index++;
		}
	}		
	// free memory
	
	// return value
}
void put_inside(cVector corner, cTetrahedron *curr)
{
	int x = 0, y = 0, z = 0;
	cVector pos = curr->cAtoms[0].v;
	cVector translate;
	
	if(pos.x > corner.x)
		x = -1;
	else if(pos.x < 0)
		x = 1;

	if(pos.y > corner.y)
		y = -1;
	else if(pos.y < 0)
		y = 1;

	if(pos.z > corner.z)
		z = -1;
	else if(pos.z < 0)
		z = 1;

	if(x != 0 || y != 0 || z != 0)
	{
		v_set(&translate, x * corner.x, y * corner.y, z * corner.z);
		mol_translate(&translate, 5, curr);
	}
		
}
void enforce_boundary(cTetrahedron *lattice, int num_mols)
{
	float x_max = 13*8.82;
	cVector corner = {x_max, x_max, x_max};
	int i;

	for(i = 0; i < num_mols; i++)
	{
		put_inside(corner, &lattice[i]);
	}
}
void make_lattice(cTetrahedron *shells, int num_shells, int num_per_shell, cTetrahedron *lattice, 
								cVector *centers, int num_centers)
{
	// host vars
	int shell_choice;
	int stime;
	long ltime;
	int i, j;
	int mol_idx;
	cTetrahedron *ptr_shells;
	cTetrahedron *ptr_lattice;
	cVector *ptr_centers;
	
	// init host vars
	shell_choice = 0;
	ltime = time(NULL);
	stime = (unsigned) ltime/2;
	srand(stime);
	mol_idx = 0;
	
	// allocate memory
	// init pointers
	ptr_centers = centers;
	ptr_lattice = lattice;
	
	// do stuff
	for(i = 0; i < num_centers; i++, ptr_centers++)
	{
		shell_choice = rand() % num_shells;
		
		ptr_shells = shells + shell_choice * 13;
		
		memcpy(ptr_lattice, ptr_shells, sizeof(cTetrahedron) * 13);

		for(j = 0; j < 13; j++, ptr_lattice++)
		{
			mol_translate(ptr_centers, 5, ptr_lattice);
			ptr_lattice->mol_idx = mol_idx++;
		}
	}
	
	// free memory
	// return value
}