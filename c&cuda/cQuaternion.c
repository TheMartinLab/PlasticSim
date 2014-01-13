
/* FUNCTIONS*/

void q_rotate(quaternion *to_rotate, cVector *axis, cVector *origin, float phi, quaternion *space, quaternion *target)
{
	quaternion *q_temp, *q_axis, *q_axis_conjugate, *q_origin, *q_temp2, *ptr_space;
	float cosPhi = cos(phi*PI/360);
	float sinPhi = sin(phi*PI/360);

	ptr_space = space;

	// point to the pre-allocated temporary memory for the local variables
	q_temp = space++;
	q_axis = space++;
	q_axis_conjugate = space++;
	q_origin = space++;
	q_temp2 = space;

	space = ptr_space;

	// normalize the axis;
	v_unit(axis, axis);
	
	// turn the axis into the rotation quaternion
	v_scale(axis, sinPhi, axis);

	// convert the vector axis into a quaternion structure
	q_set(q_axis, cosPhi, axis->x, axis->y, axis->z);
	q_conjugate(q_axis, q_axis_conjugate);
	
	// convert the origin vector into a quaternion structure
	q_set(q_origin, 0, origin->x, origin->y, origin->x);

	// subtract the origin from the quaternion and store the value in temp2
	q_subtract(to_rotate, q_origin, q_temp2);

	// post multiply and store the value in temp
	q_cross(q_temp2, q_axis_conjugate, q_temp);

	// pre multiply and store the value in temp2
	q_cross(q_axis, q_temp, q_temp2);

	// add the origin back to the rotated quaternion and store the value in the target quaternion
	q_add(q_temp2, q_origin, target);
}

void q_cross(quaternion *q1, quaternion *q2, quaternion *target)
{
	target->s = q1->s   * q2->s   - 
				q1->v.x * q2->v.x - 
				q1->v.y * q2->v.y - 
				q1->v.z * q2->v.z;
	
	target->v.x =  q1->s    * q2->v.x + 
					q1->v.x * q2->s   + 
					q1->v.y * q2->v.z - 
					q1->v.z * q2->v.y;
					
	target->v.y =  q1->s    * q2->v.y - 
					q1->v.x * q2->v.z + 
					q1->v.y * q2->s   + 
					q1->v.z * q2->v.x;
	
	target->v.z =  q1->s    * q2->v.z + 
					q1->v.x * q2->v.y - 
					q1->v.y * q2->v.x + 
					q1->v.z * q2->s;
}

void q_unit(quaternion *q, float len, quaternion *target)
{
	target->s   = q->s   / len;
	target->v.x = q->v.x / len;
	target->v.y = q->v.y / len;
	target->v.z = q->v.z / len;
}

float q_len(quaternion *q)
{
	return  sqrt((q->s * q->s)+
			(q->v.x * q->v.x) + 
			(q->v.x * q->v.x) +
			(q->v.x * q->v.x));
}

void q_mult(quaternion *q, float d, quaternion *target)
{
	target->s = q->s * d;
	target->v.x = q->v.x * d;
	target->v.y = q->v.y * d;
	target->v.z = q->v.z * d;
}

void q_add(quaternion *q1, quaternion *q2, quaternion *target)
{
	target->s = q1->s + q2->s;
	v_add(&(q1->v), &(q2->v), &(target->v));
}

void q_subtract(quaternion *q1, quaternion *q2, quaternion *target)
{
	target->s = q1->s - q2->s;
	v_subtract(&(q1->v), &(q2->v), &(target->v));
}
void q_conjugate(quaternion *q, quaternion *target)
{
	target->s = q->s;
	target->v.x = q->v.x * -1;
	target->v.y = q->v.y * -1;
	target->v.z = q->v.z * -1;
}
void q_set(quaternion *q, float a, float b, float c, float d)
{
	q->s = a;
	q->v.x=b;
	q->v.y=c;
	q->v.z=d;
}
void q_scale(quaternion *q, float d, quaternion *target)
{
	target->s = q->s * d;
	v_scale(&(q->v), d, &(target->v));
}
void q_print_t(quaternion *q)
{
	printf("\n%f\t%f\t%f\t%f", q->s, q->v.x, q->v.y, q->v.z);
}