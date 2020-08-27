//#define ENABLE_FP64

// enabling extensions
//#if defined(ENABLE_FP64)
//	#if defined(cl_khr_fp64)
//		#pragma OPENCL EXTENSION cl_khr_fp64 : enable
//	#elif defined(cl_amd_fp64)
//		#pragma OPENCL EXTENSION cl_amd_fp64 : enable
//	#else
//		#error "Double precision floating point not supported by OpenCL implementation"
//	#endif
//#endif
//
//#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
//#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable
//#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
//#pragma OPENCL EXTENSION cl_khr_local_int32_extended_atomics : enable
//
//#if defined(ENABLE_FP64)
//	#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
//	#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics : enable
//#endif

// *****************************************************************************************

// floating-point data type; depends on support from device
//#if defined(ENABLE_FP64)
//	typedef double fp_type;
//#else
//	typedef float fp_type;
//#endif

typedef float fp_type;

__constant int GAUSSIAN_2D_PARAMETERS = 6; // number of parameters
__constant int PARAM_BACKGROUND = 0;
__constant int PARAM_AMPLITUDE = 1;
__constant int PARAM_X_SPOS = 2;
__constant int PARAM_Y_SPOS = 3;
__constant int PARAM_X_SIGMA = 4;
__constant int PARAM_Y_SIGMA = 5;

__constant fp_type LAMBDA_FACTOR = 10.0f;
__constant fp_type LAMBDA_FACTOR_INV = 0.1f;

__constant fp_type MIN_EPSILON = 1e-30f;

// *****************************************************************************************

fp_type gaussian2D(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters);
fp_type gaussian2Dn(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters);

/*
//fp_type gaussian2D_derivative_0(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters);
fp_type gaussian2D_derivative_1(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters);
fp_type gaussian2D_derivative_2(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters);
fp_type gaussian2D_derivative_3(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters);
fp_type gaussian2D_derivative_4(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters);
fp_type gaussian2D_derivative_5(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters);
*/

//void gaussian2D_derivative_combined(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters, __local fp_type* const derivatives);

//void cholesky_decomposition(__global fp_type* alpha_matrix)
//void forward_substitute(__global fp_type* l_matrix, __global fp_type* temp_matrix);

// *****************************************************************************************

fp_type gaussian2D(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters)
{
	return parameters[PARAM_BACKGROUND] + parameters[PARAM_AMPLITUDE] * exp(-0.5f * ((x-parameters[PARAM_X_SPOS])*(x-parameters[PARAM_X_SPOS])/(parameters[PARAM_X_SIGMA]*parameters[PARAM_X_SIGMA]) + (y-parameters[PARAM_Y_SPOS])*(y-parameters[PARAM_Y_SPOS])/(parameters[PARAM_Y_SIGMA]*parameters[PARAM_Y_SIGMA])));
}

fp_type gaussian2Dn(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters)
{
	return (1/(2*M_PI_F*parameters[PARAM_X_SIGMA]*parameters[PARAM_Y_SIGMA])) * exp(-0.5f * ((x-parameters[PARAM_X_SPOS])*(x-parameters[PARAM_X_SPOS])/(parameters[PARAM_X_SIGMA]*parameters[PARAM_X_SIGMA]) + (y-parameters[PARAM_Y_SPOS])*(y-parameters[PARAM_Y_SPOS])/(parameters[PARAM_Y_SIGMA]*parameters[PARAM_Y_SIGMA])));
}

/*
// NOTE gaussian2D_derivative_0 always return 1.0f; optimise inline!

fp_type gaussian2D_derivative_1(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters)
{
	return exp(-0.5f * ((x-parameters[PARAM_X_SPOS])*(x-parameters[PARAM_X_SPOS])/(parameters[PARAM_X_SIGMA]*parameters[PARAM_X_SIGMA]) + (y-parameters[PARAM_Y_SPOS])*(y-parameters[PARAM_Y_SPOS])/(parameters[PARAM_Y_SIGMA]*parameters[PARAM_Y_SIGMA])));
}

fp_type gaussian2D_derivative_2(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters)
{
	return parameters[PARAM_AMPLITUDE] * exp(-0.5f * ((x-parameters[PARAM_X_SPOS])*(x-parameters[PARAM_X_SPOS])/(parameters[PARAM_X_SIGMA]*parameters[PARAM_X_SIGMA]) + (y-parameters[PARAM_Y_SPOS])*(y-parameters[PARAM_Y_SPOS])/(parameters[PARAM_Y_SIGMA]*parameters[PARAM_Y_SIGMA]))) * (x-parameters[PARAM_X_SPOS])/(parameters[PARAM_X_SIGMA]*parameters[PARAM_X_SIGMA]);
}

fp_type gaussian2D_derivative_3(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters)
{
	return parameters[PARAM_AMPLITUDE] * exp(-0.5f * ((x-parameters[PARAM_X_SPOS])*(x-parameters[PARAM_X_SPOS])/(parameters[PARAM_X_SIGMA]*parameters[PARAM_X_SIGMA]) + (y-parameters[PARAM_Y_SPOS])*(y-parameters[PARAM_Y_SPOS])/(parameters[PARAM_Y_SIGMA]*parameters[PARAM_Y_SIGMA]))) * (y-parameters[PARAM_Y_SPOS])/(parameters[PARAM_Y_SIGMA]*parameters[PARAM_Y_SIGMA]);
}

fp_type gaussian2D_derivative_4(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters)
{
	return parameters[PARAM_AMPLITUDE] * exp(-0.5f * ((x-parameters[PARAM_X_SPOS])*(x-parameters[PARAM_X_SPOS])/(parameters[PARAM_X_SIGMA]*parameters[PARAM_X_SIGMA]) + (y-parameters[PARAM_Y_SPOS])*(y-parameters[PARAM_Y_SPOS])/(parameters[PARAM_Y_SIGMA]*parameters[PARAM_Y_SIGMA]))) * (x-parameters[PARAM_X_SPOS])*(x-parameters[PARAM_X_SPOS])/(parameters[PARAM_X_SIGMA]*parameters[PARAM_X_SIGMA]*parameters[PARAM_X_SIGMA]);
}

fp_type gaussian2D_derivative_5(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters)
{
	return parameters[PARAM_AMPLITUDE] * exp(-0.5f * ((x-parameters[PARAM_X_SPOS])*(x-parameters[PARAM_X_SPOS])/(parameters[PARAM_X_SIGMA]*parameters[PARAM_X_SIGMA]) + (y-parameters[PARAM_Y_SPOS])*(y-parameters[PARAM_Y_SPOS])/(parameters[PARAM_Y_SIGMA]*parameters[PARAM_Y_SIGMA]))) * (y-parameters[PARAM_Y_SPOS])*(y-parameters[PARAM_Y_SPOS])/(parameters[PARAM_Y_SIGMA]*parameters[PARAM_Y_SIGMA]*parameters[PARAM_Y_SIGMA]);
}
*/

/*
void gaussian2D_derivative_combined(__private const fp_type x, __private const fp_type y, __global const fp_type* const parameters, __local fp_type* const derivatives)
{
	__private const fp_type xmxpos = x-parameters[PARAM_X_SPOS];
	__private const fp_type xmxpos2 = xmxpos * xmxpos;
	__private const fp_type xsig2 = parameters[PARAM_X_SIGMA] * parameters[PARAM_X_SIGMA];
	__private const fp_type xsig3 = xsig2 * parameters[PARAM_X_SIGMA];
	__private const fp_type xmxpos2dxsig2 = xmxpos2 / xsig2;
	
	__private const fp_type ymypos = y-parameters[PARAM_Y_SPOS];
	__private const fp_type ymypos2 = ymypos * ymypos;
	__private const fp_type ysig2 = parameters[PARAM_Y_SIGMA] * parameters[PARAM_Y_SIGMA];
	__private const fp_type ysig3 = ysig2 * parameters[PARAM_Y_SIGMA];
	__private const fp_type ymypos2dysig2 = ymypos2 / ysig2;
	
	__private const fp_type common_exp = exp(-0.5f * (xmxpos2dxsig2 + ymypos2dysig2));
	__private const fp_type amp_common_exp = parameters[PARAM_AMPLITUDE] * common_exp;
	
	derivatives[0] = 1.0f;
	derivatives[1] = common_exp;
	derivatives[2] = amp_common_exp * xmxpos / xsig2;
	derivatives[3] = amp_common_exp * ymypos / ysig2;
	derivatives[4] = amp_common_exp * xmxpos2 / xsig3;
	derivatives[5] = amp_common_exp * ymypos2 / ysig3;
}*/


// *****************************************************************************************

__kernel void calculate_chi2(__global const fp_type* const image_data, __global const fp_type* const x_positions, __global const fp_type* const y_positions, __private const int image_width, __private const int image_height, __global const fp_type* const parameters, __global fp_type* const chi2)
{
	// work unit
	__private const int global_id = get_global_id(0);
	__private const int image_pixel_count = image_width * image_height;
	__global const fp_type* const my_image_data = image_data+(global_id * image_pixel_count);
	__global const fp_type* const my_x_positions = x_positions+(global_id * image_pixel_count);
	__global const fp_type* const my_y_positions = y_positions+(global_id * image_pixel_count);
	__global const fp_type* const my_parameters = parameters+(global_id * GAUSSIAN_2D_PARAMETERS);
	//__global fp_type* const my_chi2 = chi2+(global_id);
	
	// loop over all pixels in image to sum the product of difference between image and model
	__private fp_type chi2_res = 0.0f;
	for(__private int y = 0; y < image_height; ++y)
	{
		for(__private int x = 0; x < image_width; ++x)
		{
			// calculate sum of squared difference
			//dy = *px - gaussian2D(x, y, my_parameters);
			__private fp_type dy = my_image_data[mad24(y,image_width,x)] - gaussian2D(my_x_positions[mad24(y,image_width,x)], my_y_positions[mad24(y,image_width,x)], my_parameters); // RSLV: inline Gaussian2D function?
			chi2_res += dy * dy;
		}
	}
	
	// store result
	chi2[global_id] = chi2_res;
}

// *****************************************************************************************

__kernel void calculate_log_mle(__global const fp_type* const image_data, __global const fp_type* const x_positions, __global const fp_type* const y_positions, __private const int image_width, __private const int image_height, __global const fp_type* const parameters, __global fp_type* const log_mle)
{
	// work unit
	__private const int global_id = get_global_id(0);
	__private const int image_pixel_count = image_width * image_height;
	__global const fp_type* const my_image_data = image_data+(global_id * image_pixel_count);
	__global const fp_type* const my_x_positions = x_positions+(global_id * image_pixel_count);
	__global const fp_type* const my_y_positions = y_positions+(global_id * image_pixel_count);
	__global const fp_type* const my_parameters = parameters+(global_id * GAUSSIAN_2D_PARAMETERS);
	//__global fp_type* const my_chi2 = chi2+(global_id);
	
	// loop over all pixels in image to sum the product of difference between image and model
	__private fp_type log_mle_res = 0.0f;
	for(__private int y = 0; y < image_height; ++y)
	{
		for(__private int x = 0; x < image_width; ++x)
		{
			// calculate log maximum likelihood
			log_mle_res += gaussian2Dn(my_x_positions[mad24(y,image_width,x)], my_y_positions[mad24(y,image_width,x)], my_parameters) * my_image_data[mad24(y,image_width,x)];
		}
	}
	
	// store result
	log_mle[global_id] = log_mle_res;
}

// *****************************************************************************************

__kernel void calculate_alpha_matrix(__global const fp_type* const x_positions, __global const fp_type* const y_positions, __private const int image_width, __private const int image_height, __global const fp_type* const parameters, __global fp_type* const alpha_matrix, __global const fp_type* const lambda)

{
	// work unit
	__private const int global_id = get_global_id(0);
	__private const int image_pixel_count = image_width * image_height;
	__global const fp_type* const my_x_positions = x_positions+(global_id * image_pixel_count);
	__global const fp_type* const my_y_positions = y_positions+(global_id * image_pixel_count);
	__global const fp_type* const my_parameters = parameters+(global_id * GAUSSIAN_2D_PARAMETERS);
	__global fp_type* my_alpha_matrix = alpha_matrix+(global_id * GAUSSIAN_2D_PARAMETERS * GAUSSIAN_2D_PARAMETERS);
	__private const fp_type my_lambda = lambda[global_id]; // NOTE: __private variable, not a __global pointer!
	
	// copy parameters to private memory
	//__private fp_type my_private_parameters[GAUSSIAN_2D_PARAMETERS];
	__private fp_type my_private_parameters[6];
	for(__private int i = 0; i < GAUSSIAN_2D_PARAMETERS; ++i)
	{
		my_private_parameters[i] = my_parameters[i];
	}
	
	// diagonal
	//__private fp_type alpha_res_00 = 0.0f; // sum of ones
	__private fp_type alpha_res_00 = image_width * image_height; // OPTIMISATION
	__private fp_type alpha_res_11 = 0.0f;
	__private fp_type alpha_res_22 = 0.0f;
	__private fp_type alpha_res_33 = 0.0f;
	__private fp_type alpha_res_44 = 0.0f;
	__private fp_type alpha_res_55 = 0.0f;
	
	// triangular matrix (mirror)
	__private fp_type alpha_res_01 = 0.0f;
	__private fp_type alpha_res_02 = 0.0f;
	__private fp_type alpha_res_03 = 0.0f;
	__private fp_type alpha_res_04 = 0.0f;
	__private fp_type alpha_res_05 = 0.0f;
	
	__private fp_type alpha_res_12 = 0.0f;
	__private fp_type alpha_res_13 = 0.0f;
	__private fp_type alpha_res_14 = 0.0f;
	__private fp_type alpha_res_15 = 0.0f;
	
	__private fp_type alpha_res_23 = 0.0f;
	__private fp_type alpha_res_24 = 0.0f;
	__private fp_type alpha_res_25 = 0.0f;
	
	__private fp_type alpha_res_34 = 0.0f;
	__private fp_type alpha_res_35 = 0.0f;
	
	__private fp_type alpha_res_45 = 0.0f;
	
	// loop over all pixel, computating all intermediate values
	for(__private int y = 0; y < image_height; ++y)
	{
		for(__private int x = 0; x < image_width; ++x)
		{
			// calculate derivates for position
			__private const fp_type xmxpos = my_x_positions[mad24(y,image_width,x)]-my_private_parameters[PARAM_X_SPOS];
			__private const fp_type xmxpos2 = xmxpos * xmxpos;
			__private const fp_type xsig2 = my_private_parameters[PARAM_X_SIGMA] * my_private_parameters[PARAM_X_SIGMA];
			__private const fp_type xsig3 = xsig2 * my_private_parameters[PARAM_X_SIGMA];
			__private const fp_type xmxpos2dxsig2 = xmxpos2 / xsig2;
			
			__private const fp_type ymypos = my_y_positions[mad24(y,image_width,x)]-my_private_parameters[PARAM_Y_SPOS];
			__private const fp_type ymypos2 = ymypos * ymypos;
			__private const fp_type ysig2 = my_private_parameters[PARAM_Y_SIGMA] * my_private_parameters[PARAM_Y_SIGMA];
			__private const fp_type ysig3 = ysig2 * my_private_parameters[PARAM_Y_SIGMA];
			__private const fp_type ymypos2dysig2 = ymypos2 / ysig2;
			
			__private const fp_type common_exp = exp(-0.5f * (xmxpos2dxsig2 + ymypos2dysig2));
			__private const fp_type amp_common_exp = my_private_parameters[PARAM_AMPLITUDE] * common_exp;
			
			//__private fp_type dy_p0 = 1.0f; // optimisation
			__private fp_type dy_p1 = common_exp;
			__private fp_type dy_p2 = amp_common_exp * xmxpos / xsig2;
			__private fp_type dy_p3 = amp_common_exp * ymypos / ysig2;
			__private fp_type dy_p4 = amp_common_exp * xmxpos2 / xsig3;
			__private fp_type dy_p5 = amp_common_exp * ymypos2 / ysig3;
			
//			// calculate derivates for position
//			//__private fp_type dy_p0 = 1.0f; // optimisation // gaussian2D_derivative_0(x, y, my_parameters);
//			__private fp_type dy_p1 = gaussian2D_derivative_1(x, y, my_parameters);
//			__private fp_type dy_p2 = gaussian2D_derivative_2(x, y, my_parameters);
//			__private fp_type dy_p3 = gaussian2D_derivative_3(x, y, my_parameters);
//			__private fp_type dy_p4 = gaussian2D_derivative_4(x, y, my_parameters);
//			__private fp_type dy_p5 = gaussian2D_derivative_5(x, y, my_parameters);
			
			// sum product
			//alpha_res_00 += dy_p0 * dy_p0; // OPTIMISATION
			alpha_res_11 += dy_p1 * dy_p1;
			alpha_res_22 += dy_p2 * dy_p2;
			alpha_res_33 += dy_p3 * dy_p3;
			alpha_res_44 += dy_p4 * dy_p4;
			alpha_res_55 += dy_p5 * dy_p5;
			
			alpha_res_01 += dy_p1; // optimisation // dy_p0 * dy_p1;
			alpha_res_02 += dy_p2; // optimisation // dy_p0 * dy_p2;
			alpha_res_03 += dy_p3; // optimisation // dy_p0 * dy_p3;
			alpha_res_04 += dy_p4; // optimisation // dy_p0 * dy_p4;
			alpha_res_05 += dy_p5; // optimisation // dy_p0 * dy_p5;
			
			alpha_res_12 += dy_p1 * dy_p2;
			alpha_res_13 += dy_p1 * dy_p3;
			alpha_res_14 += dy_p1 * dy_p4;
			alpha_res_15 += dy_p1 * dy_p5;
			
			alpha_res_23 += dy_p2 * dy_p3;
			alpha_res_24 += dy_p2 * dy_p4;
			alpha_res_25 += dy_p2 * dy_p5;
			
			alpha_res_34 += dy_p3 * dy_p4;
			alpha_res_35 += dy_p3 * dy_p5;
			
			alpha_res_45 += dy_p4 * dy_p5;
		}
	}
	
	// marquardt's lambda addition
	__private fp_type my_lambda1 = my_lambda + 1.0f;
	alpha_res_00 *= my_lambda1;
	alpha_res_11 *= my_lambda1;
	alpha_res_22 *= my_lambda1;
	alpha_res_33 *= my_lambda1;
	alpha_res_44 *= my_lambda1;
	alpha_res_55 *= my_lambda1;
	
	// store results
	*my_alpha_matrix++ = alpha_res_00;
	*my_alpha_matrix++ = alpha_res_01;
	*my_alpha_matrix++ = alpha_res_02;
	*my_alpha_matrix++ = alpha_res_03;
	*my_alpha_matrix++ = alpha_res_04;
	*my_alpha_matrix++ = alpha_res_05;
	
	*my_alpha_matrix++ = alpha_res_01;
	*my_alpha_matrix++ = alpha_res_11;
	*my_alpha_matrix++ = alpha_res_12;
	*my_alpha_matrix++ = alpha_res_13;
	*my_alpha_matrix++ = alpha_res_14;
	*my_alpha_matrix++ = alpha_res_15;
	
	*my_alpha_matrix++ = alpha_res_02;
	*my_alpha_matrix++ = alpha_res_12;
	*my_alpha_matrix++ = alpha_res_22;
	*my_alpha_matrix++ = alpha_res_23;
	*my_alpha_matrix++ = alpha_res_24;
	*my_alpha_matrix++ = alpha_res_25;
	
	*my_alpha_matrix++ = alpha_res_03;
	*my_alpha_matrix++ = alpha_res_13;
	*my_alpha_matrix++ = alpha_res_23;
	*my_alpha_matrix++ = alpha_res_33;
	*my_alpha_matrix++ = alpha_res_34;
	*my_alpha_matrix++ = alpha_res_35;
	
	*my_alpha_matrix++ = alpha_res_04;
	*my_alpha_matrix++ = alpha_res_14;
	*my_alpha_matrix++ = alpha_res_24;
	*my_alpha_matrix++ = alpha_res_34;
	*my_alpha_matrix++ = alpha_res_44;
	*my_alpha_matrix++ = alpha_res_45;
	
	*my_alpha_matrix++ = alpha_res_05;
	*my_alpha_matrix++ = alpha_res_15;
	*my_alpha_matrix++ = alpha_res_25;
	*my_alpha_matrix++ = alpha_res_35;
	*my_alpha_matrix++ = alpha_res_45;
	*my_alpha_matrix++ = alpha_res_55;
}

// *****************************************************************************************

__kernel void calculate_beta_vector(__global const fp_type* const image_data, __global const fp_type* const x_positions, __global const fp_type* const y_positions, __private const int image_width, __private const int image_height, __global fp_type* const parameters, __global fp_type* const beta_vector)
{
	// work unit
	__private const int global_id = get_global_id(0);
	__private const int image_pixel_count = image_width * image_height;
	__global const fp_type* const my_image_data = image_data+(global_id * image_pixel_count);
	__global const fp_type* const my_x_positions = x_positions+(global_id * image_pixel_count);
	__global const fp_type* const my_y_positions = y_positions+(global_id * image_pixel_count);
	__global const fp_type* const my_parameters = parameters+(global_id * GAUSSIAN_2D_PARAMETERS);
	__global fp_type* my_beta_vector = beta_vector+(global_id * GAUSSIAN_2D_PARAMETERS);
	
	// copy image data to local memory; NOTE: no point for this function, only uses each value once...
	//__local fp_type my_private_image_data[image_pixel_count];
	//for(__private int i = 0; i < image_pixel_count; ++i)
	//{
	//	my_private_image_data[i] = my_image_data[i];
	//}
	
	// copy parameters to private memory
	//__private fp_type my_private_parameters[GAUSSIAN_2D_PARAMETERS];
	__private fp_type my_private_parameters[6];
	for(__private int i = 0; i < GAUSSIAN_2D_PARAMETERS; ++i)
	{
		my_private_parameters[i] = my_parameters[i];
	}
	
	// start sums at zero
	__private fp_type my_beta_vector_0 = 0.0f;
	__private fp_type my_beta_vector_1 = 0.0f;
	__private fp_type my_beta_vector_2 = 0.0f;
	__private fp_type my_beta_vector_3 = 0.0f;
	__private fp_type my_beta_vector_4 = 0.0f;
	__private fp_type my_beta_vector_5 = 0.0f;
	
	// loop over all pixels in image
	__private fp_type px_dg = 0.0f;
	__global const fp_type* px = my_image_data;
	for(__private int y = 0; y < image_height; ++y)
	{
		for(__private int x = 0; x < image_width; ++x)
		{
			__private const fp_type xmxpos = my_x_positions[mad24(y,image_width,x)]-my_private_parameters[PARAM_X_SPOS];
			__private const fp_type xmxpos2 = xmxpos * xmxpos;
			__private const fp_type xsig2 = my_private_parameters[PARAM_X_SIGMA] * my_private_parameters[PARAM_X_SIGMA];
			__private const fp_type xsig3 = xsig2 * my_private_parameters[PARAM_X_SIGMA];
			__private const fp_type xmxpos2dxsig2 = xmxpos2 / xsig2;
			
			__private const fp_type ymypos = my_y_positions[mad24(y,image_width,x)]-my_private_parameters[PARAM_Y_SPOS];
			__private const fp_type ymypos2 = ymypos * ymypos;
			__private const fp_type ysig2 = my_private_parameters[PARAM_Y_SIGMA] * my_private_parameters[PARAM_Y_SIGMA];
			__private const fp_type ysig3 = ysig2 * my_private_parameters[PARAM_Y_SIGMA];
			__private const fp_type ymypos2dysig2 = ymypos2 / ysig2;
			
			__private const fp_type common_exp = exp(-0.5f * (xmxpos2dxsig2 + ymypos2dysig2));
			__private const fp_type amp_common_exp = my_private_parameters[PARAM_AMPLITUDE] * common_exp;
			
			//derivatives[0] = 1.0f;
			//derivatives[1] = common_exp;
			//derivatives[2] = amp_common_exp * xmxpos / xsig2;
			//derivatives[3] = amp_common_exp * ymypos / ysig2;
			//derivatives[4] = amp_common_exp * xmxpos2 / xsig3;
			//derivatives[5] = amp_common_exp * ymypos2 / ysig3;
			
			// compute intermediate shared result (difference)
			px_dg = *px - gaussian2D(my_x_positions[mad24(y,image_width,x)], my_y_positions[mad24(y,image_width,x)], my_parameters);
			my_beta_vector_0 +=  px_dg;// * gaussian2D_derivative_0(x, y, my_parameters); // NOTE: gaussian2D_derivative_0 returns 1.0f;
//			my_beta_vector_1 +=  px_dg * gaussian2D_derivative_1(x, y, my_parameters);
//			my_beta_vector_2 +=  px_dg * gaussian2D_derivative_2(x, y, my_parameters);
//			my_beta_vector_3 +=  px_dg * gaussian2D_derivative_3(x, y, my_parameters);
//			my_beta_vector_4 +=  px_dg * gaussian2D_derivative_4(x, y, my_parameters);
//			my_beta_vector_5 +=  px_dg * gaussian2D_derivative_5(x, y, my_parameters);
			my_beta_vector_1 +=  px_dg * common_exp;
			my_beta_vector_2 +=  px_dg * amp_common_exp * xmxpos / xsig2;
			my_beta_vector_3 +=  px_dg * amp_common_exp * ymypos / ysig2;
			my_beta_vector_4 +=  px_dg * amp_common_exp * xmxpos2 / xsig3;
			my_beta_vector_5 +=  px_dg * amp_common_exp * ymypos2 / ysig3;
			
			// increment pointer to pixel image data
			++px;
		}
	}
	
	// store results
	*my_beta_vector++ = my_beta_vector_0;
	*my_beta_vector++ = my_beta_vector_1;
	*my_beta_vector++ = my_beta_vector_2;
	*my_beta_vector++ = my_beta_vector_3;
	*my_beta_vector++ = my_beta_vector_4;
	*my_beta_vector++ = my_beta_vector_5;
}

// *****************************************************************************************

__kernel void cholesky_decomposition(__global fp_type* alpha_matrix)
{
	// work unit
	__private const int global_id = get_global_id(0);
	__global fp_type* const my_alpha_matrix = alpha_matrix+(global_id * GAUSSIAN_2D_PARAMETERS * GAUSSIAN_2D_PARAMETERS);
	
//	__private bool isspd = true; // is symmetric & positive-definite
	
	// loop over diagonal values
	for(__private int j = 0; j < GAUSSIAN_2D_PARAMETERS; ++j)
	{
		// STEP: calculate value of diagonal
		__private fp_type l_jj = my_alpha_matrix[mad24(j, GAUSSIAN_2D_PARAMETERS, j)];
		for(__private int k = 0; k < j; ++k)
		{
			l_jj -= my_alpha_matrix[mad24(k, GAUSSIAN_2D_PARAMETERS, j)] * my_alpha_matrix[mad24(k, GAUSSIAN_2D_PARAMETERS, j)]; /// TODO: optimize!!
		}
		
		// check if l_jj is positive; otherwise, the matrix was not positive-definite
//		if(l_jj <= 0)
//		{
//			isspd = false;
//			return false; // terminate immediately!
//		}
		
		// take square root of l_jj
		l_jj = sqrt(l_jj); // NOTE: sqrt is slow!
		
		// set new diagonal value (in-place!)
		my_alpha_matrix[mad24(j, GAUSSIAN_2D_PARAMETERS, j)] = l_jj; // TODO: optimize indexing! (with initial value index)
		
		// calculate inverse of l_jj
		__private fp_type l_jj_inv = 1.0f / l_jj;
		
		// STEP: calculate values of submatrix
		for(__private int i = j+1; i < GAUSSIAN_2D_PARAMETERS; ++i)
		{
			__private fp_type l_ij = my_alpha_matrix[mad24(j, GAUSSIAN_2D_PARAMETERS, i)];
			for(__private int k = 0; k < j; ++k)
			{
				l_ij -= my_alpha_matrix[mad24(k, GAUSSIAN_2D_PARAMETERS, i)] * my_alpha_matrix[mad24(k, GAUSSIAN_2D_PARAMETERS, j)];
			}
			l_ij *= l_jj_inv;
			my_alpha_matrix[mad24(j, GAUSSIAN_2D_PARAMETERS, i)] = l_ij;
		}
	}
	
	// reset upper triangle of matrix to zero
	for(__private int j = 1; j < GAUSSIAN_2D_PARAMETERS; ++j)
	{
		for(__private int i = 0; i < j; ++i)
		{
			my_alpha_matrix[mad24(j, GAUSSIAN_2D_PARAMETERS, i)] = 0.0f;
		}
	}
	
//	return isspd;
}

// *****************************************************************************************

__kernel void forward_substitution(__global fp_type* const l_matrix, __global fp_type* const temp_matrix)
{
	// work unit
	__private const int global_id = get_global_id(0);
	__global fp_type* const my_l_matrix = l_matrix+(global_id * GAUSSIAN_2D_PARAMETERS * GAUSSIAN_2D_PARAMETERS);
	__global fp_type* const my_temp_matrix = temp_matrix+(global_id * GAUSSIAN_2D_PARAMETERS * GAUSSIAN_2D_PARAMETERS);
	
	// initialize inverse alpha to identity matrix
	// OPTIMISED: without if(i==j) conditional statement
	for(__private int i = 0; i < GAUSSIAN_2D_PARAMETERS; ++i)
	{
		for(__private int j = i+1; j < GAUSSIAN_2D_PARAMETERS; ++j)
		{
			my_temp_matrix[mad24(i, GAUSSIAN_2D_PARAMETERS, j)] = 0.0f;
			my_temp_matrix[mad24(j, GAUSSIAN_2D_PARAMETERS, i)] = 0.0f;
		}
	}
	for(__private int d = 0; d < GAUSSIAN_2D_PARAMETERS; ++d)
	{
		my_temp_matrix[mad24(d, GAUSSIAN_2D_PARAMETERS, d)] = 1.0f;
	}
	
//	__private bool inverse_exists = true;
	
	for(__private int i = 0; i < GAUSSIAN_2D_PARAMETERS; ++i)
	{
		__private fp_type pivot_value = my_l_matrix[mad24(i, GAUSSIAN_2D_PARAMETERS, i)]; // NOTE: pivot value should not be zero; otherwise, no inverse exists!
//		if(pivot_value == 0.0f) // RSLV: set floating point margin < 1e-10?
//		{
//			inverse_exists = false;
//			return false; // terminate immediately
//		}
		
		// echelon elimination on succeeding (+1) rows
		for(__private int j = i+1; j < GAUSSIAN_2D_PARAMETERS; ++j)
		{
			// determine pivot factor
			__private fp_type pivot_factor = my_l_matrix[mad24(i, GAUSSIAN_2D_PARAMETERS, j)] / pivot_value;
			
			// subtract pivot factor multiple of pivot row from current row
			for(__private int k = 0; k < j; ++k)
			{
				my_l_matrix[mad24(k, GAUSSIAN_2D_PARAMETERS, j)] -= pivot_factor * my_l_matrix[mad24(k, GAUSSIAN_2D_PARAMETERS, i)]; // TODO: optimise!! // = 0.0f; // 
				my_temp_matrix[mad24(k, GAUSSIAN_2D_PARAMETERS, j)] -= pivot_factor * my_temp_matrix[mad24(k, GAUSSIAN_2D_PARAMETERS, i)]; // TODO: optimise!!
			}
		}
		
		// divide row by pivot value
		for(__private int k = 0; k <= i; ++k)
		{
			my_l_matrix[mad24(k, GAUSSIAN_2D_PARAMETERS, i)] /= pivot_value; // TODO: optimisation: pivot/pivot_value = 1.0f, put outside loop
			my_temp_matrix[mad24(k, GAUSSIAN_2D_PARAMETERS, i)] /= pivot_value; // TODO: optimisation: use index from previous line
//			//output[mad24(j, GAUSSIAN_2D_PARAMETERS, k)] = l_matrix[mad24(i, GAUSSIAN_2D_PARAMETERS, l)] / pivot_value; // TMP
		}
	}
	
//	return inverse_exists;
}

// *****************************************************************************************

__kernel void multiply_matrices(__global fp_type* const temp_matrix, __global fp_type* const inverse_alpha_matrix)
{
	// work unit
	__private const int global_id = get_global_id(0);
	__global fp_type* const my_temp_matrix = temp_matrix+(global_id * GAUSSIAN_2D_PARAMETERS * GAUSSIAN_2D_PARAMETERS);
	__global fp_type* const my_inverse_alpha_matrix = inverse_alpha_matrix+(global_id * GAUSSIAN_2D_PARAMETERS * GAUSSIAN_2D_PARAMETERS);
	
	// OPTIMISATION: copy temp matrix to local memory for faster computation
	// WARNING: seems to not work with the last work group of the last batch?!
//	__local fp_type my_local_temp_matrix[36];
//	for(__private int i = 0; i < GAUSSIAN_2D_PARAMETERS * GAUSSIAN_2D_PARAMETERS; ++i)
//	{
//		my_local_temp_matrix[i] = my_temp_matrix[i];
//	}
	__global fp_type* my_local_temp_matrix = my_temp_matrix;
	
	// multiply the two inverse triangular matrices (in this case a special function that transposes the alpha matrix inline for the second matrix) to get the inverse alpha matrix
	for(__private int i = 0; i < GAUSSIAN_2D_PARAMETERS; ++i)
	{
		for(__private int j = 0; j < GAUSSIAN_2D_PARAMETERS; ++j)
		//for(__private int j = 0; j <= i; ++j) // RSLV: optimisation?
		{
			__private fp_type mat_mult = 0.0f;
			for(int k = 0; k < GAUSSIAN_2D_PARAMETERS; ++k)
			{
				mat_mult += my_local_temp_matrix[mad24(j, GAUSSIAN_2D_PARAMETERS, k)] * my_local_temp_matrix[mad24(i, GAUSSIAN_2D_PARAMETERS, k)]; // NOTE in last term k,j is transposed to j,k!
			}
			my_inverse_alpha_matrix[mad24(i, GAUSSIAN_2D_PARAMETERS, j)] = mat_mult;
		}
	}
}

// *****************************************************************************************

__kernel void calculate_da_vector(__global fp_type* alpha_matrix, __global fp_type* beta_vector, __global fp_type* da_vector)
{
	// work unit
	__private const int global_id = get_global_id(0);
	__global const fp_type* const my_alpha_matrix = alpha_matrix+(global_id * GAUSSIAN_2D_PARAMETERS * GAUSSIAN_2D_PARAMETERS);
	__global const fp_type* const my_beta_vector = beta_vector+(global_id * GAUSSIAN_2D_PARAMETERS);
	__global fp_type* my_da_vector = da_vector+(global_id * GAUSSIAN_2D_PARAMETERS);
	
	// calculate da = alpha matrix * beta vector
	for(int i = 0; i < GAUSSIAN_2D_PARAMETERS; ++i)
	{
		__private fp_type da_res = 0.0f;
		for(int j = 0; j < GAUSSIAN_2D_PARAMETERS; ++j)
		{
			da_res += my_beta_vector[j] * my_alpha_matrix[mad24(i,GAUSSIAN_2D_PARAMETERS,j)];
		}
		
		// store result
		*my_da_vector = da_res;
		++my_da_vector;
	}
}

// *****************************************************************************************

__kernel void calculate_updated_parameters(__global const fp_type* const fitted_parameters, __global const fp_type* const da_vector, __global fp_type* const updated_parameters)
{
	// work unit
	__private const int param_offset = get_global_id(0) * GAUSSIAN_2D_PARAMETERS;
	__global const fp_type* my_fitted_parameters = fitted_parameters+param_offset;
	__global const fp_type* my_da_vector = da_vector+param_offset;
	__global fp_type* my_updated_parameters = updated_parameters+param_offset;
	
	// OPTIMIZED VERSION
	*my_updated_parameters++ = *my_fitted_parameters++ + *my_da_vector++; // 1
	*my_updated_parameters++ = *my_fitted_parameters++ + *my_da_vector++; // 2
	*my_updated_parameters++ = *my_fitted_parameters++ + *my_da_vector++; // 3
	*my_updated_parameters++ = *my_fitted_parameters++ + *my_da_vector++; // 4
	*my_updated_parameters++ = *my_fitted_parameters++ + *my_da_vector++; // 5
	*my_updated_parameters = *my_fitted_parameters + *my_da_vector; // 6
}

// *****************************************************************************************

// NOTE: optimisation may update current_chi2 with updated_chi2 value
__kernel void update_parameters(__global fp_type* const current_chi2, __global const fp_type* const updated_chi2, __global const fp_type* const updated_parameters, __global fp_type* const parameters, __global fp_type* const lambda)
{
	// work unit
	__private const int global_id = get_global_id(0);
	__global fp_type* const my_current_chi2 = current_chi2+(global_id);
	__global const fp_type* const my_updated_chi2 = updated_chi2+(global_id);
	__global const fp_type* my_updated_parameters = updated_parameters+(global_id * GAUSSIAN_2D_PARAMETERS);
	__global fp_type* my_parameters = parameters+(global_id * GAUSSIAN_2D_PARAMETERS);
	__global fp_type* const my_lambda = lambda+(global_id); // NOTE: __global pointer, not a __private variable!
	
	// check if new parameters are better than the old parameters
	if(*my_updated_chi2 <= *my_current_chi2) // NOTE: better model/fit means lower chi2!!
	{
		// OPTIMIZED VERSION
		*my_parameters++ = *my_updated_parameters++; // 1
		*my_parameters++ = *my_updated_parameters++; // 2
		*my_parameters++ = *my_updated_parameters++; // 3
		*my_parameters++ = *my_updated_parameters++; // 4
		*my_parameters++ = *my_updated_parameters++; // 5
		*my_parameters = *my_updated_parameters; // 6
		
		// OPTIMISATION retain updated chi^2 as current chi^2
		*my_current_chi2 = *my_updated_chi2;
		
		// decrease step size
		*my_lambda *= LAMBDA_FACTOR_INV; // divide by lambda factor
	}
	else
	{
		// increase step size
		*my_lambda *= LAMBDA_FACTOR; // multiply by lambda factor
	}
	
}

// *****************************************************************************************

__kernel void calculate_standard_error(__global fp_type* const inverse_alpha_matrix, __global const fp_type* const standard_errors)
{
	// work unit
	__private const int global_id = get_global_id(0);
	__global fp_type* const my_inverse_alpha_matrix = inverse_alpha_matrix+(global_id * GAUSSIAN_2D_PARAMETERS * GAUSSIAN_2D_PARAMETERS);
	__global fp_type* const my_standard_errors = standard_errors+(global_id * GAUSSIAN_2D_PARAMETERS);
	
	my_standard_errors[0] = sqrt(my_inverse_alpha_matrix[0]);
	my_standard_errors[1] = sqrt(my_inverse_alpha_matrix[7]);
	my_standard_errors[2] = sqrt(my_inverse_alpha_matrix[14]);
	my_standard_errors[3] = sqrt(my_inverse_alpha_matrix[21]);
	my_standard_errors[4] = sqrt(my_inverse_alpha_matrix[28]);
	my_standard_errors[5] = sqrt(my_inverse_alpha_matrix[35]);
}
