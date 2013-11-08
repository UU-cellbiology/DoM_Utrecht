
package DOM;

// ImageJ
//import static ij.IJ.*;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;

import ij.measure.ResultsTable;

// OpenCL
import static org.jocl.CL.*;
import org.jocl.*;

// Java
//import java.io.File;

/**
 *
 */
public class ParallelLMA
{
	/**
	 *	Configurable options
	 */
//	private static boolean DEBUG_MODE_ENABLED = false;
//	private static boolean PROFILING_MODE_ENABLED = false;
//	private static boolean AUTOMATIC_MODE_ENABLED = false;
//	private static boolean SILENT_MODE_ENABLED = false;
	
//	private static boolean USE_DOUBLE_PRECISION = false;
	private static int CL_DATATYPE_SIZE = Sizeof.cl_float;
	
//	private static int MAX_ITERATIONS = 30; // maximum number of iterations for fitting procedure
//	private static int BATCH_SIZE = 1024; // NOTE: rule of thumb, take number of compute units * max dimension for a work unit; 1024|2048 for MacBook Pro
//	private static int LWG_SIZE = 128; // 128 for MacBook Pro
	
	/**
	 *	Constants
	 */
	private static final int NUM_FITTING_PARAMETERS = 6;
	
	private static final double DEFAULT_X_SIGMA = 2.0f;
	private static final double DEFAULT_Y_SIGMA = 2.0f;
	private static final double DEFAULT_LAMBDA = 0.001f;
	
//	private static final int DEBUG_TRACE_FIT_INDEX = 0;
	
//	private static final boolean PRINT_RESULTS_IN_CSV_FORMAT = true;
//	private static final boolean USE_EXTENDED_CSV_FORMAT = true;
	
	/**
	 *	Hidden constructor
	 */
	private ParallelLMA() {}
	
	/**
	 *	Main program entry
	 */
	/*public static void main(String args[])
	{
		// check argument count
		if(args.length < 2)
		{
			System.err.println("Program requires argument to OpenCL kernel and a dataset");
			System.err.println("Usage: java RunParallelLMA [-debug] [-profiling] [-automatic] [-silent] [-single | -double] [-it(erations) #] [-lwg(-size) #] [-batch(-size) #] kernel dataset");
			System.err.println("Note: batch-size should be a multiple of lwg-size");
			System.exit(-1);
		}
		
		// parse arguments
		for(int i = 0; i < args.length - 2; ++i)
		{
			// debug mode
			if(args[i].toLowerCase().equals("-debug"))
			{
				DEBUG_MODE_ENABLED = true;
				System.err.println("Debug mode enabled!");
			}
			// profiling mode
			else if(args[i].toLowerCase().equals("-profiling"))
			{
				PROFILING_MODE_ENABLED = true;
				System.err.println("Profiling mode enabled!");
			}
			// automatic mode
			else if(args[i].toLowerCase().equals("-automatic"))
			{
				AUTOMATIC_MODE_ENABLED = true;
				System.err.println("Automatic mode enabled!");
			}
			// silent mode
			else if(args[i].toLowerCase().equals("-silent"))
			{
				SILENT_MODE_ENABLED = true;
				System.err.println("Silent mode enabled!");
			}
			// single precision
			else if(args[i].toLowerCase().equals("-single"))
			{
				USE_DOUBLE_PRECISION = false;
				CL_DATATYPE_SIZE = Sizeof.cl_float;
				System.err.println("Using single precision float point!");
			}
			// double precision
			else if(args[i].toLowerCase().equals("-double"))
			{
				USE_DOUBLE_PRECISION = true;
				CL_DATATYPE_SIZE = Sizeof.cl_double;
				System.err.println("Using double precision float point!");
			}
			// maximum number of iterations
			else if(args[i].toLowerCase().equals("-it") || args[i].toLowerCase().equals("-iterations"))
			{
				// get integer value
				++i;
				int val = MAX_ITERATIONS;
				try
				{
					val = Integer.parseInt(args[i]);
				}
				catch(NumberFormatException e)
				{
					System.err.println("Not a valid number for maximum number of iterations: " + args[i]);
					System.err.println("Maximum number of iterations should be greater than zero");
					System.exit(-1);
				}
				
				// constraints on value
				if(val < 1)
				{
					System.err.println("Not a valid number for maximum number of iterations: " + args[i]);
					System.err.println("Maximum number of iterations should be greater than zero");
					System.exit(-1);
				}
				
				// set value
				MAX_ITERATIONS = val;
				System.err.println("Maximum number of iterations set to " + MAX_ITERATIONS + "!");
			}
			// batch size
			else if(args[i].toLowerCase().equals("-batch") || args[i].toLowerCase().equals("-batch-size"))
			{
				// get integer value
				++i;
				int val = BATCH_SIZE;
				try
				{
					val = Integer.parseInt(args[i]);
				}
				catch(NumberFormatException e)
				{
					System.err.println("Not a valid number for batch size: " + args[i]);
					System.err.println("Batch size should be at least 64");
					System.err.println("Batch size should be a multiple of local workgroup size");
					System.exit(-1);
				}
				
				// constraints on value
				if(val < 64)
				{
					System.err.println("Not a valid number for batch size: " + args[i]);
					System.err.println("Batch size should be at least 64");
					System.err.println("Batch size should be a multiple of local workgroup size");
					System.exit(-1);
				}
				
				// set value
				BATCH_SIZE = val;
				System.err.println("Batch size set to " + BATCH_SIZE + "!");
			}
			// local work group size
			else if(args[i].toLowerCase().equals("-lwg") || args[i].toLowerCase().equals("-lwg-size"))
			{
				// get integer value
				++i;
				int val = BATCH_SIZE;
				try
				{
					val = Integer.parseInt(args[i]);
				}
				catch(NumberFormatException e)
				{
					System.err.println("Not a valid number for local work group size: " + args[i]);
					System.err.println("Local work group size should be at least 16");
					System.err.println("Batch size should be a multiple of local workgroup size");
					System.exit(-1);
				}
				
				// constraints on value
				if(val < 16)
				{
					System.err.println("Not a valid number for local work group size: " + args[i]);
					System.err.println("Local work group size should be at least 16");
					System.err.println("Batch size should be a multiple of local workgroup size");
					System.exit(-1);
				}
				
				// set value
				LWG_SIZE = val;
				System.err.println("Local work group size set to " + LWG_SIZE + "!");
			}
			else
			{
				System.err.println("Not a valid arguments: " + args[i]);
				System.exit(-1);
			}
		}
		
		// additional check on batch_size = k * lwg_size
		if(BATCH_SIZE % LWG_SIZE != 0)
		{
			System.err.println("Batch size should be a multiple of local workgroup size");
					System.exit(-1);
		}
		
		// extract arguments
		String opencl_kernel_filepath = args[args.length - 2];
		String dataset_filepath = args[args.length - 1];
		
		// check validity of arguments
		File opencl_kernel_fp = new File(opencl_kernel_filepath);
		if(!opencl_kernel_fp.exists())
		{
			System.err.println("Could not find OpenCL kernel file at " + opencl_kernel_filepath);
			System.exit(-1);
		}
		
		File dataset_fp = new File(dataset_filepath);
		if(!dataset_fp.exists())
		{
			System.err.println("Could not find dataset file at " + dataset_filepath);
			System.exit(-1);
		}
		
		// create new GPUBase object and load OpenCL kernel file
		GPUBase gpu = null;
		try
		{
			gpu = new GPUBase(USE_DOUBLE_PRECISION, AUTOMATIC_MODE_ENABLED, PROFILING_MODE_ENABLED, DEBUG_MODE_ENABLED);
		}
		catch(CLException cle)
		{
			System.err.println("Could not set up OpenCL environment. Maybe user cancelled device selection dialog?");
			if(DEBUG_MODE_ENABLED)
			{
				System.err.println(cle.getMessage());
			}
			System.exit(-1);
		}
		gpu.loadProgramFromFile(opencl_kernel_filepath);
		
		// open dataset file
		long start_time = System.nanoTime();
		ImagePlus dataset_image = openImage(dataset_filepath);//new ImagePlus(dataset_filepath);
		long end_time = System.nanoTime();
		
		if(PROFILING_MODE_ENABLED)
		{
			double time_diff = (double)(end_time - start_time);
			String[] units = new String[]{"nanoseconds", "microseconds", "milliseconds", "seconds", "minutes", "hours", "way too long!"};
			int unit_index = 0;
			while(time_diff > 1e3f && unit_index < 3)
			{
				time_diff /= 1e3f;
				++unit_index;
			}
			while(time_diff > 60f && unit_index < 5)
			{
				time_diff /= 60f;
				++unit_index;
			}
			
			if(PROFILING_MODE_ENABLED)
			{
				System.err.format("Loading dataset took %.1f %s\n", (time_diff), units[unit_index]);
			}
		}
		
		// DEBUG: print dimension information
		ImageStack img_stack = dataset_image.getImageStack();
		if(DEBUG_MODE_ENABLED)
		{
			int num_slices = img_stack.getSize();
			int spot_width = img_stack.getWidth();
			int spot_height = img_stack.getHeight();
			
			System.err.println("Number of slices in stack = " + num_slices);
			System.err.println("Width of images = " + spot_width);
			System.err.println("Height of images = " + spot_height);
		}
		
		// run fitting procedure
		start_time = System.nanoTime();
		run(gpu, img_stack);
		end_time = System.nanoTime();
		
		if(PROFILING_MODE_ENABLED)
		{
			double time_diff = (double)(end_time - start_time);
			String[] units = new String[]{"nanoseconds", "microseconds", "milliseconds", "seconds", "minutes", "hours", "way too long!"};
			int unit_index = 0;
			while(time_diff > 1e3f && unit_index < 3)
			{
				time_diff /= 1e3f;
				++unit_index;
			}
			while(time_diff > 60f && unit_index < 5)
			{
				time_diff /= 60f;
				++unit_index;
			}
			
			System.err.format("Complete execution took %.1f %s\n", (time_diff), units[unit_index]);
		}
		
		// close dataset file
		dataset_image.close();
		
		// clean exit
		System.err.println();
		System.err.println("Done!");
		System.exit(0);
	}*/
	
	// ************************************************************************************
	
	//private static void run(GPUBase gpu, ImageStack dataset) // old function header
	public static void run(GPUBase gpu, int batch_size, int lwg_size, int iterations, ImagePlus image, double[][] detected_particles, double[] frame_numbers, double psf_sigma, double pixel_size, boolean ignore_false_positives, ResultsTable res_table)
	{
		// get image stack
		ImageStack dataset = image.getStack();
		
		// set batch size
		//final int batch_size = BATCH_SIZE;
		//final int lwg_size = LWG_SIZE;
		final long[] lwgs = new long[]{lwg_size}; // set to null for auto
		
		// get dataset dimension
		final int num_spots = frame_numbers.length; //dataset.getSize();
		//System.err.println("num_spots = " + num_spots);
		//final int num_spots = BATCH_SIZE; // TMP: limit to 1000 images for testing purpose
		
		// get image dimensions
		final int image_width = image.getWidth();
		final int image_height = image.getHeight();
		
		// get spot dimensions
		final int spot_radius = (int)(psf_sigma * DOMConstants.FITRADIUS); // floor
		final int spot_width = 2 * spot_radius + 1; //dataset.getWidth();
		final int spot_height = 2 * spot_radius + 1; //dataset.getHeight();
		final int pixel_count = spot_width * spot_height;
		//System.err.println("spot_radius = " + spot_radius);
		//System.err.println("spot_width = " + spot_width);
		//System.err.println("spot_height = " + spot_height);

		// create buffers
		cl_mem image_data_buffer = clCreateBuffer(gpu._ocl_context, CL_MEM_READ_ONLY, batch_size * pixel_count * CL_DATATYPE_SIZE, null, null);
		cl_mem x_positions_buffer = clCreateBuffer(gpu._ocl_context, CL_MEM_READ_ONLY, batch_size * pixel_count * CL_DATATYPE_SIZE, null, null);
		cl_mem y_positions_buffer = clCreateBuffer(gpu._ocl_context, CL_MEM_READ_ONLY, batch_size * pixel_count * CL_DATATYPE_SIZE, null, null);
		cl_mem parameters_buffer = clCreateBuffer(gpu._ocl_context, CL_MEM_READ_WRITE, batch_size * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, null, null);
		cl_mem lambda_buffer = clCreateBuffer(gpu._ocl_context, CL_MEM_READ_WRITE, batch_size * CL_DATATYPE_SIZE, null, null);
		
		cl_mem current_chi2_buffer = clCreateBuffer(gpu._ocl_context, CL_MEM_READ_WRITE, batch_size * CL_DATATYPE_SIZE, null, null);
		cl_mem updated_chi2_buffer = clCreateBuffer(gpu._ocl_context, CL_MEM_READ_WRITE, batch_size * CL_DATATYPE_SIZE, null, null);
		
		cl_mem alpha_matrix_buffer = clCreateBuffer(gpu._ocl_context, CL_MEM_READ_WRITE, batch_size * NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, null, null);
		cl_mem beta_vector_buffer = clCreateBuffer(gpu._ocl_context, CL_MEM_READ_WRITE, batch_size * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, null, null);
		cl_mem temp_matrix_buffer = clCreateBuffer(gpu._ocl_context, CL_MEM_READ_WRITE, batch_size * NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, null, null);
		cl_mem inverse_alpha_matrix_buffer = clCreateBuffer(gpu._ocl_context, CL_MEM_READ_WRITE, batch_size * NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, null, null);
		cl_mem da_vector_buffer = clCreateBuffer(gpu._ocl_context, CL_MEM_READ_WRITE, batch_size * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, null, null);
		cl_mem updated_parameters_buffer = clCreateBuffer(gpu._ocl_context, CL_MEM_READ_WRITE, batch_size * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, null, null);
		
		cl_mem standard_errors_buffer = clCreateBuffer(gpu._ocl_context, CL_MEM_WRITE_ONLY, batch_size * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, null, null);
		
		// create kernels
		cl_kernel calculate_current_chi2 = clCreateKernel(gpu._ocl_program, "calculate_chi2", null);
		clSetKernelArg(calculate_current_chi2, 0, Sizeof.cl_mem, Pointer.to(image_data_buffer));
		clSetKernelArg(calculate_current_chi2, 1, Sizeof.cl_mem, Pointer.to(x_positions_buffer));
		clSetKernelArg(calculate_current_chi2, 2, Sizeof.cl_mem, Pointer.to(y_positions_buffer));
		clSetKernelArg(calculate_current_chi2, 3, Sizeof.cl_int, Pointer.to(new int[]{spot_width}));
		clSetKernelArg(calculate_current_chi2, 4, Sizeof.cl_int, Pointer.to(new int[]{spot_height}));
		clSetKernelArg(calculate_current_chi2, 5, Sizeof.cl_mem, Pointer.to(parameters_buffer));
		clSetKernelArg(calculate_current_chi2, 6, Sizeof.cl_mem, Pointer.to(current_chi2_buffer));
		
		cl_kernel calculate_alpha_matrix = clCreateKernel(gpu._ocl_program, "calculate_alpha_matrix", null);
		clSetKernelArg(calculate_alpha_matrix, 0, Sizeof.cl_mem, Pointer.to(x_positions_buffer));
		clSetKernelArg(calculate_alpha_matrix, 1, Sizeof.cl_mem, Pointer.to(y_positions_buffer));
		clSetKernelArg(calculate_alpha_matrix, 2, Sizeof.cl_int, Pointer.to(new int[]{spot_width}));
		clSetKernelArg(calculate_alpha_matrix, 3, Sizeof.cl_int, Pointer.to(new int[]{spot_height}));
		clSetKernelArg(calculate_alpha_matrix, 4, Sizeof.cl_mem, Pointer.to(parameters_buffer));
		clSetKernelArg(calculate_alpha_matrix, 5, Sizeof.cl_mem, Pointer.to(alpha_matrix_buffer));
		clSetKernelArg(calculate_alpha_matrix, 6, Sizeof.cl_mem, Pointer.to(lambda_buffer));
		
		cl_kernel calculate_beta_vector = clCreateKernel(gpu._ocl_program, "calculate_beta_vector", null);
		clSetKernelArg(calculate_beta_vector, 0, Sizeof.cl_mem, Pointer.to(image_data_buffer));
		clSetKernelArg(calculate_beta_vector, 1, Sizeof.cl_mem, Pointer.to(x_positions_buffer));
		clSetKernelArg(calculate_beta_vector, 2, Sizeof.cl_mem, Pointer.to(y_positions_buffer));
		clSetKernelArg(calculate_beta_vector, 3, Sizeof.cl_int, Pointer.to(new int[]{spot_width}));
		clSetKernelArg(calculate_beta_vector, 4, Sizeof.cl_int, Pointer.to(new int[]{spot_height}));
		clSetKernelArg(calculate_beta_vector, 5, Sizeof.cl_mem, Pointer.to(parameters_buffer));
		clSetKernelArg(calculate_beta_vector, 6, Sizeof.cl_mem, Pointer.to(beta_vector_buffer));
		
		cl_kernel cholesky_decomposition = clCreateKernel(gpu._ocl_program, "cholesky_decomposition", null);
		clSetKernelArg(cholesky_decomposition, 0, Sizeof.cl_mem, Pointer.to(alpha_matrix_buffer));
		
		cl_kernel forward_substitution = clCreateKernel(gpu._ocl_program, "forward_substitution", null);
		clSetKernelArg(forward_substitution, 0, Sizeof.cl_mem, Pointer.to(alpha_matrix_buffer));
		clSetKernelArg(forward_substitution, 1, Sizeof.cl_mem, Pointer.to(temp_matrix_buffer));
		
		cl_kernel multiply_matrices = clCreateKernel(gpu._ocl_program, "multiply_matrices", null);
		clSetKernelArg(multiply_matrices, 0, Sizeof.cl_mem, Pointer.to(temp_matrix_buffer));
		clSetKernelArg(multiply_matrices, 1, Sizeof.cl_mem, Pointer.to(inverse_alpha_matrix_buffer));
		
		cl_kernel calculate_da_vector = clCreateKernel(gpu._ocl_program, "calculate_da_vector", null);
		clSetKernelArg(calculate_da_vector, 0, Sizeof.cl_mem, Pointer.to(inverse_alpha_matrix_buffer));
		clSetKernelArg(calculate_da_vector, 1, Sizeof.cl_mem, Pointer.to(beta_vector_buffer));
		clSetKernelArg(calculate_da_vector, 2, Sizeof.cl_mem, Pointer.to(da_vector_buffer));
		
		cl_kernel calculate_updated_parameters = clCreateKernel(gpu._ocl_program, "calculate_updated_parameters", null);
		clSetKernelArg(calculate_updated_parameters, 0, Sizeof.cl_mem, Pointer.to(parameters_buffer));
		clSetKernelArg(calculate_updated_parameters, 1, Sizeof.cl_mem, Pointer.to(da_vector_buffer));
		clSetKernelArg(calculate_updated_parameters, 2, Sizeof.cl_mem, Pointer.to(updated_parameters_buffer));
		 
		cl_kernel calculate_updated_chi2 = clCreateKernel(gpu._ocl_program, "calculate_chi2", null);
		clSetKernelArg(calculate_updated_chi2, 0, Sizeof.cl_mem, Pointer.to(image_data_buffer));
		clSetKernelArg(calculate_updated_chi2, 1, Sizeof.cl_mem, Pointer.to(x_positions_buffer));
		clSetKernelArg(calculate_updated_chi2, 2, Sizeof.cl_mem, Pointer.to(y_positions_buffer));
		clSetKernelArg(calculate_updated_chi2, 3, Sizeof.cl_int, Pointer.to(new int[]{spot_width}));
		clSetKernelArg(calculate_updated_chi2, 4, Sizeof.cl_int, Pointer.to(new int[]{spot_height}));
		clSetKernelArg(calculate_updated_chi2, 5, Sizeof.cl_mem, Pointer.to(updated_parameters_buffer));
		clSetKernelArg(calculate_updated_chi2, 6, Sizeof.cl_mem, Pointer.to(updated_chi2_buffer));
		
		cl_kernel update_parameters = clCreateKernel(gpu._ocl_program, "update_parameters", null);
		clSetKernelArg(update_parameters, 0, Sizeof.cl_mem, Pointer.to(current_chi2_buffer));
		clSetKernelArg(update_parameters, 1, Sizeof.cl_mem, Pointer.to(updated_chi2_buffer));
		clSetKernelArg(update_parameters, 2, Sizeof.cl_mem, Pointer.to(updated_parameters_buffer));
		clSetKernelArg(update_parameters, 3, Sizeof.cl_mem, Pointer.to(parameters_buffer));
		clSetKernelArg(update_parameters, 4, Sizeof.cl_mem, Pointer.to(lambda_buffer));
		
		cl_kernel calculate_standard_error = clCreateKernel(gpu._ocl_program, "calculate_standard_error", null);
		clSetKernelArg(calculate_standard_error, 0, Sizeof.cl_mem, Pointer.to(inverse_alpha_matrix_buffer));
		clSetKernelArg(calculate_standard_error, 1, Sizeof.cl_mem, Pointer.to(standard_errors_buffer));
		
		// ******************************************************************************
		
//		// profiling
//		long start_time = 0l;
//		long end_time = 0l;
//		cl_event kernel_event = new cl_event();
//		long[] profiling_result = new long[1];
//		
//		long transform_data_profiling = 0l;
//		long write_image_data_buffer_profiling = 0l;
//		long write_parameters_buffer_profiling = 0l;
//		long write_lambda_buffer_profiling = 0l;
//		long calculate_current_chi2_profiling = 0l;
//		long calculate_alpha_matrix_profiling = 0l;
//		long calculate_beta_vector_profiling = 0l;
//		long cholesky_decomposition_profiling = 0l;
//		long forward_substitution_profiling = 0l;
//		long multiply_matrices_profiling = 0l;
//		long calculate_da_vector_profiling = 0l;
//		long calculate_updated_parameters_profiling = 0l;
//		long calculate_updated_chi2_profiling = 0l;
//		long update_parameters_profiling = 0l;
//		long read_parameters_buffer_profiling = 0l;
		
		// loop over all batches
		long start_time = System.nanoTime();
		for(int ii = 0; ii < num_spots; ii+=batch_size) // proces images in batches
		{
			// show progress bar
			IJ.showProgress(ii, num_spots);
			
			// determine actual batch size
			//int actual_batch_size = Math.min(batch_size, num_spots - ii);
			//System.err.println("Actual batch size = " + actual_batch_size);
			
//			// start profiling of data conversion
//			start_time = System.nanoTime();
			
			// fill buffers with image data and initial parameters
//			double[] image_data_d = new double[batch_size * pixel_count];
//			double[] parameters_d = new double[batch_size * NUM_FITTING_PARAMETERS];
//			double[] lambda_d = new double[batch_size];
			float[] image_data_f = new float[batch_size * pixel_count];
			float[] x_positions_f = new float[batch_size * pixel_count];
			float[] y_positions_f = new float[batch_size * pixel_count];
			float[] parameters_f = new float[batch_size * NUM_FITTING_PARAMETERS];
			float[] lambda_f = new float[batch_size];
			
			int data_index = 0;
			int parameters_index = 0;
			int lambda_index = 0;
			//for(int slice = ii+1; slice <= ii+batch_size && slice <= num_spots; ++slice) // NOTE: +1 offset for ImageJ && check if less than num_spots!!
			for(int spot = ii; spot < ii+batch_size && spot < num_spots; ++spot)
			{
				// get processor for slice where spot is present
				ImageProcessor img_proc = dataset.getProcessor((int)frame_numbers[spot]);
				
				// convert image to data buffer
				double min = Double.MAX_VALUE;
				double max = Double.MIN_VALUE;
				int sx = (int)Math.round(detected_particles[0][spot]) - spot_radius;
				int ex = (int)Math.round(detected_particles[0][spot]) + spot_radius;
				int sy = (int)Math.round(detected_particles[1][spot]) - spot_radius;
				int ey = (int)Math.round(detected_particles[1][spot]) + spot_radius;
				for(int y = sy; y <= ey; ++y)
				{
					for(int x = sx; x <= ex; ++x)
					{
						// get pixel value from image
						double px = (double)(img_proc.get(x,y));
						
						// set pixel value in data input buffer
//						if(USE_DOUBLE_PRECISION)
//						{
//							image_data_d[image_data_index] = px;
//							++image_data_index;
//						}
//						else
//						{
							image_data_f[data_index] = (float)px;
							x_positions_f[data_index] = (float)x;
							y_positions_f[data_index] = (float)y;
							++data_index;
//						}
						
						// check min/max
						if(px > max) max = px;
						if(px < min) min = px;
					}
				}
				
//				if(USE_DOUBLE_PRECISION)
//				{
//					// set initial fitting parameters
//					parameters_d[parameters_index+0] = min;
//					parameters_d[parameters_index+1] = (max - min); // amplitude = max - min
//					parameters_d[parameters_index+2] = detected_particles[0][spot]; //(spot_width >> 1);// - 1 + Math.random()*2.0;
//					parameters_d[parameters_index+3] = detected_particles[1][spot]; //(spot_height >> 1);// - 1 + Math.random()*2.0;
//					parameters_d[parameters_index+4] = psf_sigma; //DEFAULT_X_SIGMA;
//					parameters_d[parameters_index+5] = psf_sigma; //DEFAULT_Y_SIGMA;
//					parameters_index += NUM_FITTING_PARAMETERS;
//					
//					// set initial lambda parameter
//					lambda_d[lambda_index] = DEFAULT_LAMBDA;
//					++lambda_index;
//				}
//				else
//				{
					// set initial fitting parameters
					parameters_f[parameters_index+0] = (float)min;
					parameters_f[parameters_index+1] = (float)(max - min); // amplitude = max - min
					parameters_f[parameters_index+2] = (float)detected_particles[0][spot]; //(float)(spot_width >> 1);// - 1 + Math.random()*2.0;
					parameters_f[parameters_index+3] = (float)detected_particles[1][spot]; //(float)(spot_height >> 1);// - 1 + Math.random()*2.0;
					parameters_f[parameters_index+4] = (float)psf_sigma; //(float)DEFAULT_X_SIGMA;
					parameters_f[parameters_index+5] = (float)psf_sigma; //(float)DEFAULT_Y_SIGMA;
					parameters_index += NUM_FITTING_PARAMETERS;
					
					// set initial lambda parameter
					lambda_f[lambda_index] = (float)DEFAULT_LAMBDA;
					++lambda_index;
//				}
			}
			
//			end_time = System.nanoTime();
//			if(PROFILING_MODE_ENABLED)
//			{
//				transform_data_profiling += end_time - start_time;
//			}
			
			// ******************************************************************************
			
			// write data to buffers
			clEnqueueWriteBuffer(gpu._ocl_queue, image_data_buffer, true, 0, batch_size * pixel_count * CL_DATATYPE_SIZE, Pointer.to(image_data_f), 0, null, null);
			
			clEnqueueWriteBuffer(gpu._ocl_queue, x_positions_buffer, true, 0, batch_size * pixel_count * CL_DATATYPE_SIZE, Pointer.to(x_positions_f), 0, null, null);
			
			clEnqueueWriteBuffer(gpu._ocl_queue, y_positions_buffer, true, 0, batch_size * pixel_count * CL_DATATYPE_SIZE, Pointer.to(y_positions_f), 0, null, null);
			
//			if(PROFILING_MODE_ENABLED)
//			{
//				clWaitForEvents(1, new cl_event[]{kernel_event});
//				clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, Sizeof.cl_long, Pointer.to(profiling_result), null);
//				start_time = profiling_result[0];
//				clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END, Sizeof.cl_long, Pointer.to(profiling_result), null);
//				end_time = profiling_result[0];
//				write_image_data_buffer_profiling += end_time - start_time;
//			}
			
			clEnqueueWriteBuffer(gpu._ocl_queue, parameters_buffer, true, 0, batch_size * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, Pointer.to(parameters_f), 0, null, null);
			
//			if(PROFILING_MODE_ENABLED)
//			{
//				clWaitForEvents(1, new cl_event[]{kernel_event});
//				clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, Sizeof.cl_long, Pointer.to(profiling_result), null);
//				start_time = profiling_result[0];
//				clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END, Sizeof.cl_long, Pointer.to(profiling_result), null);
//				end_time = profiling_result[0];
//				write_parameters_buffer_profiling += end_time - start_time;
//			}
			
			clEnqueueWriteBuffer(gpu._ocl_queue, lambda_buffer, true, 0, batch_size * CL_DATATYPE_SIZE, Pointer.to(lambda_f), 0, null, null);
			
//			if(PROFILING_MODE_ENABLED)
//			{
//				clWaitForEvents(1, new cl_event[]{kernel_event});
//				clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, Sizeof.cl_long, Pointer.to(profiling_result), null);
//				start_time = profiling_result[0];
//				clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END, Sizeof.cl_long, Pointer.to(profiling_result), null);
//				end_time = profiling_result[0];
//				write_lambda_buffer_profiling += end_time - start_time;
//			}
			
			// ******************************************************************************
			
			// STEP: calculate chi2 for current parameters
			// NOTE: updated chi2 is retained in update_parameters
			// NOTE: updated chi2 is current chi2 in next iteration!
			clEnqueueNDRangeKernel(gpu._ocl_queue, calculate_current_chi2, 1, null, new long[]{batch_size}, lwgs, 0, null, null);
			
//			if(PROFILING_MODE_ENABLED)
//			{
//				clWaitForEvents(1, new cl_event[]{kernel_event});
//				clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, Sizeof.cl_long, Pointer.to(profiling_result), null);
//				start_time = profiling_result[0];
//				clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END, Sizeof.cl_long, Pointer.to(profiling_result), null);
//				end_time = profiling_result[0];
//				calculate_current_chi2_profiling += end_time - start_time;
//			}
			
			// iterate fixed number of times
			int iteration_count = 0;
			while(iteration_count < iterations)
			{
//				// DEBUG: print intermediate results
//				if(DEBUG_MODE_ENABLED)
//				{
//					float[] lambda_dbp = new float[batch_size];
//					clEnqueueReadBuffer(gpu._ocl_queue, lambda_buffer, CL_TRUE, 0, batch_size * CL_DATATYPE_SIZE, Pointer.to(lambda_dbp), 0, null, null);
//					
//					System.err.println("\n\n****************************************");
//					System.err.println("iteration count = " + iteration_count);
//					System.err.println("lambda = " + lambda_dbp[DEBUG_TRACE_FIT_INDEX]);
//					System.err.println("lambda factor = " + 10.0f);
//				}
				
//				// DEBUG: print intermediate results
//				if(DEBUG_MODE_ENABLED)
//				{
//					float[] params_dbp = new float[batch_size * NUM_FITTING_PARAMETERS];
//					clEnqueueReadBuffer(gpu._ocl_queue, parameters_buffer, CL_TRUE, 0, batch_size * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, Pointer.to(params_dbp), 0, null, null);
//		
//					System.err.format("current parameters = [%f, %f, %f, %f, %f, %f]\n", params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+0], params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+1], params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+2], params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+3], params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+4], params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+5]);
//				}
				
				// STEP: calculate current chi2
				// NOTE: updated chi2 is retained in update_parameters
				// NOTE: updated chi2 is current chi2 in next iteration!
				
//				// DEBUG: print intermediate results
//				float[] current_chi2_dbp = new float[batch_size]; // NOTE: outside if statement on purpose
//				//clEnqueueReadBuffer(gpu._ocl_queue, current_chi2_buffer, CL_TRUE, 0, batch_size * CL_DATATYPE_SIZE, Pointer.to(current_chi2_dbp), 0, null, null); // TMP
//				if(DEBUG_MODE_ENABLED)
//				{
//					//float[] current_chi2_dbp = new float[batch_size];
//					clEnqueueReadBuffer(gpu._ocl_queue, current_chi2_buffer, CL_TRUE, 0, batch_size * CL_DATATYPE_SIZE, Pointer.to(current_chi2_dbp), 0, null, null);
//					System.err.println("current chi2 = " + current_chi2_dbp[DEBUG_TRACE_FIT_INDEX]);
//				}
				
				// STEP: calculate alpha matrix
				clEnqueueNDRangeKernel(gpu._ocl_queue, calculate_alpha_matrix, 1, null, new long[]{batch_size}, lwgs, 0, null, null);
				
//				// DEBUG: print intermediate results
//				if(DEBUG_MODE_ENABLED)
//				{
//					float[] alpha_matrix_dbp = new float[batch_size * NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS];
//					clEnqueueReadBuffer(gpu._ocl_queue, alpha_matrix_buffer, CL_TRUE, 0, batch_size * NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, Pointer.to(alpha_matrix_dbp), 0, null, null);
//					System.err.println("alpha matrix = ");
//					for(int dbp_i = 0; dbp_i < NUM_FITTING_PARAMETERS; ++dbp_i)
//					{
//						for(int dbp_j = 0; dbp_j < NUM_FITTING_PARAMETERS; ++dbp_j)
//						{
//							System.err.print(alpha_matrix_dbp[NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS * DEBUG_TRACE_FIT_INDEX + NUM_FITTING_PARAMETERS * dbp_j + dbp_i] + " "); 
//						}
//						System.err.println();
//					}
//				}
//				
//				if(PROFILING_MODE_ENABLED)
//				{
//					clWaitForEvents(1, new cl_event[]{kernel_event});
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					start_time = profiling_result[0];
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					end_time = profiling_result[0];
//					calculate_alpha_matrix_profiling += end_time - start_time;
//				}
				
				// STEP: calculate beta vector
				clEnqueueNDRangeKernel(gpu._ocl_queue, calculate_beta_vector, 1, null, new long[]{batch_size}, lwgs, 0, null, null);
				
//				// DEBUG: print intermediate results
//				if(DEBUG_MODE_ENABLED)
//				{
//					float[] beta_vector_dbp = new float[batch_size * NUM_FITTING_PARAMETERS];
//					clEnqueueReadBuffer(gpu._ocl_queue, beta_vector_buffer, CL_TRUE, 0, batch_size * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, Pointer.to(beta_vector_dbp), 0, null, null);
//		
//					System.err.format("beta vector = [%f, %f, %f, %f, %f, %f]\n", beta_vector_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+0], beta_vector_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+1], beta_vector_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+2], beta_vector_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+3], beta_vector_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+4], beta_vector_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+5]);
//				}
//				
//				if(PROFILING_MODE_ENABLED)
//				{
//					clWaitForEvents(1, new cl_event[]{kernel_event});
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					start_time = profiling_result[0];
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					end_time = profiling_result[0];
//					calculate_beta_vector_profiling += end_time - start_time;
//				}
				
				// STEP: (a) invert alpha matrix through cholesky decomposition
				clEnqueueNDRangeKernel(gpu._ocl_queue, cholesky_decomposition, 1, null, new long[]{batch_size}, lwgs, 0, null, null);
				
//				// DEBUG: print intermediate results
//				if(DEBUG_MODE_ENABLED)
//				{
//					float[] cholesky_decomposition_dbp = new float[batch_size * NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS];
//					clEnqueueReadBuffer(gpu._ocl_queue, alpha_matrix_buffer, CL_TRUE, 0, batch_size * NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, Pointer.to(cholesky_decomposition_dbp), 0, null, null);
//					System.err.println("cholesky decomposition = ");
//					for(int dbp_i = 0; dbp_i < NUM_FITTING_PARAMETERS; ++dbp_i)
//					{
//						for(int dbp_j = 0; dbp_j < NUM_FITTING_PARAMETERS; ++dbp_j)
//						{
//							System.err.print(cholesky_decomposition_dbp[NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS * DEBUG_TRACE_FIT_INDEX + NUM_FITTING_PARAMETERS * dbp_j + dbp_i] + " "); 
//						}
//						System.err.println();
//					}
//				}
//				
//				if(PROFILING_MODE_ENABLED)
//				{
//					clWaitForEvents(1, new cl_event[]{kernel_event});
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					start_time = profiling_result[0];
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					end_time = profiling_result[0];
//					cholesky_decomposition_profiling += end_time - start_time;
//				}
				
				// STEP: (b) invert alpha matrix through cholesky decomposition
				clEnqueueNDRangeKernel(gpu._ocl_queue, forward_substitution, 1, null, new long[]{batch_size}, lwgs, 0, null, null);
				
//				// DEBUG: print intermediate results
//				if(DEBUG_MODE_ENABLED)
//				{
//					float[] identity_matrix_dbp = new float[batch_size * NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS];
//					float[] forward_substitution_dbp = new float[batch_size * NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS];
//					clEnqueueReadBuffer(gpu._ocl_queue, alpha_matrix_buffer, CL_TRUE, 0, batch_size * NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, Pointer.to(identity_matrix_dbp), 0, null, null);
//					clEnqueueReadBuffer(gpu._ocl_queue, temp_matrix_buffer, CL_TRUE, 0, batch_size * NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, Pointer.to(forward_substitution_dbp), 0, null, null);
//					System.err.println("forward substitution = ");
//					for(int dbp_i = 0; dbp_i < NUM_FITTING_PARAMETERS; ++dbp_i)
//					{
//						for(int dbp_j = 0; dbp_j < NUM_FITTING_PARAMETERS; ++dbp_j)
//						{
//							System.err.print(identity_matrix_dbp[NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS * DEBUG_TRACE_FIT_INDEX + NUM_FITTING_PARAMETERS * dbp_j + dbp_i] + " "); 
//						}
//						System.err.print("| ");
//						for(int dbp_j = 0; dbp_j < NUM_FITTING_PARAMETERS; ++dbp_j)
//						{
//							System.err.print(forward_substitution_dbp[NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS * DEBUG_TRACE_FIT_INDEX + NUM_FITTING_PARAMETERS * dbp_j + dbp_i] + " "); 
//						}
//						System.err.println();
//					}
//				}
//				
//				if(PROFILING_MODE_ENABLED)
//				{
//					clWaitForEvents(1, new cl_event[]{kernel_event});
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					start_time = profiling_result[0];
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					end_time = profiling_result[0];
//					forward_substitution_profiling += end_time - start_time;
//				}
				
				// STEP: (c) invert alpha matrix through cholesky decomposition
				clEnqueueNDRangeKernel(gpu._ocl_queue, multiply_matrices, 1, null, new long[]{batch_size}, lwgs, 0, null, null); // RSLV: single inver_alpha_matrix, or separate calls to choleksy_decomposition and forward_substitute?
				
//				// DEBUG: print intermediate results
//				if(DEBUG_MODE_ENABLED)
//				{
//					float[] inverse_alpha_matrix_dbp = new float[batch_size * NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS];
//					clEnqueueReadBuffer(gpu._ocl_queue, inverse_alpha_matrix_buffer, CL_TRUE, 0, batch_size * NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, Pointer.to(inverse_alpha_matrix_dbp), 0, null, null);
//					System.err.println("inverse alpha matrix = ");
//					for(int dbp_i = 0; dbp_i < NUM_FITTING_PARAMETERS; ++dbp_i)
//					{
//						for(int dbp_j = 0; dbp_j < NUM_FITTING_PARAMETERS; ++dbp_j)
//						{
//							System.err.print(inverse_alpha_matrix_dbp[NUM_FITTING_PARAMETERS * NUM_FITTING_PARAMETERS * DEBUG_TRACE_FIT_INDEX + NUM_FITTING_PARAMETERS * dbp_j + dbp_i] + " "); 
//						}
//						System.err.println();
//					}
//				}
//				
//				if(PROFILING_MODE_ENABLED)
//				{
//					clWaitForEvents(1, new cl_event[]{kernel_event});
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					start_time = profiling_result[0];
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					end_time = profiling_result[0];
//					multiply_matrices_profiling += end_time - start_time;
//				}
				
				// STEP: calculate da vector = alpha matrix * beta vector
				clEnqueueNDRangeKernel(gpu._ocl_queue, calculate_da_vector, 1, null, new long[]{batch_size}, lwgs, 0, null, null);
				
//				// DEBUG: print intermediate results
//				if(DEBUG_MODE_ENABLED)
//				{
//					float[] da_vector_dbp = new float[batch_size * NUM_FITTING_PARAMETERS];
//					clEnqueueReadBuffer(gpu._ocl_queue, da_vector_buffer, CL_TRUE, 0, batch_size * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, Pointer.to(da_vector_dbp), 0, null, null);
//		
//					System.err.format("da vector = [%f, %f, %f, %f, %f, %f]\n", da_vector_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+0], da_vector_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+1], da_vector_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+2], da_vector_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+3], da_vector_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+4], da_vector_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+5]);
//				}
//				
//				if(PROFILING_MODE_ENABLED)
//				{
//					clWaitForEvents(1, new cl_event[]{kernel_event});
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					start_time = profiling_result[0];
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					end_time = profiling_result[0];
//					calculate_da_vector_profiling += end_time - start_time;
//				}
				
				// STEP: calculate updated parameters = fitted parameters + da vector
				clEnqueueNDRangeKernel(gpu._ocl_queue, calculate_updated_parameters, 1, null, new long[]{batch_size}, lwgs, 0, null, null);
				
//				// DEBUG: print intermediate results
//				if(DEBUG_MODE_ENABLED)
//				{
//					float[] incremented_params_dbp = new float[batch_size * NUM_FITTING_PARAMETERS];
//					clEnqueueReadBuffer(gpu._ocl_queue, updated_parameters_buffer, CL_TRUE, 0, batch_size * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, Pointer.to(incremented_params_dbp), 0, null, null);
//					
//					System.err.format("updated parameters = [%f, %f, %f, %f, %f, %f]\n", incremented_params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+0], incremented_params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+1], incremented_params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+2], incremented_params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+3], incremented_params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+4], incremented_params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+5]);
//				}
//				
//				if(PROFILING_MODE_ENABLED)
//				{
//					clWaitForEvents(1, new cl_event[]{kernel_event});
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					start_time = profiling_result[0];
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					end_time = profiling_result[0];
//					calculate_updated_parameters_profiling += end_time - start_time;
//				}
				
				// STEP: calculate chi2 for updated parameters
				clEnqueueNDRangeKernel(gpu._ocl_queue, calculate_updated_chi2, 1, null, new long[]{batch_size}, lwgs, 0, null, null);
				
//				// DEBUG: print intermediate results
//				float[] updated_chi2_dbp = new float[batch_size]; // NOTE: outside if statement on purpose
//				//clEnqueueReadBuffer(gpu._ocl_queue, updated_chi2_buffer, CL_TRUE, 0, batch_size * CL_DATATYPE_SIZE, Pointer.to(updated_chi2_dbp), 0, null, null); // TMP
//				if(DEBUG_MODE_ENABLED)
//				{
//					//float[] updated_chi2_dbp = new float[batch_size];
//					clEnqueueReadBuffer(gpu._ocl_queue, updated_chi2_buffer, CL_TRUE, 0, batch_size * CL_DATATYPE_SIZE, Pointer.to(updated_chi2_dbp), 0, null, null);
//					System.err.println("incremented chi2 = " + updated_chi2_dbp[DEBUG_TRACE_FIT_INDEX]);
//				}
				
//				if(PROFILING_MODE_ENABLED)
//				{
//					clWaitForEvents(1, new cl_event[]{kernel_event});
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					start_time = profiling_result[0];
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					end_time = profiling_result[0];
//					calculate_updated_chi2_profiling += end_time - start_time;
//				}
				
				// STEP: update parameters (conditional)
				clEnqueueNDRangeKernel(gpu._ocl_queue, update_parameters, 1, null, new long[]{batch_size}, lwgs, 0, null, null);
				
				// DEBUG: print intermediate fitted parameters
//				if(DEBUG_MODE_ENABLED)
//				{
//					if(current_chi2_dbp[DEBUG_TRACE_FIT_INDEX] > updated_chi2_dbp[DEBUG_TRACE_FIT_INDEX])
//					{
//						System.err.println("Fit got better!");
//						System.err.println("Abs.diff chi2 = " + Math.abs(current_chi2_dbp[DEBUG_TRACE_FIT_INDEX] - updated_chi2_dbp[DEBUG_TRACE_FIT_INDEX]));
//						
//						float[] updated_params_dbp = new float[batch_size * NUM_FITTING_PARAMETERS];
//						clEnqueueReadBuffer(gpu._ocl_queue, parameters_buffer, CL_TRUE, 0, batch_size * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, Pointer.to(updated_params_dbp), 0, null, null);
//						
//						System.err.format("current parameters = [%f, %f, %f, %f, %f, %f]\n", updated_params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+0], updated_params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+1], updated_params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+2], updated_params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+3], updated_params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+4], updated_params_dbp[DEBUG_TRACE_FIT_INDEX*NUM_FITTING_PARAMETERS+5]);
//					}
//					else
//					{
//						System.err.println("fit got worse...");
//					}
//				}
				
//				if(PROFILING_MODE_ENABLED)
//				{
//					clWaitForEvents(1, new cl_event[]{kernel_event});
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					start_time = profiling_result[0];
//					clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END, Sizeof.cl_long, Pointer.to(profiling_result), null);
//					end_time = profiling_result[0];
//					update_parameters_profiling += end_time - start_time;
//				}
				
				// increment iteration count
				++iteration_count;
			}
			
			// ******************************************************************************
			
			// set lambda to zero
			lambda_f = new float[batch_size]; // shorthand for initializing all elements to zero
			clEnqueueWriteBuffer(gpu._ocl_queue, lambda_buffer, true, 0, batch_size * CL_DATATYPE_SIZE, Pointer.to(lambda_f), 0, null, null);
			
			// recompute alpha matrix (with lambda equal to zero)
			clEnqueueNDRangeKernel(gpu._ocl_queue, calculate_alpha_matrix, 1, null, new long[]{batch_size}, lwgs, 0, null, null);
			
			// invert alpha matrix
			clEnqueueNDRangeKernel(gpu._ocl_queue, cholesky_decomposition, 1, null, new long[]{batch_size}, lwgs, 0, null, null); // step (a)
			clEnqueueNDRangeKernel(gpu._ocl_queue, forward_substitution, 1, null, new long[]{batch_size}, lwgs, 0, null, null); // step (b)
			clEnqueueNDRangeKernel(gpu._ocl_queue, multiply_matrices, 1, null, new long[]{batch_size}, lwgs, 0, null, null); // step (c)
			
			// calculate standard error of fitted parameters
			clEnqueueNDRangeKernel(gpu._ocl_queue, calculate_standard_error, 1, null, new long[]{batch_size}, lwgs, 0, null, null);
			
			// ******************************************************************************

			// use fitted parameters array in order to retain initial parameters
//			double[] fitted_parameters_d = new double[batch_size * NUM_FITTING_PARAMETERS];
			float[] fitted_parameters_f = new float[batch_size * NUM_FITTING_PARAMETERS];
			
			// read results back into fitted parameters buffer (NOTE: stalling operation)
			clEnqueueReadBuffer(gpu._ocl_queue, parameters_buffer, CL_TRUE, 0, batch_size * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, Pointer.to(fitted_parameters_f), 0, null, null);
			
			// read standard error of parameters
			float[] standard_errors_f = new float[batch_size * NUM_FITTING_PARAMETERS];
			clEnqueueReadBuffer(gpu._ocl_queue, standard_errors_buffer, CL_TRUE, 0, batch_size * NUM_FITTING_PARAMETERS * CL_DATATYPE_SIZE, Pointer.to(standard_errors_f), 0, null, null);
			
			// read chi-sqaured values
			float[] chi_squared_f = new float[batch_size * NUM_FITTING_PARAMETERS];
			clEnqueueReadBuffer(gpu._ocl_queue, current_chi2_buffer, CL_TRUE, 0, batch_size * CL_DATATYPE_SIZE, Pointer.to(chi_squared_f), 0, null, null);
			
			// ******************************************************************************
			
			parameters_index = 0;
			for(int i = ii; i < ii+batch_size && i < num_spots; ++i) // NOTE: check if less than num_spots
			{
				// correction of error estimation with scaling coefficient
				float error_coefficient = (float)Math.sqrt(chi_squared_f[i - ii]/(pixel_count-6));
				for (int j = 0; j < NUM_FITTING_PARAMETERS; j++)
				{
					standard_errors_f[parameters_index + j] *= error_coefficient;
				}
				
				// determine of spot is a false positive
				double false_positive = 0.0; // higher value is higher likelihood spot is a false positive
				
				//iterations didn't converge. suspicious point; NOTE: cannot be done with fixe number of iterations!
				//if (SMLlma.iterationCount== 101)
				//{
				//	nFalsePositive = 0.5;
				//}
				
				// penalize large and small spots
				double abs_xsd = Math.abs(fitted_parameters_f[parameters_index+4]);
				double abs_ysd = Math.abs(fitted_parameters_f[parameters_index+5]);
				if(abs_xsd > 1.3 * psf_sigma || abs_xsd < 0.7 * psf_sigma || abs_ysd > 1.3 * psf_sigma || abs_ysd < 0.7 * psf_sigma)
				{
					false_positive = 1.0;
				}
				
				// penalize low localization precision (i.e. localization error higher than psf)
				if(standard_errors_f[parameters_index+2] > psf_sigma || standard_errors_f[parameters_index+3] > psf_sigma)
				{
					false_positive = 1.0;
				}
				
				// calculate integrated spot intensity and estimation of SNR
				double integrated_amplitude = 0.0;
				double integrated_noise = 0.0;
				double average_noise = 0.0;
				double stdev_noise = 0.0;
				double snr_estimate = 0.0;
				
				if(false_positive < 0.5)
				{
					// get image processor for spot
					ImageProcessor improc = image.getStack().getProcessor((int)frame_numbers[i]);
					
					// center and range of fitted spot
					int x_centroid = (int)Math.round(fitted_parameters_f[parameters_index+2] + 0.5);
					int y_centroid = (int)Math.round(fitted_parameters_f[parameters_index+3] + 0.5);
					int mxsd = (int)Math.round(DOMConstants.FITRADIUS * psf_sigma);
					int mysd = (int)Math.round(DOMConstants.FITRADIUS * psf_sigma);
					
					// check if spot if within range of original image
					if(x_centroid - mxsd - 1 >= 0 && x_centroid + mxsd + 1 < image_width && y_centroid - mysd - 1 >= 0 && y_centroid + mysd + 1 < image_height)
					{
						//
						int x, y;
						
						// integrated spot intensity
						for(x = x_centroid - mxsd; x <= x_centroid + mxsd; ++x)
						{
							for(y = y_centroid - mysd; y <= y_centroid + mysd; ++y)
							{
								integrated_amplitude += improc.get(x, y);
							}
						}
						integrated_amplitude = integrated_amplitude - (fitted_parameters_f[parameters_index+0] * (2 * mxsd + 1) * (2 * mysd + 1));
						
						// average noise around spot
						//double sum_squared = 0.0;
						y = y_centroid - mysd - 1; // left border
						for(x = x_centroid - mxsd - 1; x <= x_centroid + mxsd + 1; ++x)
						{
							double px = improc.get(x, y);
							integrated_noise += px;
							//sum_squared += px * px;
						}
						y = y_centroid - mysd + 1; // right border
						for(x = x_centroid - mxsd - 1; x <= x_centroid + mxsd + 1; ++x)
						{
							double px = improc.get(x, y);
							integrated_noise += px;
							//sum_squared += px * px;
						}
						x = x_centroid - mxsd - 1; // top border
						for(y = y_centroid - mysd - 1; y <= y_centroid + mysd + 1; ++y)
						{
							double px = improc.get(x, y);
							integrated_noise += px;
							//sum_squared += px * px;
						}
						x = x_centroid - mxsd + 1; // bottom border
						for(y = y_centroid - mysd - 1; y <= y_centroid + mysd + 1; ++y)
						{
							double px = improc.get(x, y);
							integrated_noise += px;
							//sum_squared += px * px;
						}
						average_noise = integrated_noise / (4 * mxsd + 4 * mysd + 8);
						//sum_squared = sum_squared / (4 * mxsd + 4 * mysd + 8);
						//stdev_noise = Math.sqrt(sum_squared - (average_noise * average_noise));

						// standard deviation of noise
						y = y_centroid - mysd - 1; // left border
						for(x = x_centroid - mxsd - 1; x <= x_centroid + mxsd + 1; ++x)
						{
							double px = improc.get(x, y);
							stdev_noise += Math.pow(average_noise - px,2);
							//sum_squared += px * px;
						}
						y = y_centroid - mysd + 1; // right border
						for(x = x_centroid - mxsd - 1; x <= x_centroid + mxsd + 1; ++x)
						{
							double px = improc.get(x, y);
							stdev_noise += Math.pow(average_noise - px,2);
							//sum_squared += px * px;
						}
						x = x_centroid - mxsd - 1; // top border
						for(y = y_centroid - mysd - 1; y <= y_centroid + mysd + 1; ++y)
						{
							double px = improc.get(x, y);
							stdev_noise += Math.pow(average_noise - px,2);
							//sum_squared += px * px;
						}
						x = x_centroid - mxsd + 1; // bottom border
						for(y = y_centroid - mysd - 1; y <= y_centroid + mysd + 1; ++y)
						{
							double px = improc.get(x, y);
							stdev_noise += Math.pow(average_noise - px,2);
							//sum_squared += px * px;
						}
						stdev_noise = Math.sqrt(stdev_noise / (4 * mxsd + 4 * mysd + 8));
						
						// estimate signal-to-noise ratio
						snr_estimate = fitted_parameters_f[parameters_index+1] / stdev_noise;
					}
				}
				
				// add results to table (no lock needed) if spot is within image dimenions
				if(fitted_parameters_f[parameters_index+2] > 0 && fitted_parameters_f[parameters_index+2] < image_width && fitted_parameters_f[parameters_index+3] > 0 && fitted_parameters_f[parameters_index+2] < image_height)
				{
					// ignore false positives if requested
					if(!(ignore_false_positives && false_positive > 0.3))
					{
						res_table.incrementCounter();
						
						res_table.addValue("Amplitude_fit", fitted_parameters_f[parameters_index+1]);
						
						res_table.addValue("X_(px)", fitted_parameters_f[parameters_index+2]);
						res_table.addValue("Y_(px)", fitted_parameters_f[parameters_index+3]);
						res_table.addValue("X_(nm)", fitted_parameters_f[parameters_index+2] * pixel_size);
						res_table.addValue("Y_(nm)", fitted_parameters_f[parameters_index+3] * pixel_size);
						res_table.addValue("Z_(nm)", 0.0);
						res_table.addValue("False positive", false_positive);
						res_table.addValue("X_loc_error(px)", standard_errors_f[parameters_index+2]);
						res_table.addValue("Y_loc_error(px)", standard_errors_f[parameters_index+3]);
						
						res_table.addValue("BGfit", fitted_parameters_f[parameters_index+0]);
						res_table.addValue("IntegratedInt", integrated_amplitude);
						res_table.addValue("SNR", snr_estimate);
						
						res_table.addValue("chi2_fit", chi_squared_f[i - ii]);
						res_table.addValue("Frame Number", frame_numbers[i]);
						res_table.addValue("Iterations_fit", iterations);
						res_table.addValue("SD_X_fit_(px)", fitted_parameters_f[parameters_index+4]);
						res_table.addValue("SD_Y_fit_(px)", fitted_parameters_f[parameters_index+5]);
						res_table.addValue("Amp_loc_error", standard_errors_f[parameters_index+1]);
						
						// TODO: show particles in image as ROIs
					}
				}
				
				//  increment counters
				parameters_index+=NUM_FITTING_PARAMETERS;
			}
			
//			if(PROFILING_MODE_ENABLED)
//			{
//				clWaitForEvents(1, new cl_event[]{kernel_event});
//				clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, Sizeof.cl_long, Pointer.to(profiling_result), null);
//				start_time = profiling_result[0];
//				clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END, Sizeof.cl_long, Pointer.to(profiling_result), null);
//				end_time = profiling_result[0];
//				read_parameters_buffer_profiling += end_time - start_time;
//			}
			
//			// print results in CSV format
//			if(PRINT_RESULTS_IN_CSV_FORMAT && !SILENT_MODE_ENABLED)
//			{
//				if(USE_EXTENDED_CSV_FORMAT)
//				{
//					parameters_index = 0;
//					for(int i = ii; i < ii+batch_size && i < num_spots; ++i) // NOTE: check if less than num_spots
//					{
//						if(USE_DOUBLE_PRECISION)
//						{
//							System.out.format("%d,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f\n", i+1, fitted_parameters_d[parameters_index+0], fitted_parameters_d[parameters_index+1], fitted_parameters_d[parameters_index+2], fitted_parameters_d[parameters_index+3], fitted_parameters_d[parameters_index+4], fitted_parameters_d[parameters_index+5]);
//							parameters_index+=NUM_FITTING_PARAMETERS;
//						}
//						else
//						{
//							System.out.format("%d,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f\n", i+1, fitted_parameters_f[parameters_index+0], fitted_parameters_f[parameters_index+1], fitted_parameters_f[parameters_index+2], fitted_parameters_f[parameters_index+3], fitted_parameters_f[parameters_index+4], fitted_parameters_f[parameters_index+5]);
//							parameters_index+=NUM_FITTING_PARAMETERS;
//						}
//					}
//				}
//				else
//				{
//					parameters_index = 0;
//					for(int i = ii; i < ii+batch_size && i < num_spots; ++i) // NOTE: check if less than num_spots
//					{
//						if(USE_DOUBLE_PRECISION)
//						{
//							System.out.format("%d,%.15f,%.15f,0\n", i+1, fitted_parameters_d[parameters_index+2], fitted_parameters_d[parameters_index+3]);
//							parameters_index+=NUM_FITTING_PARAMETERS;
//						}
//						else
//						{
//							System.out.format("%d,%.15f,%.15f,0\n", i+1, fitted_parameters_f[parameters_index+2], fitted_parameters_f[parameters_index+3]);
//							parameters_index+=NUM_FITTING_PARAMETERS;
//						}
//					}
//				}
//			}
		}
		
//		// print profiling results
//		if(PROFILING_MODE_ENABLED)
//		{
//			int num_batches = (int)Math.ceil(num_spots / (double)batch_size);
//			double batch_overload = (Math.ceil(num_spots / (double)batch_size) * batch_size) / num_spots;
//			
//			System.err.println();
//			System.err.println("[ Profiling results ]");
//			System.err.println("Batch size = " + batch_size);
//			System.err.println("Local workgroup size = " + lwg_size);
//			System.err.println("Number of spots = " + num_spots);
//			System.err.println("Number of batches = " + num_batches);
//			System.err.println("Batch overload = " + batch_overload);
//			System.err.println("Number of iterations = " + MAX_ITERATIONS);
//			System.err.println();
//			
//			long total_execution_time = transform_data_profiling + write_image_data_buffer_profiling + write_parameters_buffer_profiling + write_lambda_buffer_profiling + calculate_current_chi2_profiling + calculate_alpha_matrix_profiling + calculate_beta_vector_profiling + cholesky_decomposition_profiling + forward_substitution_profiling + multiply_matrices_profiling + calculate_da_vector_profiling + calculate_updated_parameters_profiling + calculate_updated_chi2_profiling + update_parameters_profiling + read_parameters_buffer_profiling;
//			System.err.format("Transform data\t\t\t%.0f%%\t%.2f\n", 100 * transform_data_profiling / (float)total_execution_time, 1e-6 * transform_data_profiling);
//			System.err.format("Write image data buffer\t\t%.0f%%\t%.2f\n", 100 * write_image_data_buffer_profiling / (float)total_execution_time, 1e-6 * write_image_data_buffer_profiling);
//			System.err.format("Write parameters buffer\t\t%.0f%%\t%.2f\n", 100 * write_parameters_buffer_profiling / (float)total_execution_time, 1e-6 * write_parameters_buffer_profiling);
//			System.err.format("Write lambda buffer\t\t%.0f%%\t%.2f\n", 100 * write_lambda_buffer_profiling / (float)total_execution_time, 1e-6 * write_lambda_buffer_profiling);
//			System.err.format("Calculate current chi2\t\t%.0f%%\t%.2f\n", 100 * calculate_current_chi2_profiling / (float)total_execution_time, 1e-6 * calculate_current_chi2_profiling);
//			System.err.format("Calculate alpha matrix\t\t%.0f%%\t%.2f\n", 100 * calculate_alpha_matrix_profiling / (float)total_execution_time, 1e-6 * calculate_alpha_matrix_profiling);
//			System.err.format("Calculate beta vector\t\t%.0f%%\t%.2f\n", 100 * calculate_beta_vector_profiling / (float)total_execution_time, 1e-6 * calculate_beta_vector_profiling);
//			System.err.format("Cholesky decomposition\t\t%.0f%%\t%.2f\n", 100 * cholesky_decomposition_profiling / (float)total_execution_time, 1e-6 * cholesky_decomposition_profiling);
//			System.err.format("Forward substitution\t\t%.0f%%\t%.2f\n", 100 * forward_substitution_profiling / (float)total_execution_time, 1e-6 * forward_substitution_profiling);
//			System.err.format("Multiply matrices\t\t%.0f%%\t%.2f\n", 100 * multiply_matrices_profiling / (float)total_execution_time, 1e-6 * multiply_matrices_profiling);
//			System.err.format("Calculate da vector\t\t%.0f%%\t%.2f\n", 100 * calculate_da_vector_profiling / (float)total_execution_time, 1e-6 * calculate_da_vector_profiling);
//			System.err.format("Calculate updated parameters\t%.0f%%\t%.2f\n", 100 * calculate_updated_parameters_profiling / (float)total_execution_time, 1e-6 * calculate_updated_parameters_profiling);
//			System.err.format("Calculate updated chi2\t\t%.0f%%\t%.2f\n", 100 * calculate_updated_chi2_profiling / (float)total_execution_time, 1e-6 * calculate_updated_chi2_profiling);
//			System.err.format("Update parameters\t\t%.0f%%\t%.2f\n", 100 * update_parameters_profiling / (float)total_execution_time, 1e-6 * update_parameters_profiling);
//			System.err.format("Read fitted parameters\t\t%.0f%%\t%.2f\n", 100 * read_parameters_buffer_profiling / (float)total_execution_time, 1e-6 * read_parameters_buffer_profiling);
//			System.err.format("Total execution time\t\t%.0f%%\t%.2f\n", 100 * total_execution_time / (float)total_execution_time, 1e-6 * total_execution_time);
//			System.err.format("Normalised execution time\t\t%.2f\n", 1e-6 * total_execution_time / batch_overload);
//			System.err.println();
//		}
		
		// print execution time
		long end_time = System.nanoTime();
		double time_diff = (double)(end_time - start_time);
		String[] units = new String[]{"nanoseconds", "microseconds", "milliseconds", "seconds", "minutes", "hours", "way too long!"};
		int unit_index = 0;
		while(time_diff > 1e3f && unit_index < 3)
		{
			time_diff /= 1e3f;
			++unit_index;
		}
		while(time_diff > 60f && unit_index < 5)
		{
			time_diff /= 60f;
			++unit_index;
		}
		System.err.format("GPU fitting took %.1f %s\n", (time_diff), units[unit_index]);
		System.err.format("Averaging %.1f spots per second\n", num_spots/((end_time-start_time)/1e9f));
		
		// ******************************************************************************
		
		// release kernels
		clReleaseKernel(calculate_standard_error);
		clReleaseKernel(update_parameters);
		clReleaseKernel(calculate_updated_chi2);
		clReleaseKernel(calculate_updated_parameters);
		clReleaseKernel(calculate_da_vector);
		clReleaseKernel(multiply_matrices);
		clReleaseKernel(forward_substitution);
		clReleaseKernel(cholesky_decomposition);
		clReleaseKernel(calculate_beta_vector);
		clReleaseKernel(calculate_alpha_matrix);
		clReleaseKernel(calculate_current_chi2);
		
		// release buffers
		clReleaseMemObject(standard_errors_buffer);
		clReleaseMemObject(updated_parameters_buffer);
		clReleaseMemObject(da_vector_buffer);
		clReleaseMemObject(inverse_alpha_matrix_buffer);
		clReleaseMemObject(beta_vector_buffer);
		clReleaseMemObject(temp_matrix_buffer);
		clReleaseMemObject(alpha_matrix_buffer);
		clReleaseMemObject(updated_chi2_buffer);
		clReleaseMemObject(current_chi2_buffer);
		clReleaseMemObject(lambda_buffer);
		clReleaseMemObject(parameters_buffer);
		clReleaseMemObject(x_positions_buffer);
		clReleaseMemObject(y_positions_buffer);
		clReleaseMemObject(image_data_buffer);
	}
}
