
package DOM;

// Java
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

import java.io.IOException;
import java.io.BufferedReader;
import java.io.FileReader;

// ImageJ
import static ij.IJ.*;
import ij.gui.GenericDialog;

// OpenCL
import static org.jocl.CL.*;
import org.jocl.*;


/**
 *	TODO: add error checking code to each OpenCL function call
 */
public class GPUBase
{
	// configurable options
	private static boolean JOCL_USE_DOUBLE_PRECISION = false; // use double precision float point computations (more accurate, but likely to be slower)
	private static boolean AUTOMATIC_MODE_ENABLED = false; // enables automatic selection of best opencl device
	private static boolean PROFILING_MODE_ENABLED = false; // enables profiling on the OpenCL device
	private static boolean DEBUG_MODE_ENABLED = false; // enables printing of debug statements to standerd error console
	
	// configurable options (fixed for now)
	private static final boolean JOCL_THROW_EXCEPTIONS = true; // JOCL will throw CLException when errors occur
	private static final LogLevel JOCL_DEBUG_LEVEL = LogLevel.LOG_WARNING;//LogLevel.LOG_TRACE; // options are { LOG_QUIET, LOG_ERROR, LOG_WARNING, LOG_INFO, LOG_TRACE, LOG_DEBUG, LOG_DEBUGTRACE }
	private static final long JOCL_DEVICE_TYPE = CL_DEVICE_TYPE_GPU; // options are { CL_DEVICE_TYPE_DEFAULT | CL_DEVICE_TYPE_ALL | CL_DEVICE_TYPE_CPU | CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_ACCELERATOR }
	
	/*
		[preprocessor options]
		-D name 				sets macro to 1
		-D name=definition		sets macro to definition
		-I dir					additional include directories
		
		[math intrinsic options]
		-cl-single-precision-constant	treat double precision floating-point constant as single precision constant
		-cl-denorms-are-zero			performance hint; flush to zero
		
		[optimization options]
		-cl-opt-disable					disable all optimizations
		-cl-strict-aliasing				assue strict aliasing rule; i.e., no two pointers of different data type point to the same data
		-cl-mad-enable					a * b + c as single intruction, reduces precision
		-cl-no-signed-zeros				no +0 or -0, just 0; only if sign is not significant
		-cl-unsafe-math-optimizations	extreme optimizations, not adviced, breaks IEEE 754 standard
		-cl-finite-math-only			no NaN, -inf en +inf in floating points and doubles
		-cl-fast-relaxed-math			sets -cl-unsafe-math-optimizations and -cl-finite-math-only
		
		[warning options]
		-w						inhibit all warning messages
		-Werror					treat all warnings as errors
		
		[opencl c version]
		-cl-std=CL1.0|CL1.1		explicitely specify opencl standard
	*/
	
	// platform, devices and context
	protected cl_platform_id _ocl_platform_id;
	protected cl_device_id _ocl_device_id;
	protected cl_context _ocl_context;
	protected cl_command_queue _ocl_queue;
	protected cl_program _ocl_program;
	
	// ******************************************************
	
	/**
	 *	Constructor to GPU Base class that sets up an OpenCL context based on the device selected by the yser from the device selection dialog
	 *	@throws: a CLException when no suiable OpenCL device is available on the system, the user cancelled the device selection dialog, or the OpenCL context could not be initialised.
	 */
	public GPUBase() throws CLException
	{
		this(false); // use single precision
	}
	
	public GPUBase(boolean double_precision) throws CLException
	{
		this(double_precision, false); // do not use autmatic mode
	}
	
	public GPUBase(boolean double_precision, boolean auto_mode) throws CLException
	{
		this(double_precision, auto_mode, false); // do not use profiling mode
	}
	
	public GPUBase(boolean double_precision, boolean auto_mode, boolean profiling_mode) throws CLException
	{
		this(double_precision, auto_mode, profiling_mode, false); // do not use debug mode
	}
	
	public GPUBase(boolean double_precision, boolean auto_mode, boolean profiling_mode, boolean debug_mode) throws CLException
	{
		JOCL_USE_DOUBLE_PRECISION = double_precision;
		AUTOMATIC_MODE_ENABLED = auto_mode;
		PROFILING_MODE_ENABLED = profiling_mode;
		DEBUG_MODE_ENABLED = debug_mode;
		setup(); // set up OpenCL context
	}
	
	/**
	 *	Gathers information on OpenCL enabled platforms and devices and presents the user with a list of devices to selected from a dialog. A OpenCL context is created for the selected device and a work unit queue is constructed.
	 *	@throws: a CLException when no suiable OpenCL device is available on the system, the user cancelled the device selection dialog, or the OpenCL context could not be initialised.
	 */
	private void setup() throws CLException
	{
		// DEBUG print
		if(DEBUG_MODE_ENABLED) System.err.println("Setting up OpenCL platform");
		
		// Enable exception when requested
		setExceptionsEnabled(JOCL_THROW_EXCEPTIONS);
		setLogLevel(JOCL_DEBUG_LEVEL);
		
		// Obtain the number of platforms
		int numPlatforms[] = new int[1];
		clGetPlatformIDs(0, null, numPlatforms);
		
		// DEBUG print
		if(DEBUG_MODE_ENABLED) System.err.println("Number of platforms detected: "+numPlatforms[0]);
		
		// Throw exception if no OpenCL-enabled platform was found on this computer
		if(numPlatforms[0] == 0) throw new CLException("No OpenCL-enabled platform found!");
		
		// Obtain the platform IDs
		cl_platform_id platformIDs[] = new cl_platform_id[numPlatforms[0]];
		clGetPlatformIDs(platformIDs.length, platformIDs, null);
		
		// Collect all suitable devices from all available platforms
		List<cl_device_id> deviceIDs = new ArrayList<cl_device_id>();
		List<String> deviceNames = new ArrayList<String>();
		for(cl_platform_id platform_id : platformIDs)
		{
			
			// Obtain the number of devices for the current platform
			int numDevices[] = new int[1];
			try{
				clGetDeviceIDs(platform_id, JOCL_DEVICE_TYPE, 0, null, numDevices);
			}
			catch(CLException e)
			{
				// clGetDeviceIDs throws CLExcpetion "CL_DEVICE_NOT_FOUND" if no OpenCL devices that matched device_type were found"
				// solution: skip to next platform
				continue;
			}
			
			// DEBUG print
			if(DEBUG_MODE_ENABLED)
			{
				String platformName = getString(platform_id, CL_PLATFORM_NAME);
				System.err.println("Number of devices in platform " + platformName + ": " + numDevices[0]);
			}
			
			// Get all device ids for the current platform
			cl_device_id devicesArray[] = new cl_device_id[numDevices[0]];
			clGetDeviceIDs(platform_id, JOCL_DEVICE_TYPE, numDevices[0], devicesArray, null);
			
			// TODO: filter devices from list with select capabilities
			// E.g. Image support,
			//		Double precision
			//		Number work units
			//		OpenCL version
			//		Out-of-order execution queues?
			
			// add all capable devices to list
			deviceIDs.addAll(Arrays.asList(devicesArray));
		}
		
		// DEBUG print
		if(DEBUG_MODE_ENABLED) System.err.println("Number of OpenCL-enabled devices detected: " + deviceIDs.size());
		
		// look up device names for all capable devices
		for(cl_device_id device_id : deviceIDs)
		{
			String deviceName = getString(device_id, CL_DEVICE_NAME);
			deviceNames.add(deviceName);
		}
		
		// DEBUG print
		if(DEBUG_MODE_ENABLED)
		{
			System.err.println("List of capable devices:");
			for(String deviceName : deviceNames)
			{
				System.err.println("  - " + deviceName);
			}
		}
		
		// Throw exception if no OpenCL-enabled device was found for any of the platforms
		if(deviceIDs.size() == 0) throw new CLException("No OpenCL-enabled devices found!");
		
		String selectedDeviceName = deviceNames.get(0); // NOTE: for now, just select first available device
		if(AUTOMATIC_MODE_ENABLED)
		{
			// TODO: determine best device automatically
		}
		else
		{
			// Show device selection dialog
			GenericDialog gd = new GenericDialog("Select OpenCL-enabled device");
			gd.addChoice("Device", deviceNames.toArray(new String[0]), deviceNames.get(0)); // NOTE: type casting using toArray(T[])
			gd.setResizable(false);
			gd.showDialog(); // NOTE: modal dialog!
			
			// Throw exception if dialog was canceled by the user
			if(!gd.wasOKed())
			{
				throw new CLException("Device selection was canceled by user"); // canceled
			}
			
			// store selected device id and its corresponding platform id
			selectedDeviceName = gd.getNextChoice();
			
		}
		
		// look up device id and platform id for selected device name
		for(int i = 0; i < deviceNames.size(); ++i)
		{
			if(deviceNames.get(i) == selectedDeviceName)
			{
				// device id
				_ocl_device_id = deviceIDs.get(i);
				
				// platform id
				cl_platform_id platform[] = new cl_platform_id[1];
				clGetDeviceInfo(_ocl_device_id, CL_DEVICE_PLATFORM, Sizeof.cl_platform_id, Pointer.to(platform), null);
				_ocl_platform_id = platform[0];
				
				// DEBUG print
				if(DEBUG_MODE_ENABLED)
				{
					System.err.println("User selected device " + deviceNames.get(i));
					System.err.println("Device id = " + _ocl_device_id);
					System.err.println("Platform id = " + _ocl_platform_id);
				}
				
				// selected device found, break from for loop
				break;
			}
		}
		
		// create context
		cl_context_properties contextProperties = new cl_context_properties();
		contextProperties.addProperty(CL_CONTEXT_PLATFORM, _ocl_platform_id);
		_ocl_context = clCreateContext(contextProperties, 1, new cl_device_id[]{_ocl_device_id}, null, null, null);
		
		// create command queue
		long JOCL_QUEUE_PROPERTIES = (PROFILING_MODE_ENABLED ? CL_QUEUE_PROFILING_ENABLE : 0); // options are { 0 | CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE } // no support for out of order execution, yet
		_ocl_queue = clCreateCommandQueue(_ocl_context, _ocl_device_id, JOCL_QUEUE_PROPERTIES, null); // NOTE: each device should have its own command queue; for now, only single device is supported, so a single command queue will suffice
	}
	
	// ******************************************************
	
	protected boolean loadProgramFromString(String program) // TODO: add throws CLException when program could not be compiled
	{
		// DEBUG print
		if(DEBUG_MODE_ENABLED)
		{
			System.err.println("Program source:");
			System.err.println(program);
		}
		
		// TODO: check if program has been compiled before and load that binary instead
		//_ocl_program = clCreateProgramWithBinary();
		
		// otherwise, create program from source
		_ocl_program = clCreateProgramWithSource(_ocl_context, 1, new String[]{program}, new long[]{program.length()}, null);
		
		String JOCL_BUILD_PROGRAM_OPTIONS = (JOCL_USE_DOUBLE_PRECISION ? "-DENABLE_FP64 " : "-cl-single-precision-constant ") + (DEBUG_MODE_ENABLED ? "-cl-opt-disable" : "-cl-strict-aliasing -cl-mad-enable -cl-no-signed-zeros -cl-fast-relaxed-math -Werror");
		
		// DEBUG print
		if(DEBUG_MODE_ENABLED) System.err.println("Compiler options: " + JOCL_BUILD_PROGRAM_OPTIONS);
		
		// build program from source or binary
		try
		{
			clBuildProgram(_ocl_program, 1, new cl_device_id[]{_ocl_device_id}, JOCL_BUILD_PROGRAM_OPTIONS, null, null);
		}
		catch(CLException e)
		{
			if(DEBUG_MODE_ENABLED)
			{
				// DEBUG print
				System.err.println("Error: could not build program!");
				System.err.println(e.toString());
				
				// early exit
				return false;
			}
		}
		
		// unload compiler
		//clUnloadCompiler(); // WARN: deprecated in OpenCL 1.2
		
		return true;
	}
	
	protected boolean loadProgramFromFile(String filepath)
	{
		// open file
		BufferedReader inp;
		try
		{
			inp = new BufferedReader(new FileReader(filepath));
		}
		catch(IOException e)
		{
			// DEBUG print
			if(DEBUG_MODE_ENABLED) System.err.println("Could not load program file from " + filepath);
			
			// early exit
			return false;
		}
		
		// load contents of file
		String program = "";
		String line = null;
		try
		{
			while((line = inp.readLine()) != null)
			{
				// add line to program
				program += line + '\n';
			}
		}
		catch(IOException e)
		{
			// DEBUG print
			if(DEBUG_MODE_ENABLED) System.err.println("Error loading program from file!");
			
			// early exit
			return false;
		}
		
		// close file
		try
		{
			inp.close();
		}
		catch(IOException e)
		{
			// ignore, do nothing
		}
		
		// load program from string
		return loadProgramFromString(program);
	}
	
	// ******************************************************
	
	/**
	 * Returns the value of the device info parameter with the given name
	 *
	 * @param device The device
	 * @param paramName The parameter name
	 * @return The value
	 */
	private static String getString(cl_device_id device, int paramName)
	{
		// Obtain the length of the string that will be queried
		long size[] = new long[1];
		clGetDeviceInfo(device, paramName, 0, null, size);

		// Create a buffer of the appropriate size and fill it with the info
		byte buffer[] = new byte[(int)size[0]];
		clGetDeviceInfo(device, paramName, buffer.length, Pointer.to(buffer), null);

		// Create a string from the buffer (excluding the trailing \0 byte)
		return new String(buffer, 0, buffer.length-1);
	}

	/**
	 * Returns the value of the platform info parameter with the given name
	 *
	 * @param platform The platform
	 * @param paramName The parameter name
	 * @return The value
	 */
	private static String getString(cl_platform_id platform, int paramName)
	{
		// Obtain the length of the string that will be queried
		long size[] = new long[1];
		clGetPlatformInfo(platform, paramName, 0, null, size);

		// Create a buffer of the appropriate size and fill it with the info
		byte buffer[] = new byte[(int)size[0]];
		clGetPlatformInfo(platform, paramName, buffer.length, Pointer.to(buffer), null);

		// Create a string from the buffer (excluding the trailing \0 byte)
		return new String(buffer, 0, buffer.length-1);
	}
	
	// ******************************************************
	
	public static void main(String args[])
	{
		// check arguments
		if(args.length < 1)
		{
			System.err.println("You need to supply a OpenCL program source file as argument");
			System.exit(-1);
		}
		
		// create new instance of GPUBase
		GPUBase gpu = null;
		try
		{
			gpu = new GPUBase();
		}
		catch(CLException e)
		{
			System.err.println("Could not create new instance of GPUBase");
			System.err.println(e.toString());
			System.exit(-1);
		}
		
		// load program file
		if(gpu.loadProgramFromFile(args[0]))
		{
			System.err.println("New instance of GPUBase created successful for program " + args[0]);
		}
		
		System.err.println();
		System.exit(0);
	}
}
