package DOM;

//import jaolho.data.lma.LMA;
import jaolho.data.lma.LMAMultiDimFunction;

//public class SMLGFit {

		/** Two-dimensional gaussian function with a form of y = a0 + a1*exp-((x1-a2)^2/(2*a4) +  (x2-a3)^2/(2*a5))*/
		// a[0] - Background
		// a[1] - Amplitude
	 	// a[2] - x0
		// a[3] - y0
		// a[4] - sigmax
		// a[5] - sigmay
	
		public class SMLTwoDGaussian extends LMAMultiDimFunction {
			@Override
			public double getY(double x[], double[] a) {
				return a[0] + a[1] * ( Math.exp(-0.5* ( (x[0]-a[2])*(x[0]-a[2])/(a[4]*a[4])    +  (x[1]-a[3])*(x[1]-a[3])/(a[5]*a[5])    ) ) );
			}
		/** Derivatives of two-dimensional gaussian function with respect to parameters **/
			@Override
			public double getPartialDerivate(double x[], double[] a, int parameterIndex) {
				switch (parameterIndex) {
					case 0: return 1;
					case 1: return (Math.exp(-0.5* ( (x[0]-a[2])*(x[0]-a[2])/(a[4]*a[4])    +  (x[1]-a[3])*(x[1]-a[3])/(a[5]*a[5])    ) ) );
					case 2: return a[1] * ( Math.exp(-0.5* ( (x[0]-a[2])*(x[0]-a[2])/(a[4]*a[4])    +  (x[1]-a[3])*(x[1]-a[3])/(a[5]*a[5])    ) ) )*(x[0]-a[2])/(a[4]*a[4]);
					case 3: return a[1] * ( Math.exp(-0.5* ( (x[0]-a[2])*(x[0]-a[2])/(a[4]*a[4])    +  (x[1]-a[3])*(x[1]-a[3])/(a[5]*a[5])    ) ) )*(x[1]-a[3])/(a[5]*a[5]);
					case 4: return a[1] * ( Math.exp(-0.5* ( (x[0]-a[2])*(x[0]-a[2])/(a[4]*a[4])    +  (x[1]-a[3])*(x[1]-a[3])/(a[5]*a[5])    ) ) )*(x[0]-a[2])*(x[0]-a[2])/(a[4]*a[4]*a[4]);
					case 5: return a[1] * ( Math.exp(-0.5* ( (x[0]-a[2])*(x[0]-a[2])/(a[4]*a[4])    +  (x[1]-a[3])*(x[1]-a[3])/(a[5]*a[5])    ) ) )*(x[1]-a[3])*(x[1]-a[3])/(a[5]*a[5]*a[5]);
				}
				throw new RuntimeException("No such parameter index: " + parameterIndex);
			}
		}
//}
	
