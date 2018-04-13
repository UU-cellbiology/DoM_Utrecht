package DOM;

import ij.IJ;
import ij.Prefs;
import ij.plugin.PlugIn;

public class Z_Calculation implements PlugIn
{
	SMLAnalysis sml = new SMLAnalysis();
	double zfitRangeMin,zfitRangeMax;
	double [] fitCoefZ,fitCoefX,fitCoefY;
	
	double []  sdx, sdy, z, x,y, xpx, ypx, sdx_err, sdy_err, z_err;
	@Override
	public void run(String arg) 
	{
		int i, nParticles;
		double dDiff, dVal, dVal1, dVal2, dVal3, dVal_Err;
		// TODO Auto-generated method stub
		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.getHeadings()[0].equals("X_(px)"))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for reconstruction, please load one.");
			return;
		}
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		IJ.log("Calculating Z values");
		fitCoefZ = new double[4];
		
		//read stored values of calibration
		
		fitCoefZ[0]=Prefs.get("SiMoLOc.ZC_polCoef0", Double.NaN);
		if(Double.isNaN(fitCoefZ[0]))
		{
			IJ.error("There is no valid Z-calibration present, please make one!");
			return;
		}
		fitCoefZ[1]=Prefs.get("SiMoLOc.ZC_polCoef1", 0);
		fitCoefZ[2]=Prefs.get("SiMoLOc.ZC_polCoef2", 0);
		fitCoefZ[3]=Prefs.get("SiMoLOc.ZC_polCoef3", 0);
		zfitRangeMin=Prefs.get("SiMoLOc.ZC_fitRangeMin", 0);
		zfitRangeMax=Prefs.get("SiMoLOc.ZC_fitRangeMax", 1000);
		IJ.log("Calibration curve range: ");
		IJ.log("Z min: " + Double.toString(zfitRangeMin) + " nm");
		IJ.log("Z max: " + Double.toString(zfitRangeMax) + " nm");
		IJ.log("Polynomial coeff 0: " + Double.toString(fitCoefZ[0]));
		IJ.log("Polynomial coeff 1: " + Double.toString(fitCoefZ[1]));
		IJ.log("Polynomial coeff 2: " + Double.toString(fitCoefZ[2]));
		IJ.log("Polynomial coeff 3: " + Double.toString(fitCoefZ[3]));
		IJ.log("Calibration creation date:" + Prefs.get("SiMoLOc.ZC_calDate","not available"));
		
		sdx = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_X);
		sdy = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_Y);
		sdx_err = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_X_err);
		sdy_err = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_Y_err);
		IJ.showStatus("Calculating Z coordinate...");
		nParticles=sdx.length;
		z = new double [nParticles];
		z_err = new double [nParticles];
		for(i=0;i<nParticles;i++)
		{
			dDiff = sdx[i]-sdy[i];
			//dVal=dDiff;
			dVal = fitCoefZ[0] + fitCoefZ[1]*dDiff +fitCoefZ[2]*dDiff*dDiff+fitCoefZ[3]*dDiff*dDiff*dDiff;
			
			if(dVal<zfitRangeMin)
			{
				dVal=zfitRangeMin;
			}
			if(dVal>zfitRangeMax)
			{
				dVal=zfitRangeMax;
			}
			z[i]=dVal;
			//error of difference
			dVal_Err = Math.sqrt(sdx_err[i]*sdx_err[i]+sdy_err[i]*sdy_err[i]);
			//summary error according to error propagation: first variance
			dVal_Err = Math.pow(fitCoefZ[1]*dVal_Err,2) + Math.pow(fitCoefZ[2]*2*dVal_Err/dDiff,2) + Math.pow(fitCoefZ[3]*3*dVal_Err/(dDiff*dDiff),2);
			z_err[i]=Math.sqrt(dVal_Err);
			if(z_err[i]==0)
			{
				z_err[i]++;
				z_err[i]--;
			}				
			IJ.showProgress(i+1, nParticles);
		}
		//account for XY wobbling
		if(Prefs.get("SiMoLOc.ZC_XYWobbling", false))
		{
			fitCoefX = new double[5];
			fitCoefY = new double[5];
			//get coefficients
			fitCoefX[0]=Prefs.get("SiMoLOc.ZC_polCoefX0", 0);
			fitCoefX[1]=Prefs.get("SiMoLOc.ZC_polCoefX1", 0);
			fitCoefX[2]=Prefs.get("SiMoLOc.ZC_polCoefX2", 0);
			fitCoefX[3]=Prefs.get("SiMoLOc.ZC_polCoefX3", 0);
			fitCoefX[4]=Prefs.get("SiMoLOc.ZC_polCoefX4", 0);
			
			fitCoefY[0]=Prefs.get("SiMoLOc.ZC_polCoefY0", 0);
			fitCoefY[1]=Prefs.get("SiMoLOc.ZC_polCoefY1", 0);
			fitCoefY[2]=Prefs.get("SiMoLOc.ZC_polCoefY2", 0);
			fitCoefY[3]=Prefs.get("SiMoLOc.ZC_polCoefY3", 0);
			fitCoefY[4]=Prefs.get("SiMoLOc.ZC_polCoefY4", 0);
			//get x and y values
			x = sml.ptable.getColumnAsDoubles(DOMConstants.Col_Xnm);
			y = sml.ptable.getColumnAsDoubles(DOMConstants.Col_Ynm);
			for(i=0;i<nParticles;i++)
			{
				dVal=(z[i]-fitCoefX[4]);
				x[i] -= fitCoefX[0] + fitCoefX[1]*dVal +fitCoefX[2]*dVal*dVal+fitCoefX[3]*dVal*dVal*dVal;
				y[i] -= fitCoefY[0] + fitCoefY[1]*dVal +fitCoefY[2]*dVal*dVal+fitCoefY[3]*dVal*dVal*dVal;;
			}
			
		}
		IJ.showStatus("Updating table with calculated Z values...");
		//update table
		sml.ptable_lock.lock();
		for(i=0;i<nParticles;i++)
		{
			sml.ptable.setValue(DOMConstants.Col_Znm, i, z[i]);
			sml.ptable.setValue(DOMConstants.Col_loc_errZ, i, z_err[i]);
			IJ.showProgress(i+1, nParticles);
		}
		if(Prefs.get("SiMoLOc.ZC_XYWobbling", false))
		{
			double pxsize =  sml.ptable.getValueAsDouble(DOMConstants.Col_X, 0)/sml.ptable.getValueAsDouble(DOMConstants.Col_Xnm, 0);
			for(i=0;i<nParticles;i++)
			{
				sml.ptable.setValue(DOMConstants.Col_Xnm, i, x[i]);
				sml.ptable.setValue(DOMConstants.Col_Ynm, i, y[i]);
				sml.ptable.setValue(DOMConstants.Col_X, i, x[i]*pxsize);
				sml.ptable.setValue(DOMConstants.Col_Y, i, y[i]*pxsize);
			}
		}
		sml.ptable_lock.unlock();
		//sml_.ptable.updateResults();
		sml.showTable();
		IJ.showStatus("Updating table with calculated Z values...done");
	} 

}
