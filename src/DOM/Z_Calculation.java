package DOM;

import ij.IJ;
import ij.Prefs;
import ij.plugin.PlugIn;

public class Z_Calculation implements PlugIn
{
	SMLAnalysis sml = new SMLAnalysis();
	double zfitRangeMin,zfitRangeMax;
	double [] fitCoef;
	
	double []  sdx, sdy, z;
	@Override
	public void run(String arg) 
	{
		int i, nParticles;
		double dVal;
		// TODO Auto-generated method stub
		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.getHeadings()[0].equals("X_(px)"))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for reconstruction, please load one.");
			return;
		}
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		IJ.log("Calculating Z values");
		fitCoef = new double[4];
		
		//read stored values of calibration
		
		fitCoef[0]=Prefs.get("SiMoLOc.ZC_polCoef0", Double.NaN);
		if(fitCoef[0]==Double.NaN)
		{
			IJ.error("There is no valid Z-calibration present, please make one!");
			return;
		}
		fitCoef[1]=Prefs.get("SiMoLOc.ZC_polCoef1", 0);
		fitCoef[2]=Prefs.get("SiMoLOc.ZC_polCoef2", 0);
		fitCoef[3]=Prefs.get("SiMoLOc.ZC_polCoef3", 0);
		zfitRangeMin=Prefs.get("SiMoLOc.ZC_fitRangeMin", 0);
		zfitRangeMax=Prefs.get("SiMoLOc.ZC_fitRangeMax", 1000);
		IJ.log("Calibration curve range: ");
		IJ.log("Z min: " + Double.toString(zfitRangeMin) + " nm");
		IJ.log("Z max: " + Double.toString(zfitRangeMax) + " nm");
		IJ.log("Polynomial coeff 0: " + Double.toString(fitCoef[0]));
		IJ.log("Polynomial coeff 1: " + Double.toString(fitCoef[1]));
		IJ.log("Polynomial coeff 2: " + Double.toString(fitCoef[2]));
		IJ.log("Polynomial coeff 3: " + Double.toString(fitCoef[3]));
		IJ.log("Calibration creation date:" + Prefs.get("SiMoLOc.ZC_calDate","not available"));
		
		sdx = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_X);
		sdy = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_Y);
		IJ.showStatus("Calculating Z coordinate...");
		nParticles=sdx.length;
		z = new double [nParticles];
		for(i=0;i<nParticles;i++)
		{
			dVal = sdx[i]-sdy[i];
			dVal = fitCoef[0] + fitCoef[1]*dVal +fitCoef[2]*dVal*dVal+fitCoef[3]*dVal*dVal*dVal;
			
			if(dVal<zfitRangeMin)
			{
				dVal=zfitRangeMin;
			}
			if(dVal>zfitRangeMax)
			{
				dVal=zfitRangeMax;
			}
			z[i]=dVal;
			IJ.showProgress(i+1, nParticles);
		}
		IJ.showStatus("Updating table with calculated Z values...");
		//update table
		sml.ptable_lock.lock();
		for(i=0;i<nParticles;i++)
		{
			sml.ptable.setValue(DOMConstants.Col_Znm, i, z[i]);
			IJ.showProgress(i+1, nParticles);
		}
		sml.ptable_lock.unlock();
		//sml_.ptable.updateResults();
		sml.showTable();
		IJ.showStatus("Updating table with calculated Z values...done");
	} 

}
