package DOM;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;

/** Export to ThunderSTORM format 
 * **/
import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.SaveDialog;
import ij.plugin.PlugIn;

public class Export_TS implements PlugIn 
{
	SMLAnalysis sml = new SMLAnalysis();
	/** Photoelectrons per A/D count **/
	double dTSExportADU;
	/** Base level (A/D count) **/
	double dTSExportBase;
	
	/** if EMGAIN used or not**/
	boolean bEMGain;
	
	/** EMGAIN value**/
	double dEMGainValue;
	
	int nPatNumber; 
	
	public void run(String arg) {
		
		String filename;
		int i;
		String sFileWrite;
		DecimalFormat df2 = new DecimalFormat ("#.##");

		/** whether z values are calculated **/
		boolean zValsPresent=false;
		
		IJ.register(Export_TS.class);
		
		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.getHeadings()[0].equals("X_(px)"))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for export, please load one.");
					return;
		}
		
		GenericDialog dgExportTS = new GenericDialog("Export to ThunderSTORM format");
		dgExportTS.addNumericField("Photoelectrons per A/D count: ",Prefs.get("SiMoLoc.TSExportADU", 0.5),2,5," ");
		dgExportTS.addNumericField("Base level: ",Prefs.get("SiMoLoc.TSExportBase", 20),2,5,"A/D count");
		dgExportTS.addCheckbox("EMGain", Prefs.get("SiMoLoc.TSExportEMGain", false));
		dgExportTS.addNumericField("Value of EMGain: ",Prefs.get("SiMoLoc.TSExportEMGainValue", 100),2,5," ");
		dgExportTS.showDialog();
		if (dgExportTS.wasCanceled())
            return;
		
		dTSExportADU =  dgExportTS.getNextNumber();
		Prefs.set("SiMoLoc.TSExportADU", dTSExportADU);
		dTSExportBase =  dgExportTS.getNextNumber();
		Prefs.set("SiMoLoc.TSExportBase", dTSExportBase);
		bEMGain = dgExportTS.getNextBoolean();
		Prefs.set("SiMoLoc.TSExportEMGain", bEMGain);
		dEMGainValue =  dgExportTS.getNextNumber();
		Prefs.set("SiMoLoc.TSExportEMGainValue", dEMGainValue);
		
		
		//if the function is called from macro and filename already provided		
		boolean bCheck = !arg.isEmpty();	
		if(bCheck)
			filename = arg;
		//otherwise let's ask user for the filename
		else
		{
			
			filename = "TSexport";
			
				
			SaveDialog sd = new SaveDialog("Export detections in ThunderSTORM format", filename, ".csv");
		    String path = sd.getDirectory();
		    if (path==null)
		       	return;
			
		    filename = path+sd.getFileName();
		}
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		IJ.log("Export to ThunderSTORM format " + filename);
		IJ.log("Photoelectrons per A/D count: "+ String.format("%.2f",dTSExportADU));
		IJ.log("Base level (A/D count): "+ String.format("%.2f",dTSExportBase));
		if(bEMGain)
		{
			IJ.log("EMGain: on");
			IJ.log("EMGain value:"+ String.format("%.2f",dEMGainValue));
		}
		else
		{
			IJ.log("EMGain: off");
			dEMGainValue=1.0;
		}
		
		
		double [] x_ =		sml.ptable.getColumnAsDoubles(DOMConstants.Col_Xnm);
		double [] y_ =		sml.ptable.getColumnAsDoubles(DOMConstants.Col_Ynm);
		double [] z_ = 		sml.ptable.getColumnAsDoubles(DOMConstants.Col_Znm);
		double [] xerr = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errX);
		double [] yerr =	sml.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errY);
		double [] zerr =	sml.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errZ);
		double [] sdx_ =	sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_X);
		double [] sdy_ =	sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_Y);
		double [] bg =	 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_BGfit);
		double [] intInt = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_IntegrInt);
		double [] chi = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_chi);
		double [] snr = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_SNR);
		double [] amp = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_AmplFit);
		double [] frame = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);

		nPatNumber = amp.length;
		
		double [] bkgstd = new double[nPatNumber];
		double [] averSD = new double[nPatNumber];
		
		if(Math.abs(z_[0])>0.0)
			zValsPresent=true;
		

		//recalculating values to photons
        for (i=0;i<nPatNumber;i++)
        {
        	//already background subtracted
        	intInt[i]=intInt[i]*dTSExportADU/dEMGainValue;
        	//correct for the offset
        	bg[i]=(bg[i]-dTSExportBase)*dTSExportADU/dEMGainValue;
        	//estimate background noise from SNR
        	if(snr[i]>0.0)
        		bkgstd[i] = (amp[i]/snr[i])*dTSExportADU/dEMGainValue;
        	else
        		bkgstd[i]=0;
        	//recalculating uncertainty (as SD sum of both values) to one value and store it at xerr
        	xerr[i]= Math.sqrt(xerr[i]*xerr[i]+yerr[i]*yerr[i]);
        }
        if(!zValsPresent)
        {
        	for (i=0;i<nPatNumber;i++)
        		averSD[i]=0.5*(sdx_[i]+sdy_[i]);
        }
        
    	//exporting 		
        try {        	
	       	File file = new File(filename);
			FileWriter writer = new FileWriter(file);
			//column headers
			if(zValsPresent)
				writer.write("\"id\",\"frame\",\"x [nm]\",\"y [nm]\",\"z [nm]\",\"sigma1 [nm]\",\"sigma2 [nm]\",\"intensity [photon]\",\"offset [photon]\",\"bkgstd [photon]\",\"chi2\",\"uncertainty_xy [nm]\",\"uncertainty_z [nm]\"\n");
			else
				writer.write("\"id\",\"frame\",\"x [nm]\",\"y [nm]\",\"sigma [nm]\",\"intensity [photon]\",\"offset [photon]\",\"bkgstd [photon]\",\"chi2\",\"uncertainty_xy [nm]\"\n");
			
			//writing values
			for (i=0;i<nPatNumber;i++)
	       	{
				//id, frame, x, y
				//sFileWrite=String.format("%d,%d,%.1f,%.1f,",i+1,(int)frame[i],x_[i],y_[i]);
				sFileWrite=Integer.toString(i+1)+","+ Integer.toString((int)frame[i])+ ","+df2.format(x_[i])+","+df2.format(y_[i])+",";
				// z and sigmas
				if(zValsPresent)
				{
					//sFileWrite=sFileWrite+String.format("%.1f,%.1f,%.1f,",z_[i],sdx_[i],sdy_[i]);
					sFileWrite=sFileWrite+df2.format(z_[i])+","+df2.format(sdx_[i])+","+df2.format(sdy_[i])+",";
				}
				else
				{
					//sFileWrite=sFileWrite+String.format("%.1f,",averSD[i]);											
					sFileWrite=sFileWrite+df2.format(averSD[i])+",";
				}
				//intensity, offset, bkgstd, chi2, xy uncertainty
				//sFileWrite=sFileWrite+String.format("%.1f,%.1f,%.1f,%.1f,%.1f",intInt[i],bg[i],bkgstd[i],chi[i],xerr[i]);
				sFileWrite=sFileWrite+df2.format(intInt[i])+","+df2.format(bg[i])+","+df2.format(bkgstd[i])+","+df2.format(chi[i])+","+df2.format(xerr[i]);
				if(zValsPresent)
					sFileWrite=sFileWrite+","+df2.format(zerr[i]);
					//sFileWrite=sFileWrite+String.format(",%.1f",zerr[i]);
				sFileWrite=sFileWrite+"\n";
				writer.write(sFileWrite);
				IJ.showProgress(i, nPatNumber);
	       	}
						
			writer.close();
        } catch (IOException e) {	
			IJ.error(e.getMessage());
			//e.printStackTrace();
		}
        IJ.showProgress(nPatNumber, nPatNumber);
        IJ.log("Exporting results... Done.");
        IJ.showStatus("Exporting results... Done.");
		return;
	}

}
