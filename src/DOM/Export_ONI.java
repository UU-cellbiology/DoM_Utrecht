package DOM;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;

import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.SaveDialog;
import ij.plugin.PlugIn;

public class Export_ONI implements PlugIn 
{
	SMLAnalysis sml = new SMLAnalysis();
	int nDataType=0;
	double dZnmPxSize=0;
	int nPatNumber; 
	
	@Override
	public void run(String arg) {
		// TODO Auto-generated method stub
		
		String filename;
		int i;
		String sFileWrite;
		DecimalFormatSymbols symbols = new DecimalFormatSymbols();
		symbols.setDecimalSeparator('.');
		DecimalFormat df2 = new DecimalFormat ("#.#####", symbols);
		
		IJ.register(Export_ONI.class);
		
		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.getHeadings()[0].equals("X_(px)"))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for export, please load one.");
					return;
		}
		String [] sDataType = new String [] {
				"ONI","DoM"};
		GenericDialog dgExportONI = new GenericDialog("Export to ONI format");
		dgExportONI.addChoice("Results table comes from:", sDataType, Prefs.get("SiMoLoc.ExportONI", "ONI"));
		dgExportONI.addNumericField("Pixel size in Z:",Prefs.get("SiMoLoc.ExportONIZnm", 0),1,5,"nm");
		dgExportONI.showDialog();
		if (dgExportONI.wasCanceled())
            return;
		
		nDataType = dgExportONI.getNextChoiceIndex();
		Prefs.set("SiMoLoc.ExportONI", sDataType[nDataType]);
		dZnmPxSize =  dgExportONI.getNextNumber();
		Prefs.set("SiMoLoc.ExportONIZnm", dZnmPxSize);
		
		//if the function is called from macro and filename already provided		
		boolean bCheck = !arg.isEmpty();	
		if(bCheck)
			filename = arg;
		//otherwise let's ask user for the filename
		else
		{
			
			filename = "ONIexport";
			
				
			SaveDialog sd = new SaveDialog("Export detections in ONI format", filename, ".csv");
		    String path = sd.getDirectory();
		    if (path==null)
		       	return;
			
		    filename = path+sd.getFileName();
		}
		
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		IJ.log("Export to ONI format " + filename);
		
		
		double [] x_ =		sml.ptable.getColumnAsDoubles(DOMConstants.Col_Xnm);
		double [] y_ =		sml.ptable.getColumnAsDoubles(DOMConstants.Col_Ynm);
		double [] z_ = 		sml.ptable.getColumnAsDoubles(DOMConstants.Col_Znm);
		double [] xpx_ =		sml.ptable.getColumnAsDoubles(DOMConstants.Col_X);
		double [] ypx_ =		sml.ptable.getColumnAsDoubles(DOMConstants.Col_Y);		
		double [] xerr = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errX);
		double [] yerr =	sml.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errY);
		double [] bg =	 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_BGfit);
		double [] intInt = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_IntegrInt);
		double [] frame = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);
		double [] channel = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_BGfit_error);
		
		nPatNumber = bg.length;
		// take into account that in ONI numbering of frames starts from zero,
		// while in DoM from 1
		for( i=0;i<nPatNumber;i++)
		{
			frame[i]=frame[i]-1;
		}
		
		double [] xprpx = new double[nPatNumber];
		double [] yprpx = new double[nPatNumber];
		double [] zpx = new double[nPatNumber];
		
		//determine pixel size for recalculation of localization precision to pixels for ONI
		double pxSizeinNm= sml.ptable.getValueAsDouble(DOMConstants.Col_Xnm, 0)/sml.ptable.getValueAsDouble(DOMConstants.Col_X, 0);
		pxSizeinNm=Math.round(pxSizeinNm);
		
    	//exporting 		
        try {        	
	       	File file = new File(filename);
			FileWriter writer = new FileWriter(file);
			//column headers
			
			writer.write("Channel,Frame,X (nm),Y (nm),Z (nm),X precision (nm),Y precision (nm),X (pix),Y (pix),Z (pix),X precision (pix),Y precision (pix),Photons,Background\n");

			
			//writing values
			for (i=0;i<nPatNumber;i++)
	       	{
				
				//channel
				//ONI imported data
				if(nDataType==0)
				{
					//channel # stored in Col_BGfit_error
					sFileWrite=Integer.toString((int)channel[i])+",";
				}
				//DoM data
				else
				{
					sFileWrite="0,";
				}
				//frame
				sFileWrite=sFileWrite+Integer.toString((int)frame[i])+",";
				
				//x, y, z in nm
				
				sFileWrite=sFileWrite+df2.format(x_[i])+","+df2.format(y_[i])+","+df2.format(z_[i])+",";
				
				//xy precision and xy in px
				sFileWrite=sFileWrite+df2.format(xerr[i])+","+df2.format(yerr[i])+","+df2.format(xpx_[i])+","+df2.format(ypx_[i])+",";
				
				//z in px
				if(dZnmPxSize==0)
				{
					sFileWrite=sFileWrite+"0,";
				}
				else
				{
					sFileWrite=sFileWrite+df2.format(z_[i]/dZnmPxSize)+",";
				}
				
				//xy precision in pixel
				sFileWrite=sFileWrite+df2.format(xerr[i]/pxSizeinNm)+","+df2.format(yerr[i]/pxSizeinNm)+",";
				
				//photons and background
				sFileWrite=sFileWrite+df2.format(intInt[i])+","+df2.format(bg[i])+"\n";
			
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
