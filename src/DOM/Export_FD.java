package DOM;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;

import ij.IJ;
import ij.io.SaveDialog;
import ij.plugin.PlugIn;
/** Export to Fairy Dust 3D webgl viewer **/
public class Export_FD implements PlugIn 
{
	SMLAnalysis sml = new SMLAnalysis();
	int nPatNumber; 
	
	public void run(String arg) {
		
		String filename;
		int i;
		String sFileWrite;
		DecimalFormat df2 = new DecimalFormat ("#.##");
		IJ.register(Export_FD.class);
		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.getHeadings()[0].equals("X_(px)"))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for export, please load one.");
					return;
		}
		
		//if the function is called from macro and filename already provided		
		boolean bCheck = !arg.isEmpty();	
		if(bCheck)
			filename = arg;
		//otherwise let's ask user for the filename
		else
		{
			
			filename = "FairyDust_export";
			
				
			SaveDialog sd = new SaveDialog("Export detections in Fairy Dust format", filename, ".xls");
		    String path = sd.getDirectory();
		    if (path==null)
		       	return;
			
		    filename = path+sd.getFileName();
		    IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
			IJ.log("Export to Fairy Dust format " + filename);
			double [] x_ =		sml.ptable.getColumnAsDoubles(DOMConstants.Col_Xnm);
			double [] y_ =		sml.ptable.getColumnAsDoubles(DOMConstants.Col_Ynm);
			double [] z_ = 		sml.ptable.getColumnAsDoubles(DOMConstants.Col_Znm);
			double [] xerr = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errX);
			double [] yerr =	sml.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errY);
			double [] zerr =	sml.ptable.getColumnAsDoubles(DOMConstants.Col_loc_errZ);
			double [] frame = 	sml.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);
			
			nPatNumber = x_.length;
			
			for (i=0;i<nPatNumber;i++)
        	{
				//recalculating uncertainity (as SD sum of both values) to one value and store it at xerr
        		xerr[i]= Math.sqrt(xerr[i]*xerr[i]+yerr[i]*yerr[i]);
        	}
	    	//exporting 		
	        try {        	
		       	File file = new File(filename);
				FileWriter writer = new FileWriter(file);
				//column headers				
				writer.write("frame\tx [nm]\ty [nm]\tz [nm]\tuncertainty_xy [nm]\tuncertainty_z [nm]\n");
				
				//writing values
				for (i=0;i<nPatNumber;i++)
		       	{
					
					//ignore uncertainty for now
					sFileWrite=Integer.toString((int)frame[i])+"\t"+df2.format(x_[i])+"\t"+df2.format(y_[i])+"\t"+df2.format(z_[i])+"\n";
					//sFileWrite=String.format("%d\t%.2f\t%.2f\t%.2f\t0.0\t0.0\n",(int)frame[i],x_[i],y_[i], z_[i]);
	
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
		}
	}

}
