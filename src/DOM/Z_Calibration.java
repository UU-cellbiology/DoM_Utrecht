package DOM;
//in detect molecules on line 329
//add width and height of original image to the results table
//sml.ptable.addValue("Original_image_size",imp.getWidth());
//sml.ptable.addValue("Original_image_size",imp.getHeight());


import java.awt.AWTEvent;
import java.awt.Button;
import java.awt.Color;
import java.awt.Panel;
import java.awt.Toolkit;
import java.awt.event.AWTEventListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
//import ij.gui.ImageCanvas;
//import ij.gui.NonBlockingGenericDialog;
//import ij.gui.Plot;
import ij.measure.CurveFitter;
import ij.plugin.PlugIn;
import ij.plugin.frame.Fitter;

public class Z_Calibration implements PlugIn{

	SMLDialog dlg = new SMLDialog();
	SMLAnalysis sml = new SMLAnalysis();
	int  maxTracklength;
	
	double []  f, sdx, sdx_err, sdy,sdy_err, trackid, particleid, tl,  sortedTracklength, sdDif, rsquared;
	double[][] longestTrack;
	
	double maxLocError = 0.5;
	
	ImagePlus plotIp;
	
	public void run(String arg) 
	{
		IJ.register(Z_Calibration.class);
	
		if (!dlg.zCalibration()) return;
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		

		switch(dlg.sZcUse)
		{
			case 0: 
				IJ.log("Creating Z calibration from image stack");
				IJ.log("Step 1: Detecting molecules");
				IJ.run("Detect Molecules");
				IJ.log("Step 2: Linking detections");
				IJ.run("Link Particles to Tracks", "for=[All particles] max=4 measure=[Initial position] maximum=10 display show");
				IJ.log("Step 3: Creating Z calibration from Results table");
				zCal_particleTable(); 
				//zCal_imageStack(); 
				break;
			case 1: 
				IJ.log("Creating Z calibration from Results table");
				zCal_particleTable(); 
				
				break;
			case 2:
				IJ.log("Creating Z calibration from user provided coefficients");
				zCal_polCoef(); 
				break;
		}
	}
		
	
	/** making z-calibration from Particles table,
	 * link particles if data is not provided*/
	public void zCal_particleTable()
	{
		
		int i;
		Plot p;
		ImagePlus ip;
		
		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.getHeadings()[0].equals("X_(px)"))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for reconstruction, please load one.");
			return;
		}
		//check that the particles are linked
		//TODO add linking stage here
		if ( !sml.ptable.columnExists(DOMConstants.Col_TrackID))
		{
			//IJ.error("Particle table is probably not linked.");
			IJ.log("WARNING! Particles are not linked in the table, performing automatic linking.");
			IJ.run("Link Particles to Tracks", "for=[All particles] max=4 measure=[Initial position] maximum=10");;
			sml = new SMLAnalysis();
			//return;
		}
		
		IJ.log("Polynomial degree: "+Integer.toString(dlg.fitPolynomialDegree));
		IJ.log("Distance between Z planes: "+Double.toString(dlg.zCalDistBetweenPlanes)+" nm");
		IJ.log("R^2 threshold: "+ Double.toString(dlg.zCalRsquareThreshold));
		
		//first, sort data by particle and trackN
		Sort_Results.sorting_external_silent(sml, DOMConstants.Col_ParticleID, true);
		Sort_Results.sorting_external_silent(sml, DOMConstants.Col_TrackID, true);
		
		//get values
		sdx = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_X);
		sdx_err = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_X_err);
		sdy = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_Y);
		sdy_err = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_Y_err);
		f = sml.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);
		trackid = sml.ptable.getColumnAsDoubles(DOMConstants.Col_TrackID);
		particleid = sml.ptable.getColumnAsDoubles(DOMConstants.Col_ParticleID);
		rsquared = sml.ptable.getColumnAsDoubles(DOMConstants.Col_chi);
		
		//find longest track taking into account r2 threshold
		int dLongestTrackID = (int)trackid[0];
		maxTracklength = 0;
		int dCurrentLength = 0;
		double dCurrentTrack = trackid[0];
		
		for(i=0;i<sdx.length;i++)
		{
			if((int)trackid[i]==(int)dCurrentTrack)
			{
				if(rsquared[i]>=dlg.zCalRsquareThreshold)
					{dCurrentLength++;}
			}
			else
			{
				if(dCurrentLength>maxTracklength)
				{
					maxTracklength=dCurrentLength;
					dLongestTrackID=(int)dCurrentTrack;
				}
				dCurrentLength=1;
				dCurrentTrack=trackid[i];
			}
		}
		//check that track has enough point
		if(maxTracklength<10)
		{
			IJ.error("Error! Longest track is less than 10 points. Adjust linking, decrease z-step or choose different calibration stack.");
			return;			
		}
		IJ.log("Longest track length after filtering: "+Integer.toString(maxTracklength) + " particles");
		
		//fill in data from the longest track
		longestTrack = new double[5][maxTracklength];
		dCurrentLength = 0;
		for(i=0;i<sdx.length;i++)
		{
				if((int)trackid[i]==(int)dLongestTrackID)
				{
					if(rsquared[i]>=dlg.zCalRsquareThreshold)
					{				
						longestTrack[0][dCurrentLength]=sdx[i];
						longestTrack[1][dCurrentLength]=sdy[i];
						longestTrack[2][dCurrentLength]=f[i];
						longestTrack[3][dCurrentLength]=sdx_err[i];
						longestTrack[4][dCurrentLength]=sdy_err[i];
						dCurrentLength++;
					}
				}
		}
		//subtract first frame
		for(i=1;i<maxTracklength;i++)
		{
			longestTrack[2][i]=(longestTrack[2][i]-longestTrack[2][0])*dlg.zCalDistBetweenPlanes;
		}
		//first frame
		longestTrack[2][0]=0;
		
		//calculate difference between SDx and SDy
		sdDif = new double[maxTracklength];

		for(i = 0; i<maxTracklength;i++)
		{
			sdDif[i] = longestTrack[0][i] - longestTrack[1][i];
			
		}
		
		//Estimate position of zero
		//let's take middle third
		int nBeg = (int)Math.round(maxTracklength/3.0);
		int nEnd = (int)Math.round(2.0*maxTracklength/3.0);
		double zOffset = 0;
		double dMin = Double.MAX_VALUE;
		for (i=nBeg;i<=nEnd;i++)
		{
			
			if(Math.abs(sdDif[i])<dMin)
			{
				zOffset = longestTrack[2][i];
				dMin = sdDif[i];
			}
			
		}
		//ok, let's see if user-provided estimation of Z required

		if(dlg.zOverride)
		{
			p = new Plot("Width and height of PSF","Z position (nm)","size of PSF (nm)");
			//plot sdx
			p.setColor(Color.red);
			//p.addPoints(longestTrack[2], longestTrack[0], Plot.CONNECTED_CIRCLES);	
			p.addPoints(longestTrack[2], longestTrack[0],longestTrack[3], Plot.CIRCLE);	
			//plot sdy
			p.setColor(Color.blue);
			p.addPoints(longestTrack[2], longestTrack[1],longestTrack[4], Plot.CIRCLE);
			p.setColor(Color.black);
	
			p.setLegend("width	height", Plot.AUTO_POSITION);
			ip = p.getImagePlus();
			ip.show();
	

			zOffset = dlgSetZero(zOffset);
			ip.close();
		
		}
		IJ.log("Z offset shift: "+Double.toString(zOffset) + " nm");
		
		//subtract z-offset from all z-values
		for(i = 0; i < maxTracklength; i++)
		{
			longestTrack[2][i] = longestTrack[2][i] - zOffset;  
		}

		IJ.log("Fitting results: ");
		//let's show plot and define range of fitting
		p = new Plot("Difference in PSF width and height","Z position (nm)","width-height (nm)");
		p.addPoints(longestTrack[2], sdDif, Plot.CONNECTED_CIRCLES);
		ip = p.getImagePlus();
		ip.show();
		

		//fitPars = fitRangeMin, fitRangeMax, polDegree
		int[] fitPars = dlgFitRange((int)longestTrack[2][0],(int)longestTrack[2][maxTracklength-1]);
		ip.close();
		if(fitPars[0] == 0 & fitPars[1] == 0)
			return;
		
		int indexMin = findIndex(longestTrack[2],fitPars[0]);
		int indexMax = findIndex(longestTrack[2],fitPars[1]);
		
		
		CurveFitter cf = polyFit(Arrays.copyOfRange(sdDif,indexMin,indexMax),Arrays.copyOfRange(longestTrack[2],indexMin,indexMax),dlg.fitPolynomialDegree);
		//get plot window (for closing later)
		ImagePlus impFitPlot = IJ.getImage();
		
		double[] fitCoef = cf.getParams();


		if(dlgStoreCal()){
			Prefs.set("SiMoLOc.ZC_polDegree", dlg.fitPolynomialDegree);
			switch(dlg.fitPolynomialDegree){
				case 1: Prefs.set("SiMoLOc.ZC_polCoef0", fitCoef[0]);
						Prefs.set("SiMoLOc.ZC_polCoef1", fitCoef[1]);
						Prefs.set("SiMoLOc.ZC_polCoef2", 0);
						Prefs.set("SiMoLOc.ZC_polCoef3", 0);
						Prefs.set("SiMoLOc.ZC_fitRangeMin", fitPars[0]);
						Prefs.set("SiMoLOc.ZC_fitRangeMax", fitPars[1]);
						break;
				case 2: Prefs.set("SiMoLOc.ZC_polCoef0", fitCoef[0]);
						Prefs.set("SiMoLOc.ZC_polCoef1", fitCoef[1]);
						Prefs.set("SiMoLOc.ZC_polCoef2", fitCoef[2]);
						Prefs.set("SiMoLOc.ZC_polCoef3", 0);
						Prefs.set("SiMoLOc.ZC_fitRangeMin", fitPars[0]);
						Prefs.set("SiMoLOc.ZC_fitRangeMax", fitPars[1]);
						break;
				case 3: Prefs.set("SiMoLOc.ZC_polCoef0", fitCoef[0]);
						Prefs.set("SiMoLOc.ZC_polCoef1", fitCoef[1]);
						Prefs.set("SiMoLOc.ZC_polCoef2", fitCoef[2]);
						Prefs.set("SiMoLOc.ZC_polCoef3", fitCoef[3]);
						Prefs.set("SiMoLOc.ZC_fitRangeMin", fitPars[0]);
						Prefs.set("SiMoLOc.ZC_fitRangeMax", fitPars[1]);
						break;

			}
			msg("Calibration stored!");
		}
		
		impFitPlot.close();
		
	}
	
	Plot plot;
	
	public void zCal_polCoef()
	{
		plot = new Plot("title","x-label","y-label");
		plot.addPoints(new double[]{1,2,3,4,5,6,7,8,9,10}, new double[]{1,2,3,4,5,6,7,8,9,10}, 2);
		
		final ImagePlus imp = new ImagePlus();//plot.getImagePlus();
		plot.setImagePlus(imp);
		plot.draw();
		
		final NonBlockingGenericDialog nb = new NonBlockingGenericDialog("title");
		nb.addCheckbox("clabel", true);
		//nb.addImage(imp);
		ImageCanvas ic = new ImageCanvas(imp);
		final Panel jp = new Panel();
		jp.add(ic);
		nb.add(jp);
		
		Button b = new Button("Do Fit");
		b.addActionListener(new ActionListener(){

			@Override
			public void actionPerformed(ActionEvent arg0) {
				Random rand = new Random();
				double[] y = new double[10];
				 for ( int j = 0 ; j < y.length ; j++ ) { 
		              y[j] = rand.nextInt(101);
		           }
				plot = new Plot("title","x-label","y-label");
				plot.addPoints(new double[]{1,2,3,4,5,6,7,8,9,10}, y, 2);				
				plot.setImagePlus(imp);
				plot.draw();
				for(int i=0;i<jp.getComponentCount();i++){
					if(jp.getComponent(i) instanceof ImageCanvas){
						jp.remove(i);
					}
				}
				
				ImageCanvas ic = new ImageCanvas(imp);
				
				jp.add(ic);
				
				nb.pack();
				
				///WHERE IS THIS METHOD??? EUGENE
				//nb.revalidate();
				
				
			}
			
		});
		nb.add(b);
		
		if(Toolkit.getDefaultToolkit().getAWTEventListeners(AWTEvent.MOUSE_EVENT_MASK).length == 0){
			IJ.log("no listeneres yet");
			Toolkit.getDefaultToolkit().addAWTEventListener(new AWTEventListener(){

				@Override
				public void eventDispatched(AWTEvent arg0) {
					// TODO Auto-generated method stub
					MouseEvent me = (MouseEvent) arg0;
					if(arg0.getSource() instanceof ImageCanvas){
						if(me.getID() == MouseEvent.MOUSE_CLICKED){
							
							IJ.log("x = "+String.valueOf(plot.descaleX(me.getX()))+" y = "+String.valueOf(plot.descaleY(me.getY())));
						}
					}
					
				}
			
			}
			,AWTEvent.MOUSE_EVENT_MASK);

		}
		else{
			IJ.log("listerner already in place");
		}
		
		
		nb.showDialog();
		if(Toolkit.getDefaultToolkit().getAWTEventListeners(AWTEvent.MOUSE_EVENT_MASK).length != 0){
			AWTEventListener l = Toolkit.getDefaultToolkit().getAWTEventListeners()[0];
			Toolkit.getDefaultToolkit().removeAWTEventListener(l);
			IJ.log("Listener removed");
		}
	}
	
	
	class nbdlgZCal extends NonBlockingGenericDialog{
		
		String title = "Z-calibration NonBlocking Dialog";
		
		
		public nbdlgZCal(String title) {
			super(title);
		}

		NonBlockingGenericDialog nb = new NonBlockingGenericDialog("Z-calibration NonBlocking Dialog");
		
		public void toggleFitOptions(boolean on){
			if(on){
				
			}
		}
		
		
	}
	
	public CurveFitter polyFit(double[]x, double[]y, int degree){
		
		CurveFitter cf = new CurveFitter(x, y);
		int polDegree;
		
		switch(degree){
			case 1: polDegree = ij.measure.Calibration.STRAIGHT_LINE;break;
			case 2: polDegree = ij.measure.Calibration.POLY2;break;
			case 3: polDegree = ij.measure.Calibration.POLY3;break;
			default: polDegree = ij.measure.Calibration.STRAIGHT_LINE;break;
		}
		cf.doFit(polDegree); 
		//double [] initialParams= new double[5];
		//cf.doCustomFit("y = b + c*(x-a) + d*(x-a)*(x-a)+e*(x-a)*(x-a)*(x-a)", initialParams, false);
		Fitter.plot(cf);
		
		IJ.log(cf.getResultString());
		return cf;
	}
	
	
	public int findIndex(double[] array, int value){
		for(int i=0; i<array.length; i++){
			if((int)array[i] >= value)
				return i;
		}
		return -1;
	}
	
	
	public int dlgSetZero(double zOffset) {
		GenericDialog dlg = new GenericDialog("Set new zero Z position");
		dlg.addNumericField("New zero at Z = ", (int)zOffset, 1,6," nm");
		
		
		dlg.setResizable(false);
		dlg.showDialog();
		if (dlg.wasCanceled())
            return -1;
		
		return (int)dlg.getNextNumber();
	}
	
	public int[] dlgFitRange(int nMin, int nMax) {
		GenericDialog dlg = new GenericDialog("Define Fit & Detection Z-Range");
		dlg.addNumericField("Range Z min = ", nMin, 0,5," nm");
		dlg.addNumericField("Range Z max = ", nMax, 0,5," nm");
		
		dlg.setResizable(false);
		dlg.showDialog();
		if (dlg.wasCanceled())
            return new int[]{0,0};
		
		return new int[]{(int)dlg.getNextNumber(),(int)dlg.getNextNumber()};
	}
	
	public boolean dlgStoreCal() {
		GenericDialog dlg = new GenericDialog("Store calibration?");
		dlg.addMessage("Store the new calibration?");

		dlg.setResizable(false);
		dlg.showDialog();
		return dlg.wasOKed();
	}
	/*
	double[][] transpose (double[][] array) {
	  if (array == null || array.length == 0)//empty or unset array, nothing do to here
	    return array;

	  int width = array.length;
	  int height = array[0].length;

	  double[][] array_new = new double[height][width];

	  for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
	      array_new[y][x] = array[x][y];
	    }
	  }
	  return array_new;
	}
*/
	
	void msg(String message){
		GenericDialog gd = new GenericDialog("Message");
		gd.addMessage(message);
		gd.showDialog();
	}
	void msg(double message){
		GenericDialog gd = new GenericDialog("Message");
		gd.addMessage(String.valueOf(message));
		gd.showDialog();
	}
	void msg(int message){
		GenericDialog gd = new GenericDialog("Message");
		gd.addMessage(String.valueOf(message));
		gd.showDialog();
	}

}
