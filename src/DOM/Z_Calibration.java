package DOM;
//in detect molecules on line 329
//add width and height of original image to the results table
//sml.ptable.addValue("Original_image_size",imp.getWidth());
//sml.ptable.addValue("Original_image_size",imp.getHeight());


import java.awt.AWTEvent;
import java.awt.Button;
import java.awt.Choice;
import java.awt.Color;
import java.awt.GridLayout;
import java.awt.Panel;
import java.awt.TextField;
import java.awt.Toolkit;
import java.awt.event.AWTEventListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.Date;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Vector;

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
	/** array storing longest track used for calibration */
	double[][] longestTrack;
	/** coefficients of polynomial fittes to z curve */
	double[] fitCoefZ;
	/** fit range of Z curve in nm*/
	int[] fitRange;
	/** zero Z offset value */
	double zOffset;
	/** fit range of Z curve in indexes of longestTrack array*/
	int [] indexFitRange;
	/** Main fitting dialog */
	final NonBlockingGenericDialog nb = new NonBlockingGenericDialog("Fit Z-calibration");  
	
	/** original range min */
	int nRefMinZ;
	/** original range max */
	int nRefMaxZ;
	/** range min after offset*/
	int nUpdMinZ;
	/** range max after offset */
	int nUpdMaxZ;
	
	int nPolDegree=3;
	
	/** whether fit was performed or not */
	boolean bFitDone = false;
	
	Plot PlotSDdiff;
	Plot PlotSDxy;
	ImagePlus plotIp;
	
	/** main curve fitter */
	CurveFitter cf;
	
	public static ImageCanvas icDiff;
	public static ImageCanvas icSDxy;
	
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
				zCal_particleTableInteractive();
				//zCal_imageStack(); 
				break;
			case 1: 
				IJ.log("Creating Z calibration from Results table");
				zCal_particleTableInteractive(); 				
				break;
			//case 2:
			//	IJ.log("Creating Z calibration from user provided coefficients");
			//	zCal_polCoef(); 
			//	break;
		}
	}
		
	
	/** making z-calibration from Particles table,
	 * link particles if data is not provided
	 * 
	 * OLD VERSION, DEPRECATED
	 * */
	
	/*
	
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
		
		get_values_from_table();		
		
		
		if(!calculateLongestTrack())
			return;
		
		
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
		fitRange = dlgFitRange((int)longestTrack[2][0],(int)longestTrack[2][maxTracklength-1]);
		ip.close();
		if(fitRange[0] == 0 & fitRange[1] == 0)
			return;
		
		indexFitRange = new int [2];
		indexFitRange[0]=findIndex(longestTrack[2],fitRange[0]);
		indexFitRange[1]=findIndex(longestTrack[2],fitRange[1]);

		
		
		CurveFitter cf = polyFit(Arrays.copyOfRange(sdDif,indexFitRange[0],indexFitRange[1]),Arrays.copyOfRange(longestTrack[2],indexFitRange[0],indexFitRange[1]),dlg.fitPolynomialDegree);
		//get plot window (for closing later)
		ImagePlus impFitPlot = IJ.getImage();
		
		fitCoefZ = cf.getParams();


		if(dlgStoreCal())
		{
			zCal_storeCalibration();
		}
		
		impFitPlot.close();
		
	}
	
*/
	
	public void zCal_particleTableInteractive()
	{
		
		
		//check that the table is present
		if (sml.ptable.getCounter()==0 || !sml.ptable.getHeadings()[0].equals("X_(px)"))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for reconstruction, please load one.");
			return;
		}
		
		//check that the particles are linked
		if ( !sml.ptable.columnExists(DOMConstants.Col_TrackID))
		{
			//IJ.error("Particle table is probably not linked.");
			IJ.log("WARNING! Particles are not linked in the table, performing automatic linking.");
			IJ.run("Link Particles to Tracks", "for=[All particles] max=4 measure=[Initial position] maximum=10");;
			sml = new SMLAnalysis();
			//return;
		}
		
		IJ.log("Distance between Z planes: "+Double.toString(dlg.zCalDistBetweenPlanes)+" nm");
		IJ.log("R^2 threshold: "+ Double.toString(dlg.zCalRsquareThreshold));
		
		//first, sort data by particle and trackN
		Sort_Results.sorting_external_silent(sml, DOMConstants.Col_ParticleID, true);
		Sort_Results.sorting_external_silent(sml, DOMConstants.Col_TrackID, true);
		
		//load values
		get_values_from_table();
		if(!calculateLongestTrack())
			return;
		estimateZzeroPosition();
		indexFitRange = new int [2];
		/*
		fitRange[0] = (int)longestTrack[2][0];
		fitRange[1] = (int)longestTrack[2][maxTracklength-1]);
		*/
		PlotSDdiff = new Plot("Difference in PSF width and height (unfitted)","Z position (nm)","width-height (nm)");
		PlotSDdiff.addPoints(longestTrack[5], sdDif, Plot.CONNECTED_CIRCLES);
		//PlotSDxy.addPoints(new double[]{1,2,3,4,5,6,7,8,9,10}, new double[]{1,2,3,4,5,6,7,8,9,10}, 2);
		PlotSDdiff.setLegend("not fitted", Plot.AUTO_POSITION);
		final ImagePlus impDiffSD = new ImagePlus();//plot.getImagePlus();
		PlotSDdiff.setImagePlus(impDiffSD);
		PlotSDdiff.draw();
		
		PlotSDxy = new Plot("Width and height of PSF","Z position (nm)","size of PSF (nm)");
		//plot sdx
		PlotSDxy.setColor(Color.red);
		//p.addPoints(longestTrack[2], longestTrack[0], Plot.CONNECTED_CIRCLES);	
		PlotSDxy.addPoints(longestTrack[2], longestTrack[0],longestTrack[3], Plot.CIRCLE);	
		//plot sdy
		PlotSDxy.setColor(Color.blue);
		PlotSDxy.addPoints(longestTrack[2], longestTrack[1],longestTrack[4], Plot.CIRCLE);
		PlotSDxy.setColor(Color.black);

		PlotSDxy.setLegend("width	height", Plot.AUTO_POSITION);
		
		//PlotSDdiff.addPoints(longestTrack[2], sdDif, Plot.CONNECTED_CIRCLES);
		//PlotSDxy.addPoints(new double[]{1,2,3,4,5,6,7,8,9,10}, new double[]{1,2,3,4,5,6,7,8,9,10}, 2);
		
		final ImagePlus impSDxy = new ImagePlus();//plot.getImagePlus();
		PlotSDxy.setImagePlus(impSDxy);
		PlotSDxy.draw();
		
		//final NonBlockingGenericDialog 
		//nb = new NonBlockingGenericDialog("Fit Z-calibration");
		nb.setOKLabel("Save calibration");
		//nb.setLayout(new GridLayout(2,2));
		//nb.addNumericField("test value", 10, 4);  //.addCheckbox("clabel", true);
		//nb.addImage(impDiffSD);
		icDiff = new ImageCanvas(impDiffSD);
		icSDxy = new ImageCanvas(impSDxy);
		final Panel jp = new Panel();
		//final Panel jpParams = new Panel();
		
		jp.setLayout(new GridLayout(2,1));
		//jp.add(icSDxy,0,0);
		//jp.add(icDiff,1,0);
	
		jp.add(icDiff,0,0);
		jp.add(icSDxy,1,0);
	
		//nb.add(jp);
		
		
		Button b = new Button("Perform Fit");
		b.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent arg0) 
			{
				
				double [] dRangeZ;
				double [] fittedZ;
				double [] dRangeDiff;
				double dVal;
				//Random rand = new Random();
				//double[] y = new double[10];
				int nTemp,k,m;
				
				bFitDone = true;
				//reading values
				Vector parameters;
				TextField tf1, tf2;
				Choice wPol;
				parameters=nb.getChoices();
				wPol = (Choice) parameters.get(0);
				nPolDegree=wPol.getSelectedIndex()+1;
				parameters= nb.getNumericFields();				
				
				//Z-offset
				tf1 = (TextField) parameters.get(0);
				zOffset = borderCheckOffset((int)Double.parseDouble(tf1.getText()));
				tf1.setText(Integer.toString((int)zOffset));
				//update borders
				nUpdMinZ = (int) (nRefMinZ-zOffset);
				nUpdMaxZ = (int) (nRefMaxZ-zOffset);
				//left border
				tf1 = (TextField) parameters.get(1);
				fitRange[0]=borderCheckRange((int)Double.parseDouble(tf1.getText()));
				tf2 = (TextField) parameters.get(2);
				fitRange[1]=borderCheckRange((int)Double.parseDouble(tf2.getText()));
				if (fitRange[0]>fitRange[1])
				{
					nTemp=fitRange[0];
					fitRange[0] = fitRange[1];
					fitRange[1] = nTemp;
				}
				tf1.setText(Integer.toString(fitRange[0]));
				tf2.setText(Integer.toString(fitRange[1]));
				indexFitRange[0]=findIndexIni(longestTrack[5],fitRange[0]);
				indexFitRange[1]=findIndexLast(longestTrack[5],fitRange[1]);
				
				for ( int j = 0 ; j < maxTracklength ; j++ ) 
				 { 
		             // y[j] = rand.nextInt((int)Double.parseDouble(tf1.getText()));
					longestTrack[5][j]=longestTrack[2][j]-zOffset;
		         }
				
				dRangeZ=Arrays.copyOfRange(longestTrack[5],indexFitRange[0],indexFitRange[1]);
				dRangeDiff=Arrays.copyOfRange(sdDif,indexFitRange[0],indexFitRange[1]);
				//do the fitting
				int polDegree;
				cf = new CurveFitter(dRangeDiff,dRangeZ);
				switch(nPolDegree)
				{
					case 1: polDegree = ij.measure.Calibration.STRAIGHT_LINE;break;
					case 2: polDegree = ij.measure.Calibration.POLY2;break;
					case 3: polDegree = ij.measure.Calibration.POLY3;break;
					default: polDegree = ij.measure.Calibration.STRAIGHT_LINE;break;
				}
				cf.doFit(polDegree); 
				//String delimsn = "[\n]+";
				String [] sFitResultsSplit;
				String sFitResults=cf.getResultString();
				sFitResultsSplit=sFitResults.split("[\n]+");
				fitCoefZ = cf.getParams();
				fittedZ=new double[dRangeZ.length];
				//int teeest = dRangeZ.length;
				//calculate fitted approximations
				for(k=0;k<dRangeZ.length;k++)
				{
					dVal =fitCoefZ[0];
					for (m=1;m<=nPolDegree;m++)
						dVal+=fitCoefZ[m]*Math.pow(dRangeDiff[k], m);
					fittedZ[k] = dVal;
				}
				
				PlotSDdiff = new Plot("Difference in PSF width and height","Z position (nm)","width-height (nm)");
				//PlotSDdiff.setColor(Color.red);
				PlotSDdiff.addPoints(dRangeZ, dRangeDiff, Plot.CIRCLE);	
				PlotSDdiff.setColor(Color.blue);
				PlotSDdiff.addPoints(fittedZ, dRangeDiff, Plot.LINE);
				PlotSDdiff.setColor(Color.black);
				PlotSDdiff.setLegend("data	fit ("+sFitResultsSplit[8]+")", Plot.AUTO_POSITION);
				PlotSDdiff.setLimitsToFit(true);
				PlotSDdiff.setImagePlus(impDiffSD);
				PlotSDdiff.draw();
				
				int nPlot=0;
				for(int i=0;i<jp.getComponentCount();i++)
				{
					if(jp.getComponent(i) instanceof ImageCanvas)
					{
						nPlot++;
						if(nPlot==2)
							jp.remove(i);
					}
				}
				
				icDiff = new ImageCanvas(impDiffSD);
				
				//jp.add(icDiff,0,0);
				
				jp.add(icDiff,1);
				//jp.add(icSDxy,1,0);
				
				nb.pack();
				
				nb.revalidate();
				
				
			}
			
		});
		//jpParams.add(comp)
		String [] zcPolDegreeOptions = new String [] {"1","2","3"};
		nb.addChoice("Polynomial degree", zcPolDegreeOptions, "3");
		nb.addNumericField("Z zero position: ", zOffset, 0, 4," nm    <- top plot ");  //.addCheckbox("clabel", true);
		nb.addNumericField("Range Z min = ", nUpdMinZ, 0,4," nm    <- bottom plot");  //.addCheckbox("clabel", true);
		nb.addNumericField("Range Z max = ", nUpdMaxZ, 0,4," nm    <- bottom plot");  
		//Label RFit= new Label("Goodness of fit (R^2): undefined");
		//nb.add(RFit,2);
		nb.add(b,1);
		
		
		nb.addPanel(jp);
		
		//I have no idea what this part is doing (EUGENE)
		if(Toolkit.getDefaultToolkit().getAWTEventListeners(AWTEvent.MOUSE_EVENT_MASK).length == 0)
		{
			//IJ.log("no listeneres yet");
			Toolkit.getDefaultToolkit().addAWTEventListener(new AWTEventListener(){

				@Override
				public void eventDispatched(AWTEvent arg0) {
					// TODO Auto-generated method stub
					MouseEvent me = (MouseEvent) arg0;
					if(arg0.getSource() instanceof ImageCanvas){
						if(me.getID() == MouseEvent.MOUSE_CLICKED){
							
							IJ.log("x = "+String.valueOf(PlotSDdiff.descaleX(me.getX()))+" y = "+String.valueOf(PlotSDdiff.descaleY(me.getY())));
						}
					}
					
				}
			
			}
			,AWTEvent.MOUSE_EVENT_MASK);

		}
		else
		{
			//IJ.log("listerner already in place");
		}
		
		
		nb.showDialog();
		if(nb.wasOKed() && bFitDone)
		{
			//print out full fitting information
			IJ.log(cf.getResultString());
			zCal_storeCalibration();
		}
		else
		{
			IJ.log("Z calibration WAS NOT stored!");
		}
		if(Toolkit.getDefaultToolkit().getAWTEventListeners(AWTEvent.MOUSE_EVENT_MASK).length != 0){
			AWTEventListener l = Toolkit.getDefaultToolkit().getAWTEventListeners()[0];
			Toolkit.getDefaultToolkit().removeAWTEventListener(l);
			//.IJ.log("Listener removed");
		}
	}
	
	/*
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
	}*/
	
	
	public int findIndexIni(double[] array, int value){
		for(int i=0; i<array.length; i++){
			if((int)array[i] >= value)
				return i;
		}
		return -1;
	}
	public int findIndexLast(double[] array, int value){
		for(int i=0; i<array.length; i++)
		{
			if((int)array[i] == value)
				return i+1;
			if((int)array[i] > value)
				return i;
		}
		return -1;
	}
	
	/*
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
	*/
	
	public boolean calculateLongestTrack()
	{
		
		int i;
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
		//case only one track is present
		
		if(dCurrentTrack == trackid[0])
			maxTracklength=dCurrentLength;
		//check that track has enough point
		if(maxTracklength<10)
		{
			IJ.error("Error! Longest track is less than 10 points. Adjust linking, decrease z-step or choose different calibration stack.");
			return false;			
		}
		IJ.log("Longest track length after filtering: "+Integer.toString(maxTracklength) + " particles");
		
		//fill in data from the longest track
		longestTrack = new double[6][maxTracklength];
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
		return true;
	}
	
	public void zCal_storeCalibration()
	{
		DateFormat df = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss");

		// Get the date today using Calendar object.
		Date today = Calendar.getInstance().getTime();        
		// Using DateFormat format method we can create a string 
		// representation of a date with the defined format.
		String reportDate = df.format(today);
		
		Prefs.set("SiMoLOc.ZC_polDegree", nPolDegree);
		Prefs.set("SiMoLOc.ZC_calDate",reportDate);
		switch(nPolDegree)
		{
			case 1: Prefs.set("SiMoLOc.ZC_polCoef0", fitCoefZ[0]);
					Prefs.set("SiMoLOc.ZC_polCoef1", fitCoefZ[1]);
					Prefs.set("SiMoLOc.ZC_polCoef2", 0);
					Prefs.set("SiMoLOc.ZC_polCoef3", 0);
					Prefs.set("SiMoLOc.ZC_fitRangeMin", fitRange[0]);
					Prefs.set("SiMoLOc.ZC_fitRangeMax", fitRange[1]);
					break;
			case 2: Prefs.set("SiMoLOc.ZC_polCoef0", fitCoefZ[0]);
					Prefs.set("SiMoLOc.ZC_polCoef1", fitCoefZ[1]);
					Prefs.set("SiMoLOc.ZC_polCoef2", fitCoefZ[2]);
					Prefs.set("SiMoLOc.ZC_polCoef3", 0);
					Prefs.set("SiMoLOc.ZC_fitRangeMin", fitRange[0]);
					Prefs.set("SiMoLOc.ZC_fitRangeMax", fitRange[1]);
					break;
			case 3: Prefs.set("SiMoLOc.ZC_polCoef0", fitCoefZ[0]);
					Prefs.set("SiMoLOc.ZC_polCoef1", fitCoefZ[1]);
					Prefs.set("SiMoLOc.ZC_polCoef2", fitCoefZ[2]);
					Prefs.set("SiMoLOc.ZC_polCoef3", fitCoefZ[3]);
					Prefs.set("SiMoLOc.ZC_fitRangeMin", fitRange[0]);
					Prefs.set("SiMoLOc.ZC_fitRangeMax", fitRange[1]);
					break;

		}
		IJ.log("Z calibration stored!");
		//msg("Calibration stored!");
	}
	public void get_values_from_table()
	{
		//get values
		sdx = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_X);
		sdx_err = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_X_err);
		sdy = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_Y);
		sdy_err = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SD_Y_err);
		f = sml.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);
		trackid = sml.ptable.getColumnAsDoubles(DOMConstants.Col_TrackID);
		particleid = sml.ptable.getColumnAsDoubles(DOMConstants.Col_ParticleID);
		rsquared = sml.ptable.getColumnAsDoubles(DOMConstants.Col_chi);

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
	/** check border position for offset */
	public int borderCheckOffset(int inParam)
	{
		if (inParam<nRefMinZ)
			return nRefMinZ;
		if(inParam>nRefMaxZ)
			return nRefMaxZ;
		
		return inParam;
	}
	/** check border position for range */
	public int borderCheckRange(int inParam)
	{
		if (inParam<nUpdMinZ)
			return nUpdMinZ;
		if(inParam>nUpdMaxZ)
			return nUpdMaxZ;
		
		return inParam;
	}

	/** function estimates position of focal plane
	 * by minimum of SDx-SDy and total range of calibration */
	
	void estimateZzeroPosition()
	{
		//Estimate position of zero
		//let's take middle third
		int nBeg = (int)Math.round(maxTracklength/3.0);
		int nEnd = (int)Math.round(2.0*maxTracklength/3.0);
		int i;
		zOffset = 0;
		double dMin = Double.MAX_VALUE;
		for (i=nBeg;i<=nEnd;i++)
		{
			
			if(Math.abs(sdDif[i])<dMin)
			{
				zOffset = longestTrack[2][i];
				dMin = Math.abs(sdDif[i]);
			}
			
		}
		fitRange = new int[2];
		nRefMinZ=(int)longestTrack[2][0];
		nRefMaxZ=(int)longestTrack[2][maxTracklength-1];
		fitRange[0] = nRefMinZ;
		fitRange[1] = nRefMaxZ;
		
		//subtract z-offset from all z-values and store it separately
		for(i = 0; i < maxTracklength; i++)
		{
			longestTrack[5][i] = longestTrack[2][i] - zOffset;  
		}
		nUpdMinZ= (int) (nRefMinZ-zOffset);
		nUpdMaxZ= (int) (nRefMaxZ-zOffset);
	}
	
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
