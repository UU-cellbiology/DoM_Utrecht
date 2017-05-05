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
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Date;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Vector;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
import ij.io.SaveDialog;
//import ij.gui.ImageCanvas;
//import ij.gui.NonBlockingGenericDialog;
//import ij.gui.Plot;
import ij.measure.CurveFitter;
import ij.plugin.PlugIn;

public class Z_Calibration implements PlugIn{

	SMLDialog dlg = new SMLDialog();
	SMLAnalysis sml = new SMLAnalysis();

	/** dynamic array containing all tracks **/
	ArrayList<ArrayList<Double[]>> all_tracks = new ArrayList<ArrayList<Double[]>>();
	
	/** maximul length of track in frames **/
	int  maxTracklength;
	
	double []  f, sdx, sdx_err, sdy,sdy_err, trackid, particleid, tl,  sortedTracklength, sdDif, rsquared,x,y;
	/** array storing longest track used for calibration */
	double[][] longestTrack;
	/** coefficients of polynomial fit to Z curve */
	double[] fitCoefZ;
	/** coefficients of polynomial fit to X curve */
	double[] fitCoefX;
	/** coefficients of polynomial fit to Y curve */
	double[] fitCoefY;
	/** fit range of Z curve in nm*/
	int[] fitRange;
	/** zero Z offset value */
	double zOffset;
	/** X offset value */
	int zOffind;
	/** X offset value */
	//double xOffset;
	/** Y offset value */
	//double yOffset;
	/** fit range of Z curve in indexes of longestTrack array*/
	int [] indexFitRange;
	/** Main fitting dialog */
	final NonBlockingGenericDialog nb = new NonBlockingGenericDialog("Fit Z-calibration");  
	/** fitting button */
	Button b;
	
	/** original range min */
	int nRefMinZ;
	/** original range max */
	int nRefMaxZ;
	/** range min after offset*/
	int nUpdMinZ;
	/** range max after offset */
	int nUpdMaxZ;
	
	int nPolDegree=3;
	
	int nFrameMin, nFrameMax;
	
	/** whether fit was performed or not */
	boolean bFitDone = false;
	
	double [][] dDiff;
	
	Plot PlotSDdiff;
	Plot PlotSDxy;
	Plot PlotX;
	Plot PlotY;
	ImagePlus plotIp;

	
	/** main curve fitter for Z */
	CurveFitter cf;
	/** main curve fitter for X and Y */
	CurveFitter cfXY;
	
	public static ImageCanvas icDiff;
	public static ImageCanvas icSDxy;
	public static ImageCanvas icX;
	public static ImageCanvas icY;
	
	public void run(String arg) 
	{
		IJ.register(Z_Calibration.class);
		
		//save calibration
		if(arg.equals("save"))
		{
			saveZcalibration();
			return;
		}
		//save calibration
		if(arg.equals("load"))
		{
			loadZcalibration();
			return;
		}
		
		//output calibration to IJ.Log
		if(arg.equals("show"))
		{
			logZcalibration();
			return;
		}
		
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
	 * */
	
	public void zCal_particleTableInteractive()
	{
		int i,j;
		ArrayList<Double []> inner;
		
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
			IJ.run("Link Particles to Tracks", "for=[All particles] max=4 measure=[Initial position] maximum=10");
			sml = new SMLAnalysis();
			//return;
		}
		
		IJ.log("Distance between Z planes: "+Double.toString(dlg.zCalDistBetweenPlanes)+" nm");
		IJ.log("R^2 threshold: "+ Double.toString(dlg.zCalRsquareThreshold));
		if(dlg.bZCalWobbling)	
			IJ.log("Account for \"wobbling\" in X and Y: on");
		else
			IJ.log("Account for \"wobbling\" in X and Y: off");
		
		//first, sort data by particle and trackN
		Sort_Results.sorting_external_silent(sml, DOMConstants.Col_ParticleID, true);
		Sort_Results.sorting_external_silent(sml, DOMConstants.Col_TrackID, true);
		
		//load values
		get_values_from_table();
		
		
		
		if(!calculateLongestTrack())
			return;
		
		calculateAverageTracks();
		
		estimateZzeroPosition();
		zCal_calcMeanSD();
		indexFitRange = new int [2];


		PlotSDdiff = new Plot("Difference in PSF width and height (unfitted)","Z position (nm)","width-height (nm)");
		//PlotSDdiff.addPoints(dDiff[0], dDiff[1],dDiff[2], Plot.DOT);
		
		/*PlotSDdiff.addPoints(dDiff[0], dDiff[1], Plot.DOT);
		PlotSDdiff.setColor(Color.gray);
		double [] x,y,y2;
		for(i=0;i<all_tracks.size();i++)
		{
			inner = all_tracks.get(i);
			x = new double [inner.size()];
			y = new double [inner.size()];
			for(j=0;j<inner.size();j++)
			{
				x[j]=(inner.get(j)[2]*dlg.zCalDistBetweenPlanes)-zOffset;
				y[j]=inner.get(j)[8];
			}
			PlotSDdiff.addPoints(x, y, Plot.CIRCLE);
		
		}
		PlotSDdiff.setColor(Color.black);
		PlotSDdiff.addPoints(dDiff[0], dDiff[1],dDiff[2], Plot.CIRCLE);*/
		PlotSDdiff.addPoints(longestTrack[5], sdDif, Plot.CIRCLE);
		PlotSDdiff.setLegend("not fitted", Plot.AUTO_POSITION);
		//PlotSDdiff.setLimitsToFit(true);
		final ImagePlus impDiffSD = new ImagePlus();//plot.getImagePlus();
		PlotSDdiff.setImagePlus(impDiffSD);
		PlotSDdiff.draw();
		
		
		PlotSDxy = new Plot("Width and height of PSF","Z position (nm)","size of PSF (nm)");
		//PlotSDxy.addPoints(dDiff[0], dDiff[3],dDiff[4], Plot.DOT);
		//PlotSDxy.addPoints(dDiff[0], dDiff[5],dDiff[6], Plot.DOT);
		//PlotSDxy.addPoints(dDiff[0], dDiff[3], Plot.DOT);
		//PlotSDxy.addPoints(dDiff[0], dDiff[5], Plot.DOT);

		//PlotSDxy.setLimitsToFit(true);
		//PlotSDxy.setLimitsToFit(false);
		
		
		/*
		for(i=0;i<all_tracks.size();i++)
		{
			inner = all_tracks.get(i);
			x = new double [inner.size()];
			y = new double [inner.size()];
			y2 = new double [inner.size()];
			for(j=0;j<inner.size();j++)
			{
				x[j]=(inner.get(j)[2]*dlg.zCalDistBetweenPlanes)-zOffset;
				y[j]=inner.get(j)[0];
				y2[j]=inner.get(j)[1];
				
			}
			PlotSDxy.setColor(Color.red);
			PlotSDxy.addPoints(x, y, Plot.CIRCLE);
			PlotSDxy.setColor(Color.blue);
			PlotSDxy.addPoints(x, y2, Plot.CIRCLE);

		}*/
		//plot sdx
		PlotSDxy.setColor(Color.red);
		PlotSDxy.addPoints(longestTrack[2], longestTrack[0],longestTrack[3], Plot.CIRCLE);	
		//plot sdy		
		PlotSDxy.setColor(Color.blue);
		PlotSDxy.addPoints(longestTrack[2], longestTrack[1],longestTrack[4], Plot.CIRCLE);
		PlotSDxy.setColor(Color.black);
		//PlotSDxy.addPoints(dDiff[0], dDiff[3],dDiff[4], Plot.CIRCLE);
		//PlotSDxy.addPoints(dDiff[0], dDiff[5],dDiff[6], Plot.CIRCLE);
		PlotSDxy.setLegend("width	height", Plot.AUTO_POSITION);
		//PlotSDxy.setLimitsToFit(true);
		final ImagePlus impSDxy = new ImagePlus();//plot.getImagePlus();
		PlotSDxy.setImagePlus(impSDxy);
		PlotSDxy.draw();
		
		PlotX = new Plot("X wobbling","Z position (nm)","X shift (nm)");
		
		boolean bVal;
		double xyOffset;
		int nInd=0;
		/*for(i=0;i<all_tracks.size();i++)
		{
			inner = all_tracks.get(i);
			bVal = false;
			x = new double [inner.size()];
			y = new double [inner.size()];
			for(j=0;j<inner.size();j++)
			{
				x[j]=inner.get(j)[2]-zOffset;
				if(x[j]<0.5*dlg.zCalDistBetweenPlanes)
				{
					nInd =j;
					bVal=true;
				}
				y[j]=inner.get(j)[6];
			}
			if(bVal)
			{
				//finding x offset
				xyOffset = y[nInd];
				for(j=0;j<inner.size();j++)
				{
					y[j]-=xyOffset;
				}				
				
				PlotX.addPoints(x, y, Plot.CIRCLE);
			}

		}
		*/
		PlotX.addPoints(longestTrack[5], longestTrack[6], Plot.CIRCLE);
		PlotX.setLegend("not fitted", Plot.AUTO_POSITION);
		//PlotX.setLimitsToFit(true);
		final ImagePlus impX = new ImagePlus();
		PlotX.setImagePlus(impX);
		PlotX.draw();

		PlotY = new Plot("Y wobbling","Z position (nm)","Y shift (nm)");
		/*for(i=0;i<all_tracks.size();i++)
		{
			inner = all_tracks.get(i);
			bVal = false;
			x = new double [inner.size()];
			y = new double [inner.size()];
			for(j=0;j<inner.size();j++)
			{
				x[j]=inner.get(j)[2]-zOffset;
				if(x[j]<0.5*dlg.zCalDistBetweenPlanes)
				{
					nInd =j;
					bVal=true;
				}
				y[j]=inner.get(j)[7];
			}
			if(bVal)
			{
				//finding x offset
				xyOffset = y[nInd];
				for(j=0;j<inner.size();j++)
				{
					y[j]-=xyOffset;
				}				
				
				PlotY.addPoints(x, y, Plot.CIRCLE);
			}
		}
		*/
		PlotY.addPoints(longestTrack[5], longestTrack[7], Plot.CIRCLE);
		
		//PlotY.setLimitsToFit(true);
		PlotY.setLegend("not fitted", Plot.AUTO_POSITION);
		final ImagePlus impY = new ImagePlus();
		PlotY.setImagePlus(impY);
		PlotY.draw();
		
		nb.setOKLabel("Store calibration");
		//nb.setLayout(new GridLayout(2,2));
		//nb.addNumericField("test value", 10, 4);  //.addCheckbox("clabel", true);
		//nb.addImage(impDiffSD);
		icDiff = new ImageCanvas(impDiffSD);
		icSDxy = new ImageCanvas(impSDxy);
		icX = new ImageCanvas(impX);
		icY = new ImageCanvas(impY);
		
		final Panel jp = new Panel();
		
		if(dlg.bZCalWobbling)
		{	
			jp.setLayout(new GridLayout(2,2));
			jp.add(icY,0,0);
			jp.add(icSDxy,1,0);
			
			jp.add(icDiff,1,1);
			jp.add(icX,0,1);
		}
		else
		{
			jp.setLayout(new GridLayout(2,1));		
			jp.add(icDiff,0,0);
			jp.add(icSDxy,1,0);
		}
	
		b = new Button("Perform Fit");
		
		b.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent arg0) 
			{
				
				double [] dRangeZ, dRangeX, dRangeY;
				double [] fittedZ, fittedX, fittedY;
				double [] dRangeDiff;

				//Random rand = new Random();
				//double[] y = new double[10];
				int nTemp;
				
				bFitDone = true;
				b.setLabel("Fitting...");
				b.setEnabled(false);
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
				
				fittedZ=calcFittedVals(dRangeDiff,fitCoefZ);
			
				
				PlotSDdiff = new Plot("Difference in PSF width and height","Z position (nm)","width-height (nm)");
				PlotSDdiff.addPoints(dRangeZ, dRangeDiff, Plot.CIRCLE);	
				PlotSDdiff.setColor(Color.blue);
				PlotSDdiff.addPoints(fittedZ, dRangeDiff, Plot.LINE);
				PlotSDdiff.setColor(Color.black);
				PlotSDdiff.setLegend("data	fit ("+sFitResultsSplit[8]+")", Plot.AUTO_POSITION);
				PlotSDdiff.setLimitsToFit(true);
				PlotSDdiff.setImagePlus(impDiffSD);
				PlotSDdiff.draw();
				icDiff = new ImageCanvas(impDiffSD);
				
				if(dlg.bZCalWobbling)
				{
					
					
					int nOffInd = findIndexIni(longestTrack[2],(int)zOffset);
					
					double dXShift=longestTrack[6][nOffInd];
					double dYShift=longestTrack[7][nOffInd];
					for ( int j = 0 ; j < maxTracklength ; j++ ) 
					 { 
						longestTrack[6][j]=longestTrack[6][j]-dXShift;
						longestTrack[7][j]=longestTrack[7][j]-dYShift;
			         }
					
					dRangeX=Arrays.copyOfRange(longestTrack[6],indexFitRange[0],indexFitRange[1]);
					dRangeY=Arrays.copyOfRange(longestTrack[7],indexFitRange[0],indexFitRange[1]);
					
					double [] iniParam = new double[5];
					
					cfXY = new CurveFitter(dRangeZ,dRangeX);
					cfXY.doCustomFit("y = a + b*(x-e) + c*(x-e)*(x-e)+d*(x-e)*(x-e)*(x-e)",iniParam,false); 
					//cfXY.doFit(nPolDegree); 
					fitCoefX = cfXY.getParams();
					fittedX=calcFittedValsXY(dRangeZ,fitCoefX);
					sFitResults=cfXY.getResultString();
					sFitResultsSplit=sFitResults.split("[\n]+");
					PlotX = new Plot("X wobbling","Z position (nm)","X shift (nm)");
					PlotX.addPoints(dRangeZ, dRangeX, Plot.CIRCLE);
					PlotX.setColor(Color.blue);
					PlotX.addPoints(dRangeZ, fittedX, Plot.LINE);
					PlotX.setColor(Color.black);
					PlotX.setLegend("data	X fit ("+sFitResultsSplit[8]+")", Plot.AUTO_POSITION);
					PlotX.setLimitsToFit(true);
					PlotX.setImagePlus(impX);
					PlotX.draw();
					icX = new ImageCanvas(impX);
					
					
					
					cfXY = new CurveFitter(dRangeZ,dRangeY);
					//cfXY.doFit(nPolDegree); 
					
					cfXY.doCustomFit("y = a + b*(x-e) + c*(x-e)*(x-e)+d*(x-e)*(x-e)*(x-e)",iniParam,false); 
					fitCoefY = cfXY.getParams();
					fittedY=calcFittedValsXY(dRangeZ,fitCoefY);
					sFitResults=cfXY.getResultString();
					sFitResultsSplit=sFitResults.split("[\n]+");
					
					PlotY = new Plot("Y wobbling","Z position (nm)","Y shift (nm)");
					PlotY.addPoints(dRangeZ, dRangeY, Plot.CIRCLE);
					PlotY.setColor(Color.blue);
					PlotY.addPoints(dRangeZ, fittedY, Plot.LINE);
					PlotY.setColor(Color.black);
					PlotY.setLegend("data	Y fit ("+sFitResultsSplit[8]+")", Plot.AUTO_POSITION);
					PlotY.setLimitsToFit(true);
					PlotY.setImagePlus(impY);
					PlotY.draw();
					icY = new ImageCanvas(impY);
				}
				
				jp.removeAll();
				
				if(dlg.bZCalWobbling)
				{	
					jp.add(icY,0,0);
					jp.add(icSDxy,1,0);
				
					jp.add(icDiff,1,1);
					jp.add(icX,0,1);
				}
				else
				{
					jp.add(icDiff,0,0);
					jp.add(icSDxy,1,0);
				}
				b.setLabel("Perform Fit");
				b.setEnabled(true);
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
	
	/** function finds the longest track **/
	
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
		longestTrack = new double[8][maxTracklength];
		dCurrentLength = 0;
		for(i=0;i<sdx.length;i++)
		{
				if((int)trackid[i]==(int)dLongestTrackID)
				{
					if(rsquared[i]>=dlg.zCalRsquareThreshold)
					{				
						longestTrack[0][dCurrentLength]=sdx[i];
						longestTrack[1][dCurrentLength]=sdy[i];
						longestTrack[2][dCurrentLength]=(f[i]-nFrameMin)*dlg.zCalDistBetweenPlanes;
						longestTrack[3][dCurrentLength]=sdx_err[i];
						longestTrack[4][dCurrentLength]=sdy_err[i];
						/// five is reserved for shifted z axis
						longestTrack[6][dCurrentLength]=x[i];
						longestTrack[7][dCurrentLength]=y[i];
						dCurrentLength++;
					}
				}
		}
		//subtract first frame
		/*
		for(i=1;i<maxTracklength;i++)
		{
			longestTrack[2][i]=(longestTrack[2][i]-longestTrack[2][0])*dlg.zCalDistBetweenPlanes;
		}*/
		//first frame
		//longestTrack[2][0]=0;
		
		//calculate difference between SDx and SDy
		sdDif = new double[maxTracklength];

		for(i = 0; i<maxTracklength;i++)
		{
			sdDif[i] = longestTrack[0][i] - longestTrack[1][i];
			
		}
		return true;
	}
	/** function pulls together all tracks to all_tracks structure **/
	public boolean calculateAverageTracks()
	{
		int i;		
		//find min and max frame
		nFrameMin = Integer.MAX_VALUE;
		nFrameMax = Integer.MIN_VALUE;
		for(i=0;i<f.length;i++)
		{
			if(f[i]<nFrameMin)
				nFrameMin=(int)f[i];
			if(f[i]>nFrameMax)
				nFrameMax=(int)f[i];
			
			
		}
		//all tracks
		
		
	    ArrayList<Double []> inner = new ArrayList<Double []>();
	    
	    double dCurrentTrack = trackid[0];
	    Double [] dPoint;
		
		for(i=0;i<sdx.length;i++)
		{
			if((int)trackid[i]==(int)dCurrentTrack)
			{
				if(rsquared[i]>=dlg.zCalRsquareThreshold)
					{
						dPoint = new Double [9];
						dPoint[0]=sdx[i];
						dPoint[1]=sdy[i];
						//dPoint[2]=(f[i]-nFrameMin)*dlg.zCalDistBetweenPlanes;
						dPoint[2]=(f[i]-nFrameMin);
						dPoint[3]=sdx_err[i];
						dPoint[4]=sdy_err[i];
						/// five is reserved for shifted z axis
						dPoint[6]=x[i];
						dPoint[7]=y[i];
						dPoint[8]=sdx[i]-sdy[i];
						//dCurrentLength++;
						inner.add(dPoint);
					}
			}
			else
			//track is finished, add it
			{
				//arbitrary threshold to remove noise for now
				if(inner.size()>(maxTracklength-3))
					all_tracks.add(inner);
				
				inner = new ArrayList<Double []>();
				/*if(dCurrentLength>maxTracklength)
				{
					maxTracklength=dCurrentLength;
					dLongestTrackID=(int)dCurrentTrack;
				}
				dCurrentLength=1;*/
				
				dCurrentTrack=trackid[i];
				dPoint = new Double [9];
				dPoint[0]=sdx[i];
				dPoint[1]=sdy[i];
				//dPoint[2]=(f[i]-nFrameMin)*dlg.zCalDistBetweenPlanes;
				dPoint[2]=(f[i]-nFrameMin);
				dPoint[3]=sdx_err[i];
				dPoint[4]=sdy_err[i];
				/// five is reserved for shifted z axis
				dPoint[6]=x[i];
				dPoint[7]=y[i];
				dPoint[8]=sdx[i]-sdy[i];
					
				//dCurrentLength++;
				inner.add(dPoint);
				
			}
		}
		if(inner.size()>(maxTracklength-3))
			all_tracks.add(inner);
		
		
		return true;
		
	}
	
	/**function calculates average and SD values for SD,X,Y,etc**/
	
	public void zCal_calcMeanSD()
	{
		int i,j;
		int nInd;
		double [][] dZ = new double [7][nFrameMax-nFrameMin+1];
		ArrayList<Double []> inner;
		
		double [] x;
		double [] y;
		//calculate all points that are present
		for(i=0;i<all_tracks.size();i++)
		{
			inner = all_tracks.get(i);
			for(j=0;j<inner.size();j++)
			{
				nInd = inner.get(j)[2].intValue();
				dZ[0][nInd]++;
				dZ[1][nInd]+=inner.get(j)[8];//sdx-sdy
				dZ[3][nInd]+=inner.get(j)[0];//sdx
				dZ[5][nInd]+=inner.get(j)[1];//sdy
				//x[j]=inner.get(j)[2];
				///y[j]=inner.get(j)[8];
			}					
		}		
		//calculate average
		for(i=0;i<(nFrameMax-nFrameMin+1);i++)
		{
			if(dZ[0][i]>0)
			{
				dZ[1][i]=dZ[1][i]/dZ[0][i];
				dZ[3][i]=dZ[3][i]/dZ[0][i];
				dZ[5][i]=dZ[5][i]/dZ[0][i];
			}
		}
		//calculate variance
		for(i=0;i<all_tracks.size();i++)
		{
			inner = all_tracks.get(i);
			for(j=0;j<inner.size();j++)
			{
				nInd = inner.get(j)[2].intValue();
				dZ[2][nInd]+=Math.pow(inner.get(j)[8]-dZ[1][nInd],2);		
				dZ[4][nInd]+=Math.pow(inner.get(j)[0]-dZ[3][nInd],2);
				dZ[6][nInd]+=Math.pow(inner.get(j)[1]-dZ[5][nInd],2);
			}
		}
		nInd= 0;
		//calculate SD
		for(i=0;i<(nFrameMax-nFrameMin+1);i++)
		{
			if(dZ[0][i]>0)
			{
				dZ[2][i]=Math.sqrt(dZ[2][i]/dZ[0][i]);
				dZ[4][i]=Math.sqrt(dZ[4][i]/dZ[0][i]);
				dZ[6][i]=Math.sqrt(dZ[6][i]/dZ[0][i]);
				nInd++;
			}
		}		
		dDiff = new double [7][nInd];
		nInd=-1;
		for(i=0;i<(nFrameMax-nFrameMin+1);i++)
		{
			if(dZ[0][i]>0)
			{
				nInd++;
				dDiff[0][nInd] = i*dlg.zCalDistBetweenPlanes-zOffset;
				dDiff[1][nInd] =  dZ[1][i];
				dDiff[2][nInd] =  dZ[2][i];
				dDiff[3][nInd] =  dZ[3][i];
				dDiff[4][nInd] =  dZ[4][i];
				dDiff[5][nInd] =  dZ[5][i];
				dDiff[6][nInd] =  dZ[6][i];				
			}
		
		}
		
		return;
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
		Prefs.set("SiMoLOc.ZC_XYWobbling", dlg.bZCalWobbling);
		if(dlg.bZCalWobbling)
		{
			//fitted values of X polynomial
			Prefs.set("SiMoLOc.ZC_polCoefX0", fitCoefX[0]);
			Prefs.set("SiMoLOc.ZC_polCoefX1", fitCoefX[1]);
			Prefs.set("SiMoLOc.ZC_polCoefX2", fitCoefX[2]);
			Prefs.set("SiMoLOc.ZC_polCoefX3", fitCoefX[3]);
			Prefs.set("SiMoLOc.ZC_polCoefX4", fitCoefX[4]);
			//fitted values of Y polynomial
			Prefs.set("SiMoLOc.ZC_polCoefY0", fitCoefY[0]);
			Prefs.set("SiMoLOc.ZC_polCoefY1", fitCoefY[1]);
			Prefs.set("SiMoLOc.ZC_polCoefY2", fitCoefY[2]);
			Prefs.set("SiMoLOc.ZC_polCoefY3", fitCoefY[3]);
			Prefs.set("SiMoLOc.ZC_polCoefY4", fitCoefY[4]);

		}
			
		
		IJ.log("Z calibration stored!");
		//msg("Calibration stored!");
	}
	/** Saves Z calibration to file 
	 * */
	public void saveZcalibration()
	{
		SaveDialog sd = new SaveDialog("Save Z calibration", "Z_calibration", ".txt");
        String path = sd.getDirectory();
        boolean bXYWobble;
        String calDate;
        if (path==null)
        	return;
        String filename = path+sd.getFileName();
        
        
		
		//read stored values of calibration
        fitCoefZ = new double[4];
        fitRange = new int[2];
		fitCoefZ[0]=Prefs.get("SiMoLOc.ZC_polCoef0", Double.NaN);
		if(Double.isNaN(fitCoefZ[0]))
		{
			IJ.error("There is no valid Z-calibration stored, please make one!");
			return;
		}
		
		fitCoefZ[1]=Prefs.get("SiMoLOc.ZC_polCoef1", 0);
		fitCoefZ[2]=Prefs.get("SiMoLOc.ZC_polCoef2", 0);
		fitCoefZ[3]=Prefs.get("SiMoLOc.ZC_polCoef3", 0);
		fitRange[0]=(int)Prefs.get("SiMoLOc.ZC_fitRangeMin", 0);
		fitRange[1]=(int)Prefs.get("SiMoLOc.ZC_fitRangeMax", 1000);
		bXYWobble=Prefs.get("SiMoLOc.ZC_XYWobbling", false);
		if(bXYWobble)
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
		}
		
		//date
		calDate=Prefs.get("SiMoLOc.ZC_calDate","not available");
		PrintWriter writer;
		try {
			writer = new PrintWriter(filename, "UTF-8");
			writer.println("Zdate\t"+calDate);
			writer.println("XYWobble\t"+Boolean.toString(bXYWobble));
			writer.println("Zmin\t"+Integer.toString(fitRange[0]));
			writer.println("Zmax\t"+Integer.toString(fitRange[1]));
			writer.println("Zpolycoeff0\t"+Double.toString(fitCoefZ[0]));
			writer.println("Zpolycoeff1\t"+Double.toString(fitCoefZ[1]));
			writer.println("Zpolycoeff2\t"+Double.toString(fitCoefZ[2]));
			writer.println("Zpolycoeff3\t"+Double.toString(fitCoefZ[3]));
			if(bXYWobble)
			{
				writer.println("Xpolycoeff0\t"+Double.toString(fitCoefX[0]));
				writer.println("Xpolycoeff1\t"+Double.toString(fitCoefX[1]));
				writer.println("Xpolycoeff2\t"+Double.toString(fitCoefX[2]));
				writer.println("Xpolycoeff3\t"+Double.toString(fitCoefX[3]));
				writer.println("Xpolycoeff4\t"+Double.toString(fitCoefX[4]));
				
				writer.println("Ypolycoeff0\t"+Double.toString(fitCoefY[0]));
				writer.println("Ypolycoeff1\t"+Double.toString(fitCoefY[1]));
				writer.println("Ypolycoeff2\t"+Double.toString(fitCoefY[2]));
				writer.println("Ypolycoeff3\t"+Double.toString(fitCoefY[3]));
				writer.println("Ypolycoeff4\t"+Double.toString(fitCoefY[4]));
			}
			writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			IJ.log(e.getMessage());
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			//e.printStackTrace();
			IJ.log(e.getMessage());
		}

        
        IJ.log("Z calibration saved.");
		return;
	}
	
	/** Loads Z calibration from file 
	 * */
	public void loadZcalibration()
	{
		String sCalibration;
		String [] sCalSplittedRows;
		String [] sCalSplittedVals;
		String delimsn = "[\n]+";
		String delimst = "[\t]+";
		boolean bXYWobble;
		sCalibration = IJ.openAsString("");
		
		if(sCalibration==null || sCalibration.equals(""))
			return;
		
		sCalSplittedRows = sCalibration.split(delimsn);
		
		sCalSplittedVals =sCalSplittedRows[0].split(delimst);
		if(!(sCalSplittedVals[0].equals("Zdate")))
		{
			IJ.error("Not able to detect a valid Z calibration, try another file.");
			return;
		}
		Prefs.set("SiMoLOc.ZC_calDate",sCalSplittedVals[1]);
		sCalSplittedVals =sCalSplittedRows[1].split(delimst);
		bXYWobble = Boolean.parseBoolean(sCalSplittedVals[1]);
		Prefs.set("SiMoLOc.ZC_XYWobbling", bXYWobble);
		sCalSplittedVals =sCalSplittedRows[2].split(delimst);
		Prefs.set("SiMoLOc.ZC_fitRangeMin",Integer.parseInt(sCalSplittedVals[1]));
		sCalSplittedVals =sCalSplittedRows[3].split(delimst);
		Prefs.set("SiMoLOc.ZC_fitRangeMax",Integer.parseInt(sCalSplittedVals[1]));
		sCalSplittedVals =sCalSplittedRows[4].split(delimst);
		Prefs.set("SiMoLOc.ZC_polCoef0",Double.parseDouble(sCalSplittedVals[1]));
		sCalSplittedVals =sCalSplittedRows[5].split(delimst);
		Prefs.set("SiMoLOc.ZC_polCoef1",Double.parseDouble(sCalSplittedVals[1]));
		sCalSplittedVals =sCalSplittedRows[6].split(delimst);
		Prefs.set("SiMoLOc.ZC_polCoef2",Double.parseDouble(sCalSplittedVals[1]));
		sCalSplittedVals =sCalSplittedRows[7].split(delimst);
		Prefs.set("SiMoLOc.ZC_polCoef3",Double.parseDouble(sCalSplittedVals[1]));
		if(bXYWobble)
		{
			sCalSplittedVals =sCalSplittedRows[8].split(delimst);
			Prefs.set("SiMoLOc.ZC_polCoefX0",Double.parseDouble(sCalSplittedVals[1]));
			sCalSplittedVals =sCalSplittedRows[9].split(delimst);
			Prefs.set("SiMoLOc.ZC_polCoefX1",Double.parseDouble(sCalSplittedVals[1]));
			sCalSplittedVals =sCalSplittedRows[10].split(delimst);
			Prefs.set("SiMoLOc.ZC_polCoefX2",Double.parseDouble(sCalSplittedVals[1]));
			sCalSplittedVals =sCalSplittedRows[11].split(delimst);
			Prefs.set("SiMoLOc.ZC_polCoefX3",Double.parseDouble(sCalSplittedVals[1]));
			sCalSplittedVals =sCalSplittedRows[12].split(delimst);
			Prefs.set("SiMoLOc.ZC_polCoefX4",Double.parseDouble(sCalSplittedVals[1]));
			
			sCalSplittedVals =sCalSplittedRows[13].split(delimst);
			Prefs.set("SiMoLOc.ZC_polCoefY0",Double.parseDouble(sCalSplittedVals[1]));
			sCalSplittedVals =sCalSplittedRows[14].split(delimst);
			Prefs.set("SiMoLOc.ZC_polCoefY1",Double.parseDouble(sCalSplittedVals[1]));
			sCalSplittedVals =sCalSplittedRows[15].split(delimst);
			Prefs.set("SiMoLOc.ZC_polCoefY2",Double.parseDouble(sCalSplittedVals[1]));
			sCalSplittedVals =sCalSplittedRows[16].split(delimst);
			Prefs.set("SiMoLOc.ZC_polCoefY3",Double.parseDouble(sCalSplittedVals[1]));
			sCalSplittedVals =sCalSplittedRows[17].split(delimst);
			Prefs.set("SiMoLOc.ZC_polCoefY4",Double.parseDouble(sCalSplittedVals[1]));
		}
		IJ.log("Z calibration loaded.");
		return;
	}
	
	/** outputs stored Z calibration to IJ.Log window
	 * */
	public void logZcalibration()
	{
		boolean bXYWobble;
		String calDate;
		
		//read stored values of calibration
        fitCoefZ = new double[4];
        fitRange = new int[2];
		fitCoefZ[0]=Prefs.get("SiMoLOc.ZC_polCoef0", Double.NaN);
		if(Double.isNaN(fitCoefZ[0]))
		{
			IJ.error("There is no valid Z-calibration stored, please make one!");
			return;
		}
		
		fitCoefZ[1]=Prefs.get("SiMoLOc.ZC_polCoef1", 0);
		fitCoefZ[2]=Prefs.get("SiMoLOc.ZC_polCoef2", 0);
		fitCoefZ[3]=Prefs.get("SiMoLOc.ZC_polCoef3", 0);
		fitRange[0]=(int)Prefs.get("SiMoLOc.ZC_fitRangeMin", 0);
		fitRange[1]=(int)Prefs.get("SiMoLOc.ZC_fitRangeMax", 1000);
		bXYWobble=Prefs.get("SiMoLOc.ZC_XYWobbling", false);
		if(bXYWobble)
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
		}
		
		//date
		calDate=Prefs.get("SiMoLOc.ZC_calDate","not available");
		
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		IJ.log("Calibration creation date: "+ calDate);
		if(bXYWobble)	
			IJ.log("Account for \"wobbling\" in X and Y: on");
		else
			IJ.log("Account for \"wobbling\" in X and Y: off");
		
		IJ.log("Zmin: " + Integer.toString(fitRange[0]) + " nm");
		IJ.log("Zmax: " + Integer.toString(fitRange[1]) + " nm");
		IJ.log("Z polynomial coeff 0: " + Double.toString(fitCoefZ[0]));
		IJ.log("Z polynomial coeff 1: " + Double.toString(fitCoefZ[1]));
		IJ.log("Z polynomial coeff 2: " + Double.toString(fitCoefZ[2]));
		IJ.log("Z polynomial coeff 3: " + Double.toString(fitCoefZ[3]));
		if(bXYWobble)
		{
			IJ.log("X polynomial coeff 0: " + Double.toString(fitCoefX[0]));
			IJ.log("X polynomial coeff 1: " + Double.toString(fitCoefX[1]));
			IJ.log("X polynomial coeff 2: " + Double.toString(fitCoefX[2]));
			IJ.log("X polynomial coeff 3: " + Double.toString(fitCoefX[3]));
			IJ.log("X polynomial coeff 4: " + Double.toString(fitCoefX[4]));
			
			IJ.log("Y polynomial coeff 0: " + Double.toString(fitCoefY[0]));
			IJ.log("Y polynomial coeff 1: " + Double.toString(fitCoefY[1]));
			IJ.log("Y polynomial coeff 2: " + Double.toString(fitCoefY[2]));
			IJ.log("Y polynomial coeff 3: " + Double.toString(fitCoefY[3]));
			IJ.log("Y polynomial coeff 4: " + Double.toString(fitCoefY[4]));
		}
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
		x=sml.ptable.getColumnAsDoubles(DOMConstants.Col_Xnm);
		y=sml.ptable.getColumnAsDoubles(DOMConstants.Col_Ynm);

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
	double [] calcFittedVals(double [] dRangeDiff, double fitCoef[])
	{
		double [] fittedZ = new double [dRangeDiff.length];
		int k, m;
		double dVal;
		for(k=0;k<dRangeDiff.length;k++)
		{
			dVal =fitCoef[0];
			for (m=1;m<=nPolDegree;m++)
				dVal+=fitCoef[m]*Math.pow(dRangeDiff[k], m);
			fittedZ[k] = dVal;
		}
		return fittedZ;
	}
	
	double [] calcFittedValsXY(double [] dRangeDiff, double fitCoef[])
	{
		double [] fittedZ = new double [dRangeDiff.length];
		int k, m;
		double dVal;
		for(k=0;k<dRangeDiff.length;k++)
		{
			dVal =fitCoef[0];
			for (m=1;m<=3;m++)
				dVal+=fitCoef[m]*Math.pow(dRangeDiff[k]-fitCoef[4], m);
			fittedZ[k] = dVal;
		}
		return fittedZ;
	}
	
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
		//int zOffind;
		double xOffset;
		double yOffset;
		
		zOffset = 0;
		zOffind=nBeg;
		double dMin = Double.MAX_VALUE;
		for (i=nBeg;i<=nEnd;i++)
		{
			
			if(Math.abs(sdDif[i])<dMin)
			{
				zOffset = longestTrack[2][i];
				dMin = Math.abs(sdDif[i]);
				zOffind = i;
			}
			
		}
		fitRange = new int[2];
		nRefMinZ=(int)longestTrack[2][0];
		nRefMaxZ=(int)longestTrack[2][maxTracklength-1];
		fitRange[0] = nRefMinZ;
		fitRange[1] = nRefMaxZ;
		xOffset = longestTrack[6][zOffind];
		yOffset = longestTrack[7][zOffind];
		//subtract z-offset from all z-values and store it separately
		for(i = 0; i < maxTracklength; i++)
		{
			longestTrack[5][i] = longestTrack[2][i] - zOffset;  
			//correct X and Y
			longestTrack[6][i] = longestTrack[6][i] - xOffset;  
			longestTrack[7][i] = longestTrack[7][i] - yOffset;  
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
