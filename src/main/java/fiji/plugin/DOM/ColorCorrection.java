package fiji.plugin.DOM;


import java.awt.Color;
import java.awt.Rectangle;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.Arrow;
import ij.gui.Line;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.io.SaveDialog;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;


public class ColorCorrection implements PlugIn
{
	
	ImagePlus imp;
	ImagePlus imp_refMax;
	ImagePlus imp_warpMax;
	
	//ImagePlus mark_ref_imp;
	//ImagePlus mark_warp_imp;

	
	ImageStack is_temp;

	SMLAnalysis sml;
	SMLDialog dlg = new SMLDialog();
	
	double [][] xyref; 
	double [][] xystat;
	double [][] xymov;
	double [] SNRref, fref;
	double [][] xywarp; 
	double [] SNRwarp, fwarp;
	/** Shift between channels in pixels */
	double [] dXYShift;
	/** dimension of the B-spline grid*/
	int dimX, dimY;
	/** B-Spline grid coefficients */
	double [][][] O_trans;
	/** list of colocalized particles coordinates together */
	double [] xy_col;
	/** Max SNR threshold*/
	double dCCSNR;
	/** Max distance between particles in pixels */
	double dCCDist;
	/**number of colocalized particles */
	int nColocPatN;
	/** image width and height*/
	double dMheight, dMwidth;
	/** pixel size in nm*/
	double dPxtoNm;
	
	int [] Spacing;
	/** Main fitting dialog */
	final NonBlockingGenericDialog nbStore = new NonBlockingGenericDialog("Save/Store color calibration");  
	
	/** date when calibration was created*/	
	String reportDate;
	
	
	@Override
	public void run(String arg) 
	{
	
		int i;
		
		Roi refROI, warpedROI;
		boolean bRefRoi = false;
		boolean bWarpRoi = false;
		Rectangle nRect;
		
		
		IJ.register(ColorCorrection.class);
		
		
		//output calibration to IJ.Log
		if(arg.equals("show"))
		{
			ColocCal_showCalibration();
			return;
		}
		
		//correct Results table
		if(arg.equals("table"))
		{
			if(ColocCal_getCalibration())
			{
				sml = new SMLAnalysis();
				ColocCal_transformTable();
			}
			return;
		}
		
		//correct Results table
		if(arg.equals("image"))
		{
			if(ColocCal_getCalibration())
			{
				sml = new SMLAnalysis();
				ColocCal_transformImage();
			}
			return;
		}
		
		//save calibration
		if(arg.equals("save"))
		{
			if(ColocCal_getCalibration())
			{
				ColocCal_Save();
			}
			return;
			
		}
		//save calibration
		if(arg.equals("load"))
		{

			ColocCal_Load();
			return;
		}
		
		//show parameters dialog
		if(!dlg.dglColorCalibration())
			return;
		dCCSNR= dlg.dCCSNR;
		dCCDist = dlg.dCCDist;
		
		NonBlockingGenericDialog nb = new NonBlockingGenericDialog("Choose image of REFERENCE channel");
		
		nb.addMessage("Choose calibration image/stack (or ROI) for REFERENCE channel");
		nb.showDialog();
		if(!nb.wasOKed())
			return;
		//let's start detection
		
		IJ.log("Creating color chromatic calibration from image stack");
		IJ.log("Step 1: Detecting molecules in reference channel");
		IJ.run("Detect Molecules");
		
		imp = IJ.getImage();
		
		if (imp==null)
		{
		    IJ.noImage();
		    return;
		}
		else if (imp.getType() != ImagePlus.GRAY8 && imp.getType() != ImagePlus.GRAY16 ) 
		{
		    IJ.error("8 or 16 bit greyscale image required");
		    return;
		}
			
		refROI=imp.getRoi();
		
		if(!(refROI==null))
		{
			bRefRoi = true;
			imp=imp.duplicate();
		}
		
				
		ZProjector zProj = new ZProjector(imp);
		zProj.setMethod(ZProjector.MAX_METHOD);
	
		zProj.doProjection();
		imp_refMax = zProj.getProjection();
		
		//imp_refMax.show();
		
		sml = new SMLAnalysis();
		//get detection results
		SNRref = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SNR);
		fref = sml.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);
		xyref = new double [2][fref.length];
		xyref[0] = sml.ptable.getColumnAsDoubles(DOMConstants.Col_X);
		xyref[1] = sml.ptable.getColumnAsDoubles(DOMConstants.Col_Y);
		
		if(bRefRoi)
		{
			//if ROI is selected, let's move point of origin of coordinates
			nRect = refROI.getBounds();
			for(i=0;i<xyref[0].length;i++)
			{
				xyref[0][i]-=nRect.getX();
				xyref[1][i]-=nRect.getY();
			}
		}

		dPxtoNm =  sml.ptable.getValueAsDouble(DOMConstants.Col_Xnm, 0)/sml.ptable.getValueAsDouble(DOMConstants.Col_X, 0);
		
		
		nb = new NonBlockingGenericDialog("Choose image of WARPED channel");
		nb.addMessage("Choose calibration image/stack (or ROI) for WARPED channel");
		nb.showDialog();
		if(!nb.wasOKed())
			return;

		
		imp = IJ.getImage();
		if (imp==null)
		{
		    IJ.noImage();
		    return;
		}
		else if (imp.getType() != ImagePlus.GRAY8 && imp.getType() != ImagePlus.GRAY16 ) 
		{
		    IJ.error("8 or 16 bit greyscale image required");
		    return;
		}
		//let's start detection	
		IJ.log("Step 2: Detecting molecules in warped channel");
		IJ.run("Detect Molecules");
		
		warpedROI=imp.getRoi();
		if(!(warpedROI==null))
		{
			bWarpRoi = true;
			imp=imp.duplicate();
		}
		zProj = new ZProjector(imp);
		zProj.setMethod(ZProjector.MAX_METHOD);
		zProj.doProjection();
		imp_warpMax = zProj.getProjection();
		
		//imp_warpMax.show();
		sml = new SMLAnalysis();
		//get detection results
		SNRwarp = sml.ptable.getColumnAsDoubles(DOMConstants.Col_SNR);
		fwarp = sml.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);
		xywarp = new double [2][fwarp.length];
		xywarp[0] = sml.ptable.getColumnAsDoubles(DOMConstants.Col_X);
		xywarp[1] = sml.ptable.getColumnAsDoubles(DOMConstants.Col_Y);
		
		
		if(bWarpRoi)
		{
			//if ROI is selected, let's move point of origin of coordinates
			nRect = warpedROI.getBounds();
			for(i=0;i<xyref[0].length;i++)
			{
				xywarp[0][i]-=nRect.getX();
				xywarp[1][i]-=nRect.getY();
			}
		}
		
		IJ.log("Aligning maximum projection pictures...");
		
		AdjustDimensions();
		ImCrossCorrelation imCrCorr = new ImCrossCorrelation();
		//dXYShift = imCrCorr.calcShiftFFTCorrelation(imp_refMax.getProcessor(),imp_warpMax.getProcessor());
		dXYShift = imCrCorr.calcShiftFFTCorrelationDouble(imp_refMax.getProcessor(),imp_warpMax.getProcessor(),1);
				
		IJ.showProgress(1.0);
		
		IJ.log("X shift: "+Double.toString(dXYShift[0])+ " px");
		IJ.log("Y shift: "+Double.toString(dXYShift[1]) + " px");
		IJ.log("Step 3: Finding corresponding images of particles in both channels...");
		
		//update warped particles coordinates
		for(i=0;i<fwarp.length;i++)
		{
			
			xywarp[0][i]+=dXYShift[0];
			xywarp[1][i]+=dXYShift[1];

		}	
		
		//find images of particles in both channels
		if(CCfindParticles()<3)
		{
			IJ.error("Less than 3 particles are found in both channels. Try to increase distance or lower SNR threshold.");
		    return;
		}
				
				
		IJ.log("In total "+Integer.toString(nColocPatN)+" particles found.");
		IJ.log("Density: " + Double.toString(nColocPatN*1000000/(dMheight*dMwidth*dPxtoNm*dPxtoNm))+" particles/um2");
		
		
		IJ.log("Step 4: Calculating B-spline transform...");
		O_trans = calculateBsplineTransform(5.0);
		IJ.log("..done.");
		
		if(dlg.bCCShowPoints)
		{
			ShowRefPoints();
		}
		if(dlg.bCCShowMap)
		{
			GenerateMap();
		}
		if(dlg.bCCShowResults)
		{
			OutputResults();
		}
		
		DateFormat df = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss");
		// Get the date today using Calendar object.
		Date today = Calendar.getInstance().getTime();        
		// Using DateFormat format method we can create a string 
		// representation of a date with the defined format.
		reportDate = df.format(today);
		
		String [] ccSaveOptions = new String [] {"Store and save calibration","Only save calibration","Only store calibration"};
		nbStore.addChoice("I'd like to:", ccSaveOptions, "Store and save calibration");
		nbStore.showDialog();
		if(nbStore.wasOKed())
		{
			int nChoice = nbStore.getNextChoiceIndex(); 
			if(nChoice==0 || nChoice==2)
			{
				ColocCal_storeCalibration();
			};
			if(nChoice==0 || nChoice==1)
			{
				ColocCal_Save();
			};
						
			
		}
		//nbStore.setOKLabel("Store calibration");
		//Button b = new Button("Store and save calibration");
		
		return;
	}
	
	/** function calculating b-spline transform
	 * using provided coordinates
	 * based on matlab code of
	 *  D.Kroon University of Twente (August 2010)
	 * */
	double [][][] calculateBsplineTransform(double MaxRef)
	{
		int MaxItt;
		Spacing = new int [2];
		int dx,dy;
		
		double [][][] O_ref;
		double [][][] O_add;
		double [][][] O_int;
		double [][] R;
		double [][] xyreg;
		double dDistMean;
		int i,j,k, nIterCount;
		
		 // Calculate max refinements steps
	     MaxItt=Math.min((int)Math.floor(Math.log(dMheight*0.5)/Math.log(2.0)), (int)Math.floor(Math.log(dMwidth*0.5)/Math.log(2.0))); 
	     // set b-spline grid spacing in x and y direction
	     
	     Spacing[0]=(int) Math.pow(2,MaxItt);
	     Spacing[1]=(int) Math.pow(2,MaxItt);
	     // Make an initial uniform b-spline grid
		
	     // Determine grid spacing
	     dx=Spacing[0]; 
	     dy=Spacing[1];
	     dimX=(int) ((dMwidth+(dx*3))/dx)+1;
	     dimY=(int) ((dMheight+(dy*3))/dy)+1;
	     
	     
	     //[X,Y]=ndgrid(-dx:dx:(sizeI(1)+(dx*2)),-dy:dy:(sizeI(2)+(dy*2)));
	     O_ref =new double [2][dimX][dimY];
	     O_add =new double [2][dimX][dimY];
	     O_int  =new double [2][dimX][dimY];
	     for(i=0;i<dimX;i++)
	    	 for(j=0;j<dimY;j++)
	    	 {
	    		 O_ref[0][i][j]=(i-1)*dx;
	    		 O_ref[1][i][j]=(j-1)*dy;
	    	 }
	     
	     R= new double[2][nColocPatN];
	     dDistMean=0;
	     for(i=0;i<nColocPatN;i++)
	     {
	    	 R[0][i]=xystat[0][i]-xymov[0][i];
	    	 R[1][i]=xystat[1][i]-xymov[1][i];
	    	 dDistMean+=Math.sqrt(R[0][i]*R[0][i]+R[1][i]*R[1][i]);
	     }
	     dDistMean/=nColocPatN;
	     //display new average
	     IJ.log("Initial average distance (minus translation): " + Double.toString(dDistMean*dPxtoNm) + " nm");

	     for(nIterCount = 0;nIterCount<(int)MaxRef;nIterCount++)
	     {
		     // Make a b-spline grid which minimizes the difference between the
		     // corresponding points
		     O_add= bspline_grid_fitting_2d(O_add, Spacing, R, xymov);
		     
		     O_int  =new double [2][dimX][dimY];
		     for(k=0;k<2;k++)
			     for(i=0;i<dimX;i++)
			    	 for(j=0;j<dimY;j++)
			    	 {
			    		 O_int[k][i][j]=O_ref[k][i][j]+O_add[k][i][j];
			    	 }
		     xyreg = bspline_transform_slow_2d(O_int, Spacing, xymov);
		     
		     // Calculate the remaining difference between the points
		     R= new double[2][nColocPatN];
		     dDistMean=0;
		     for(i=0;i<nColocPatN;i++)
		     {
		    	 R[0][i]=xystat[0][i]-xyreg[0][i];
		    	 R[1][i]=xystat[1][i]-xyreg[1][i];
		    	 dDistMean+=Math.sqrt(R[0][i]*R[0][i]+R[1][i]*R[1][i]);
		     }
		     dDistMean/=nColocPatN;
		     //display new average
		     IJ.log("Iteration "+ Integer.toString(nIterCount)+ " average distance: " + Double.toString(dDistMean*dPxtoNm) + " nm");
		     
		     if(nIterCount!=(int)(MaxRef-1))
		     {
		    	 // Refine the update-grid and reference grid	
		    	  O_add = refine_grid(O_add,Spacing);
		    	  Spacing[0]*=2;
		    	  Spacing[1]*=2;
		    	  O_ref = refine_grid(O_ref,Spacing);
		    	  dimX=(int) ((dMwidth+(Spacing[0]*3))/Spacing[0])+1;
		 	      dimY=(int) ((dMheight+(Spacing[1]*3))/Spacing[1])+1;
		     }
	     }	     
	     
	     return O_int;	     
	     
	}
	/** function that refines B-spline grid
	 * */
	
	double [][][] refine_grid(double [][][] O_trans, int[] Spacing)
	{
		int i,j,h,iter;
		double [][][]  O_newA, O_newB, O_new;
		double [] O_trans_linear,O_newA_linear,O_newB_linear;
		int sizeX, sizeY;
		int [] I,J,H;
		int [] ind;
		double [] P0;
		double [] P1;
		double [] P2;
		double [] P3;
		double [][] Pnew;
		int O_newA_szx,O_newB_szx,O_newB_szy;
		
		
		for (i=0;i<2;i++)
			Spacing[i]=Spacing[i]/2;
		sizeX=O_trans[0].length;
		sizeY=O_trans[0][0].length;
		
		// Refine B-spline grid in the x-direction
		O_newA_szx = ((sizeX-2)*2-1)+2;
		O_newA = new double[2][O_newA_szx][sizeY];
		
		
		
		O_newA_linear = new double [2*O_newA_szx*sizeY]; 
		I = new int[(sizeX-3)*sizeY*2];
		J = new int[(sizeX-3)*sizeY*2];
		H = new int[(sizeX-3)*sizeY*2];
		
		for(i=0;i<(sizeX-3);i++)
			for(j=0;j<sizeY;j++)
				for(h=0;h<2;h++)
				{
					I[i+j*(sizeX-3)+(sizeX-3)*sizeY*h]=i+1;
					J[i+j*(sizeX-3)+(sizeX-3)*sizeY*h]=j+1;
					H[i+j*(sizeX-3)+(sizeX-3)*sizeY*h]=h+1;
				}
		ind = sub2ind(O_trans,I,J,H);
		P0= new double[(sizeX-3)*sizeY*2];
		P1= new double[(sizeX-3)*sizeY*2];
		P2= new double[(sizeX-3)*sizeY*2];
		P3= new double[(sizeX-3)*sizeY*2];

		O_trans_linear = new double [2*sizeX*sizeY];
		for(i=0;i<sizeX;i++)
			for(j=0;j<sizeY;j++)
				for(h=0;h<2;h++)
				{

					O_trans_linear[i+j*(sizeX)+(sizeX)*sizeY*h]=O_trans[h][i][j];
				}
		
		for(iter=0;iter<(sizeX-3)*sizeY*2;iter++)
		{
			P0[iter]=O_trans_linear[ind[iter]-1];
			P1[iter]=O_trans_linear[ind[iter]];
			P2[iter]=O_trans_linear[ind[iter]+1];
			P3[iter]=O_trans_linear[ind[iter]+2];
			
		}
		Pnew = split_knots(P0,P1,P2,P3);
		
		for(i=0;i<(sizeX-3)*sizeY*2;i++)
			I[i]=1+(I[i]-1)*2;
		ind = sub2ind(O_newA,I,J,H);
		
		for(iter=0;iter<(sizeX-3)*sizeY*2;iter++)
		{
			O_newA_linear[ind[iter]-1]= Pnew[0][iter];
			O_newA_linear[ind[iter]]  = Pnew[1][iter];
			O_newA_linear[ind[iter]+1]= Pnew[2][iter];
			O_newA_linear[ind[iter]+2]= Pnew[3][iter];
			O_newA_linear[ind[iter]+3]= Pnew[4][iter];			
		}
		O_newB_szx = O_newA_szx;
		O_newB_szy = ((sizeY-2)*2-1)+2;
		O_newB=new double [2][O_newB_szx][O_newB_szy];
		O_newB_linear = new double[O_newB_szx*O_newB_szy*2];
		
		I = new int[(sizeY-3)*O_newA_szx*2];
		J = new int[(sizeY-3)*O_newA_szx*2];
		H = new int[(sizeY-3)*O_newA_szx*2];		
		for(j=0;j<O_newA_szx;j++)
			for(i=0;i<(sizeY-3);i++)
				for(h=0;h<2;h++)
				{
					J[j+i*O_newA_szx+(sizeY-3)*O_newA_szx*h]=j+1;
					I[j+i*O_newA_szx+(sizeY-3)*O_newA_szx*h]=i+1;
					H[j+i*O_newA_szx+(sizeY-3)*O_newA_szx*h]=h+1;
				}		
		ind = sub2ind(O_newA,J,I,H);
		P0= new double[(sizeY-3)*O_newA_szx*2];
		P1= new double[(sizeY-3)*O_newA_szx*2];
		P2= new double[(sizeY-3)*O_newA_szx*2];
		P3= new double[(sizeY-3)*O_newA_szx*2];

		for(iter=0;iter<(sizeY-3)*O_newA_szx*2;iter++)
		{
			P0[iter]=O_newA_linear[ind[iter]-1];
			P1[iter]=O_newA_linear[ind[iter]+O_newA_szx-1];
			P2[iter]=O_newA_linear[ind[iter]+O_newA_szx*2-1];
			P3[iter]=O_newA_linear[ind[iter]+O_newA_szx*2-1];

		}
		Pnew = split_knots(P0,P1,P2,P3);
		for(i=0;i<(sizeY-3)*O_newA_szx*2;i++)
			I[i]=1+(I[i]-1)*2;
		ind = sub2ind(O_newB,J,I,H);
		for(iter=0;iter<(sizeY-3)*O_newA_szx*2;iter++)
		{
			O_newB_linear[ind[iter]-1]= Pnew[0][iter];
			O_newB_linear[ind[iter] + O_newB_szx -1]  = Pnew[1][iter];
			O_newB_linear[ind[iter] + 2*O_newB_szx -1]= Pnew[2][iter];
			O_newB_linear[ind[iter] + 3*O_newB_szx -1]= Pnew[3][iter];
			O_newB_linear[ind[iter] + 4*O_newB_szx -1]= Pnew[4][iter];					
		}		
		int lsizeX,lsizeY;


		lsizeX=(int) ((dMwidth+(Spacing[0]*3))/Spacing[0])+1;
		lsizeY=(int) ((dMheight+(Spacing[1]*3))/Spacing[1])+1;
		O_new = new double[2][lsizeX][lsizeY];
		 // Set the final refined matrix

		for(i=0;i<lsizeX;i++)
			for(j=0;j<lsizeY;j++)
				for(h=0;h<2;h++)
					O_new[h][i][j] = O_newB_linear[i + j*O_newB_szx +h*O_newB_szx*O_newB_szy];				
		
		return O_new;
	}
	
	/**
	 * Splitting knots of b-spline
	 * */
	double [][] split_knots(double [] P0,double [] P1,double [] P2,double [] P3)
	{
		int i, nL=P0.length;
		double [][] Pnew = new double [5][nL];
		
		for (i=0;i<nL;i++)
		{
			Pnew[0][i]=(4.0/8.0)*(P0[i]+P1[i]);
			Pnew[1][i]=(1.0/8.0)*(P0[i]+6.0*P1[i]+P2[i]);
			Pnew[2][i]=(4.0/8.0)*(P1[i]+P2[i]);
			Pnew[3][i]=(1.0/8.0)*(P1[i]+6.0*P2[i]+P3[i]);
			Pnew[4][i]=(4.0/8.0)*(P2[i]+P3[i]);
		}
		return Pnew;
	}
	/**  Linear index from multiple subscripts.
%   SUB2IND is used to determine the equivalent single index
%   corresponding to a given set of subscript values.
%   Just a  simplified copy of corresponding matlab function
	 * */
	 
	int[] sub2ind(double [][][] grid, int [] I, int [] J, int [] H)
	{
		int l1,l2;
		l1 = grid[0].length;
		l2 = grid[0][0].length;
		int totN = I.length;
		int[] finX = new int[totN];
		int [] k = new int[3];
		int i;
		
		k[0] = 1;
		k[1] = l1;
		k[2] = l1*l2;
		for (i=0;i<totN;i++)
		{
			
			finX[i] += 1 + (I[i]-1)*k[0];
			finX[i] += (J[i]-1)*k[1];
			finX[i] += (H[i]-1)*k[2];
		}
		return finX;
	}
	/**
	 * Function transforming coordinates X in 2D using O_trans B-spline transform
	 * */
	double[][] bspline_transform_slow_2d(double [][][] O_trans, int [] Spacing, double [][] X)
	{
		int nPoints=X[0].length;
		double [][] preg = new double[2][nPoints];
		double [] x2 = new double [nPoints];
		double [] y2 = new double [nPoints];
		double [] u = new double [nPoints];
		double [] v = new double [nPoints];
		int [] ixs = new int [nPoints];
		int [] iys = new int [nPoints];
		int [] ix = new int [16*nPoints];
		int [] iy = new int [16*nPoints];
		double [][] Cx = new double[nPoints][16];
		double [][] Cy = new double[nPoints][16];
		int i,j;
		double [] CheckBound = new double [16*nPoints];	
		double [][] W;
		
		// Make row vectors of input coordinates
		x2=X[0]; 
		y2=X[1];
		
		// This code calculates for every coordinate in X, the indices of all
		// b-spline knots which have influence on the transformation value of
		// this point
		// Make indices of all neighborgh knots to every point
		
		int [] m = new int [16];
		int [] l = new int [16];
		
		for (i=0;i<4;i++)
			for (j=0;j<4;j++)
			{
				m[i+j*4]=i;
				l[i*4+j]=i;
			}
		
		for(i=0;i<nPoints;i++)
		{
			ixs[i]=(int)Math.floor(x2[i]/Spacing[0]);	
			iys[i]=(int)Math.floor(y2[i]/Spacing[1]);
		}
		for(i=0;i<nPoints;i++)
			for(j=0;j<16;j++)
			{
				ix[j+(i*16)] = ixs[i]+m[j];
				iy[j+(i*16)] = iys[i]+l[j];
			}
		// Points outside the bspline grid are set to the upper corner
		for(i=0;i<nPoints*16;i++)
		{
			if(ix[i]<0 | ix[i]>(dimX-1)|iy[i]<0|iy[i]>(dimY-1))
			{
				ix[i]=1;
				iy[i]=1;
				CheckBound[i]=0;
			}
			else
				CheckBound[i]=1;
		}
		// Look up the b-spline knot values in neighborhood of the points in (x2,y2)
		for(i=0;i<nPoints;i++)
			for(j=0;j<16;j++)
			{
				Cx[i][j] = O_trans[0][ix[j+(i*16)]][iy[j+(i*16)]]*CheckBound[j+(i*16)]; 
				Cy[i][j] = O_trans[1][ix[j+(i*16)]][iy[j+(i*16)]]*CheckBound[j+(i*16)]; 
			}
		// Calculate the b-spline interpolation constants u,v in the center cell
		// range between 0 and 1
		for(i=0;i<nPoints;i++)
		{
			v[i]  = (x2[i]-ixs[i]*Spacing[0])/Spacing[0];
			u[i]  = (y2[i]-iys[i]*Spacing[1])/Spacing[1];
		}
		// Get the b-spline coefficients in amatrix W, which contains
		// the influence of all knots on the points in (x2,y2)
		W=bspline_coefficients_2d(v,u);
		
		// Calculate the transformation of the points in (x2,y2) by the b-spline grid
		for(i=0;i<nPoints;i++)
			for(j=0;j<16;j++)
			{
				preg[0][i]+=W[j][i]*Cx[i][j];
				preg[1][i]+=W[j][i]*Cy[i][j];
			}			
		
		return preg;
	}
	
	/** Function makes a b-spline grid which minimizes the difference between the
	      corresponding points
	 * */
	double[][][] bspline_grid_fitting_2d(double[][][] O, int [] Spacing, double [][] R, double [][] X)
	{
		double[][][] O_trans;
		int [] gx = new int [nColocPatN];
		int [] gy = new int [nColocPatN];
		double [] ax = new double [nColocPatN];
		double [] ay = new double [nColocPatN];
		int [] ix, iy;
		//int [][] indexx = new int [16][nColocPatN];
		int [] indexx = new int [16*nColocPatN];
		//int [][] indexy = new int [16][nColocPatN];
		int [] indexy = new int [16*nColocPatN];
		int [][] index =new int [2][16*nColocPatN];
		double [][] W;
		double [] W2;
		double [][] WT;
		double [] WNx;
		double [] WNy;
		double [] S;
		double dTemp;
		double [][] numx;
		double [][] numy;
		double [][] dnum;
		int [] siz = new int [2];
		
		int i,j;
		// calculate which is the closest point on the lattice to the top-left
		// corner and find ratio's of influence between lattice point.
		for(i=0;i<nColocPatN;i++)
		{
			gx[i]=(int) Math.floor(X[0][i]/Spacing[0]);
			gy[i]=(int) Math.floor(X[1][i]/Spacing[1]);
		}

		// Calculate b-spline coordinate within b-spline cell, range 0..1
		for(i=0;i<nColocPatN;i++)
		{
			ax[i]=(X[0][i]-gx[i]*Spacing[0])/Spacing[0];
			ay[i]=(X[1][i]-gy[i]*Spacing[1])/Spacing[1]; 
		}
		for(i=0;i<nColocPatN;i++)
		{
			gx[i]+=2;
			gy[i]+=2;
		}
		//TODO add verification
		//if(any(ax<0)||any(ax>1)||any(ay<0)||any(ay>1)), error('grid error'), end;
	
		W = bspline_coefficients_2d(ax, ay);
		
		
		
		// Make indices of all neighborgh knots to every point
		ix = new int [16];
		iy = new int [16];
		
		for (i=-1;i<3;i++)
			for (j=0;j<4;j++)
			{
				ix[i+1+j*4]=i;
				iy[(i+1)*4+j]=i;
			}

		for(i=0;i<nColocPatN;i++)
			for(j=0;j<16;j++)
			{
				indexx[j+(i*16)] = gx[i]+ix[j];
				indexy[j+(i*16)] = gy[i]+iy[j];
			}
		//Limit to boundaries grid
		for(i=0;i<nColocPatN*16;i++)
		{
			if(indexx[i]<1)
				indexx[i]=1;
			if(indexx[i]>dimX)
				indexx[i]=dimX;
			if(indexy[i]<1)
				indexy[i]=1;
			if(indexy[i]>dimY)
				indexy[i]=dimY;	
			index[0][i]=indexx[i];
			index[1][i]=indexy[i];
		}
		
		// according too Lee et al. we update a numerator and a denumerator for
		// each knot. In our case we need two numerators, because our value is a
		// vector dy,dx. If we want to be able to add/remove keypoints, we need 
		// to store the numerators in seperate arrays.
		
		W2 = new double [16*nColocPatN];
		WT = new double [16][nColocPatN];
		WNx = new double [16*nColocPatN];
		WNy = new double [16*nColocPatN];
		S = new double [nColocPatN];
		for(i=0;i<nColocPatN;i++)
		{
			S[i]=0;
			for(j=0;j<16;j++)
			{
				dTemp=W[j][i]*W[j][i];
				W2[j+(i*16)]=dTemp;
				S[i]+=dTemp;
				WT[j][i]=dTemp*W[j][i];
			}
		}
		for(i=0;i<nColocPatN;i++)
			for(j=0;j<16;j++)
			{
				WNx[j+(i*16)]=WT[j][i]*R[0][i]/S[i];
				WNy[j+(i*16)]=WT[j][i]*R[1][i]/S[i];
			}
		siz[0]=dimX;
		siz[1]=dimY;
		numx=accumarray(index,WNx,siz);	
		numy=accumarray(index,WNy,siz);
		dnum=accumarray(index,W2,siz);
		
		// calculate actual values of knots from the numerator and denumerator that
		// we calculated previously
		// and update the b-spline transformation grid
		O_trans = new double [2][dimX][dimY];
		for(i=0;i<dimX;i++)
			for(j=0;j<dimY;j++)
			{
				//numx[i][j]=numx[i][j]/(dnum[i][j]+Double.MIN_VALUE);
				//numy[i][j]=numy[i][j]/(dnum[i][j]+Double.MIN_VALUE);
				
				O_trans[0][i][j] = O[0][i][j]+(numx[i][j]/(dnum[i][j]+Double.MIN_VALUE));
				O_trans[1][i][j] = O[1][i][j]+(numy[i][j]/(dnum[i][j]+Double.MIN_VALUE));
			}
		
		// Update the b-spline transformation grid
		//O_trans(:,:,1)=ux+O(:,:,1);
		//O_trans(:,:,2)=uy+O(:,:,2);
		
		
		return O_trans;
	}
	
	 double[][]  bspline_coefficients_2d(double[] u, double [] v)
	{
		double [][] Bu;
		double [][] Bv;
		double [][] W;
		int i, N, k, m;
		
		Bu=bspline_coefficients_1d(u);
		Bv=bspline_coefficients_1d(v);
		 
		N=u.length;
		
		W = new double [16][N];
		for(i=0;i<N;i++)
		{
				
			/*
			//probably smart cycle would be better
			W[0][i]=Bu[0][i]*Bv[0][i];
			W[1][i]=Bu[1][i]*Bv[0][i];
			W[2][i]=Bu[2][i]*Bv[0][i];
			W[3][i]=Bu[3][i]*Bv[0][i];
			
			W[4][i]=Bu[0][i]*Bv[1][i];
			W[5][i]=Bu[1][i]*Bv[1][i];
			W[6][i]=Bu[2][i]*Bv[1][i];
			W[7][i]=Bu[3][i]*Bv[1][i];
			
			W[8][i]=Bu[0][i]*Bv[2][i];
			W[9][i]=Bu[1][i]*Bv[2][i];
			W[10][i]=Bu[2][i]*Bv[2][i];
			W[11][i]=Bu[3][i]*Bv[2][i];
			
			W[12][i]=Bu[0][i]*Bv[3][i];
			W[13][i]=Bu[1][i]*Bv[3][i];
			W[14][i]=Bu[2][i]*Bv[3][i];
			W[15][i]=Bu[3][i]*Bv[3][i];
			*/
			
			for(k=0;k<4;k++)
				for(m=0;m<4;m++)
				{
					W[m+(k*4)][i]=Bu[m][i]*Bv[k][i];
				}

		}
		return W;
		
	}
	 /** Matlab's accumarray function */
	 double [][] accumarray (int [][] subs, double [] val, int [] sz)
	 {
		 double [][] retval;
		 int i;
	
		 //let's start accumulation
		 retval = new double [sz[0]][sz[1]];
		 
		 for(i=0;i<16*nColocPatN;i++)
		 {
			 retval[subs[0][i]-1][subs[1][i]-1]+=val[i];
		 }
		 
		 return retval;
	 }
	 /** the title of this function is self explanatory
	  * */
	 double[][]  bspline_coefficients_1d(double [] u)
	 {
		 double [][] W;
		 W = new double[4][u.length];		 
		 
		 for (int i=0;i<u.length;i++)
		 {
			W[0][i] = Math.pow(1-u[i],3)/6; 
			W[1][i] = ( 3*Math.pow(u[i],3) - 6*u[i]*u[i]+ 4)/6;
			W[2][i] = (-3*Math.pow(u[i],3) + 3*Math.pow(u[i],2) + 3*u[i] + 1)/6;
			W[3][i] = Math.pow(u[i],3)/6;
		 }
		 
		 return W;
	 }
	
	
	/** function locating particles images in both channels
	 *  returns number of colocalized particles  
	 * */
	int CCfindParticles()
	{
		ArrayList<double[]> colocTable = new ArrayList<double[]>();
		double [] xyfound;
		double dMinDist, dDist;
		int ind_warped;
		int i,j;
		int nRef, nWarp;
		nRef = xyref[0].length;
		nWarp = xywarp[0].length;
		int [] fwarpmark = new int [nWarp];
		//looking for colocalization
		for (i=0;i<nRef;i++)
		{
			if(SNRref[i]>dCCSNR)
			{
				dMinDist = Double.MAX_VALUE;
				ind_warped = -1;
				for (j=0;j<nWarp;j++)
				{
					if(SNRwarp[j]>dCCSNR)
					{
						//not marked as used yet
						if(fwarpmark[j]<3)
						{
							//on the same frame
							if((int)fref[i]== (int)fwarp[j])
							{
								dDist = Math.sqrt(Math.pow(xyref[0][i]-xywarp[0][j], 2)+Math.pow(xyref[1][i]-xywarp[1][j], 2));
								if (dDist<dCCDist)
								{
									if(dDist<dMinDist)
									{
										dMinDist=dDist;
										ind_warped = j;
										
									}
								}
							}
						}
					}
				}
				//let's see if found anything
				if(ind_warped>=0)
				{

					//mark particle as used
					fwarpmark[ind_warped]=5;
					xyfound = new double [4];
					xyfound[0]=xyref[0][i];
					xyfound[1]=xyref[1][i];
					xyfound[2]=xywarp[0][ind_warped];
					xyfound[3]=xywarp[1][ind_warped];
					colocTable.add(xyfound);
				}
			}
		}
		nColocPatN = colocTable.size();
		xymov = new double[2][nColocPatN];
		xystat = new double[2][nColocPatN];
		//let's add found point to new array
		for(i=0;i<nColocPatN;i++)
		{
			xystat[0][i] = colocTable.get(i)[0];
			xystat[1][i] = colocTable.get(i)[1];
			 xymov[0][i] = colocTable.get(i)[2];
			 xymov[1][i] = colocTable.get(i)[3];
		}
		return nColocPatN;
	}
    /**
     * function adjusting canvas size of two images to biggest dimension (separately in x and y)
     * 
     * */  
    void AdjustDimensions()
    {
		
		if(imp_refMax.getHeight() == imp_warpMax.getHeight())
		{
			dMheight = imp_refMax.getHeight();
		}
		else
		{
			if(imp_refMax.getHeight() < imp_warpMax.getHeight())
			{
				ImageStatistics stats = ImageStatistics.getStatistics(imp_refMax.getProcessor(), ImageStatistics.MEAN, null);
		        ImageProcessor ip2 = imp_refMax.getProcessor().createProcessor(imp_refMax.getWidth(), imp_warpMax.getHeight());
		        ip2.setValue(stats.mean);
		        ip2.fill();
		        ip2.insert(imp_refMax.getProcessor(), 0, 0);
		        
				imp_refMax.setProcessor(ip2);
				dMheight=imp_warpMax.getHeight();
			}
			else
			{
				ImageStatistics stats = ImageStatistics.getStatistics(imp_warpMax.getProcessor(), ImageStatistics.MEAN, null);
		        ImageProcessor ip2 = imp_warpMax.getProcessor().createProcessor(imp_warpMax.getWidth(), imp_refMax.getHeight());
		        ip2.setValue(stats.mean);
		        ip2.fill();
		        ip2.insert(imp_warpMax.getProcessor(), 0, 0);
		        
				imp_warpMax.setProcessor(ip2);
				dMheight=imp_refMax.getHeight();
			}
		}
		
		if(imp_refMax.getWidth() == imp_warpMax.getWidth())
		{
			dMwidth = imp_refMax.getWidth();	
		}
		else
		{
			if(imp_refMax.getWidth() < imp_warpMax.getWidth())
			{
				ImageStatistics stats = ImageStatistics.getStatistics(imp_refMax.getProcessor(), ImageStatistics.MEAN, null);
		        ImageProcessor ip2 = imp_refMax.getProcessor().createProcessor(imp_warpMax.getWidth(), imp_refMax.getHeight());
		        ip2.setValue(stats.mean);
		        ip2.fill();
		        ip2.insert(imp_refMax.getProcessor(), 0, 0);
		        
				imp_refMax.setProcessor(ip2);
				dMwidth=imp_warpMax.getWidth();
			}
			else
			{
				ImageStatistics stats = ImageStatistics.getStatistics(imp_warpMax.getProcessor(), ImageStatistics.MEAN, null);
		        ImageProcessor ip2 = imp_warpMax.getProcessor().createProcessor(imp_refMax.getWidth(), imp_warpMax.getHeight());
		        ip2.setValue(stats.mean);
		        ip2.fill();
		        ip2.insert(imp_warpMax.getProcessor(), 0, 0);
		        
				imp_warpMax.setProcessor(ip2);
				dMwidth=imp_refMax.getWidth();
				
			}
		}
		
		return;
    }    
    
    /** 
     * Function generates distortion map image (after translation correction)
     * 
     * */
    void GenerateMap()
    {
    	//transform coordinates
		xywarp = bspline_transform_slow_2d(O_trans, Spacing, xymov);
		//double [][][] vectorDir;
		//vectorDir = new double [2][2][nColocPatN];
		double xmin = Double.MAX_VALUE;
		double xmax = (-1)*Double.MAX_VALUE;
		double ymin = Double.MAX_VALUE;
		double ymax =(-1)*Double.MAX_VALUE;
		int i;
		
		double rmax = (-1)*Double.MAX_VALUE;
		double dDistance;
		double rScale;
		Overlay Directions = new Overlay();
		Arrow arrow;
		//get max displacemnt length
		// bounding rectangle
		for (i=0;i<nColocPatN;i++)
		{
			dDistance = Math.sqrt(Math.pow(xywarp[0][i]-xymov[0][i],2) + Math.pow(xywarp[1][i]-xymov[1][i],2));
			if(dDistance >rmax)
				rmax = dDistance;
		}
		
		rScale = 50/rmax;
		for (i=0;i<nColocPatN;i++)
		{
			xywarp[0][i] = xymov[0][i]+(xywarp[0][i]-xymov[0][i])*rScale;
			xywarp[1][i] = xymov[1][i]+(xywarp[1][i]-xymov[1][i])*rScale;
			if(xmin>xywarp[0][i])
				xmin = xywarp[0][i];
			if(xmin>xymov[0][i])
				xmin = xymov[0][i];
			if(ymin>xywarp[1][i])
				ymin = xywarp[1][i];
			if(ymin>xymov[1][i])
				ymin = xymov[1][i];
			if(xmax<xywarp[0][i])
				xmax = xywarp[0][i];
			if(xmax<xymov[0][i])
				xmax = xymov[0][i];
			if(ymax<xywarp[1][i])
				ymax = xywarp[1][i];
			if(ymax<xymov[1][i])
				ymax = xymov[1][i];
		}
					
		
		//xmin = xmin*rScale;
		//ymin = ymin*rScale;
		//xmax = xmax*rScale-xmin;
		//ymax = ymax*rScale-ymin;
		
		ColorProcessor imDirection = new ColorProcessor((int)Math.ceil(xmax), (int)Math.ceil(ymax));
		for (i=0;i<nColocPatN;i++)
		{
			//arrow = new Arrow(xymov[0][i]*rScale-xmin,xymov[1][i]*rScale-ymin, xywarp[0][i]*rScale-xmin,xywarp[1][i]*rScale-ymin);
			arrow = new Arrow(xymov[0][i]-xmin,xymov[1][i]-ymin, xywarp[0][i]-xmin,xywarp[1][i]-ymin);
			arrow.setStyle(Arrow.OPEN);
			arrow.setHeadSize(3.0);
			arrow.setStrokeWidth(1.0);
			arrow.setStrokeColor(Color.white);
			Directions.add(arrow);
		}
		ImagePlus map;
		map = new ImagePlus("Distortion map (minus rigid translation)", imDirection);
		map.setOverlay(Directions);
		map.show();
    }
    /** Shows composite image of two channels with overlayed
     * detected spots
     * */    
    void ShowRefPoints()
    {
    	ImagePlus [] images = new ImagePlus[2];    	
    	ImagePlus imFin;
    	Overlay imOV = new Overlay();
    	PointRoi imPointsMov;
    	PointRoi imPointsReg;
    	int i;
    	
    	images[1]=imp_refMax;
    	images[0]=imp_warpMax;
    	imFin = RGBStackMerge.mergeChannels(images, false);
    	imFin.setTitle("Reference channel(2) is green(o), warped channel(1) is red(+)");
    	
    	imPointsMov = new PointRoi(xymov[0][0]-dXYShift[0],xymov[1][0]-dXYShift[1]);
    	imPointsMov.setPointType(1);
    	//imPointsMov.setStrokeColor(Color.gr);
    	for(i=1;i<nColocPatN;i++)
    	{
    		imPointsMov.addPoint(xymov[0][i]-dXYShift[0], xymov[1][i]-dXYShift[1]);
    	}
    	imOV.add(imPointsMov);
    	
    	imPointsReg = new PointRoi(xystat[0][0],xystat[1][0]);
    	imPointsReg.setPointType(3);
    	//imPointsReg.setStrokeColor(Color.green);
    	for(i=1;i<nColocPatN;i++)
    	{
    		imPointsReg.addPoint(xystat[0][i], xystat[1][i]);
    	}
    	imOV.add(imPointsReg);
    	for(i=0;i<nColocPatN;i++)
    	{
    		imOV.add(new Line(xystat[0][i], xystat[1][i],xymov[0][i]-dXYShift[0], xymov[1][i]-dXYShift[1]));
    		
    	}
    	
    	
    	imFin.setOverlay(imOV);
    	imFin.show();
    	
    }
    /** functions stores current color calibration to prefs */
    void ColocCal_storeCalibration()
    {   
    	int i,j,h;
    	String sFieldName;
    	
    	Prefs.set("SiMoLOc.CC_dPxtoNm",dPxtoNm);
		Prefs.set("SiMoLOc.CC_calDate",reportDate);		
		Prefs.set("SiMoLOc.CC_dMheight",dMheight);
		Prefs.set("SiMoLOc.CC_dMwidth",dMwidth);
		Prefs.set("SiMoLOc.CC_Spacing0",Spacing[0]);
		Prefs.set("SiMoLOc.CC_Spacing1",Spacing[1]);
		Prefs.set("SiMoLOc.CC_dimX",dimX);
		Prefs.set("SiMoLOc.CC_dimY",dimY);
		Prefs.set("SiMoLOc.CC_dXYShift0",dXYShift[0]);
		Prefs.set("SiMoLOc.CC_dXYShift1",dXYShift[1]);
		
		for(i=0;i<dimX;i++)
			for(j=0;j<dimY;j++)
				for(h=0;h<2;h++)
				{
					sFieldName = "SiMoLOc.O_trans_" + Integer.toString(h)+"_"+ Integer.toString(i)+"_"+ Integer.toString(j);
					Prefs.set(sFieldName,O_trans[h][i][j]);
				}
		
		IJ.log("Calibration grid "+Integer.toString(dimX)+"x"+Integer.toString(dimY)+"x2 stored.");
		IJ.log("Calibration date: "+ reportDate);
    }
    
    /** functions loads current color calibration from prefs to the object */
    boolean ColocCal_getCalibration()
    {   
    	int i,j,h;
    	String sFieldName;
    	
    	
    	//read stored values of calibration
		
    	dPxtoNm=Prefs.get("SiMoLOc.CC_dPxtoNm", Double.NaN);
    	
		if(Double.isNaN(dPxtoNm))
		{
			IJ.error("There is no valid chromatic calibration present, please make one!");
			return false;
		}
		
		reportDate = Prefs.get("SiMoLOc.CC_calDate","");
    	dMheight = Prefs.get("SiMoLOc.CC_dMheight",0);
    	dMwidth = Prefs.get("SiMoLOc.CC_dMwidth",0);
    	Spacing = new int[2];
		Spacing[0] = (int) Prefs.get("SiMoLOc.CC_Spacing0",0);
		Spacing[1] = (int) Prefs.get("SiMoLOc.CC_Spacing1",0);
		dimX = (int) Prefs.get("SiMoLOc.CC_dimX",0);
		dimY = (int) Prefs.get("SiMoLOc.CC_dimY",0);
		dXYShift = new double [2];
		dXYShift[0] = Prefs.get("SiMoLOc.CC_dXYShift0",0);
		dXYShift[1] = Prefs.get("SiMoLOc.CC_dXYShift1",0);
		O_trans = new double[2][dimX][dimY];
		for(i=0;i<dimX;i++)
			for(j=0;j<dimY;j++)
				for(h=0;h<2;h++)
				{
					sFieldName = "SiMoLOc.O_trans_" + Integer.toString(h)+"_"+ Integer.toString(i)+"_"+ Integer.toString(j);
					O_trans[h][i][j] = Prefs.get(sFieldName,0);
				}
		
		return true;
    }
    
    /** functions prints current color calibration from prefs to the IJ log window*/
    void ColocCal_showCalibration()
    {   
    	int i,j,h;
    	String sFieldName;
    	String sVals;
    	
    	
    	//read stored values of calibration
		
    	dPxtoNm=Prefs.get("SiMoLOc.CC_dPxtoNm", Double.NaN);
    	
		if(Double.isNaN(dPxtoNm))
		{
			IJ.error("There is no valid chromatic calibration present, please make one!");
			return;
		}
		IJ.log("Pixel size: "+Double.toString(dPxtoNm)+" nm");
		reportDate = Prefs.get("SiMoLOc.CC_calDate","");
		IJ.log("Calibration date: "+ reportDate);
    	dMheight = Prefs.get("SiMoLOc.CC_dMheight",0);    	
    	dMwidth = Prefs.get("SiMoLOc.CC_dMwidth",0);
    	IJ.log("Image size: "+Integer.toString((int)dMwidth)+"x"+Integer.toString((int)dMheight)+" px");
    	Spacing = new int[2];
		Spacing[0] = (int) Prefs.get("SiMoLOc.CC_Spacing0",0);
		Spacing[1] = (int) Prefs.get("SiMoLOc.CC_Spacing1",0);
		IJ.log("Spacing: "+Integer.toString(Spacing[0])+"x"+Integer.toString(Spacing[1]));
		dimX = (int) Prefs.get("SiMoLOc.CC_dimX",0);
		dimY = (int) Prefs.get("SiMoLOc.CC_dimY",0);
		dXYShift = new double [2];
		dXYShift[0] = Prefs.get("SiMoLOc.CC_dXYShift0",0);
		dXYShift[1] = Prefs.get("SiMoLOc.CC_dXYShift1",0);	
		
		IJ.log("Shift in X: "+Integer.toString((int)dXYShift[0])+" px and Y:"+Integer.toString((int)dXYShift[1])+" px");
		IJ.log("Calibration grid size: "+Integer.toString(dimX)+"x"+Integer.toString(dimY)+"x2");
		O_trans = new double[2][dimX][dimY];
		for(h=0;h<2;h++)
		{
			IJ.log("Layer " + Integer.toString(h));
			
			for(j=0;j<dimY;j++)
			{
				sVals = "";
				for(i=0;i<dimX;i++)

				{
					sFieldName = "SiMoLOc.O_trans_" + Integer.toString(h)+"_"+ Integer.toString(i)+"_"+ Integer.toString(j);
					O_trans[h][i][j] = Prefs.get(sFieldName,0);
					sVals = sVals +Double.toString(O_trans[h][i][j])+" "; 
				}
				IJ.log(sVals);
			}
		}		
		return;
    }
    /** 
     * function applying b-spline grid transform to the table data
     * */
    void ColocCal_transformTable()
    {
    	double  dPxtoNmTable;
    	int i, nPatN;
    	
    	if (sml.ptable.getCounter()==0 || !sml.ptable.getHeadings()[0].equals("X_(px)"))
		{
			IJ.error("Not able to detect a valid 'Particles Table' for reconstruction, please load one.");
			return;
		}
    	IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		IJ.log("Applying chromatic correction to the table...");
		IJ.log("Calibration date: "+ reportDate);
		IJ.log("Shift in X: "+Integer.toString((int)dXYShift[0])+" px and Y:"+Integer.toString((int)dXYShift[1])+" px");
		IJ.log("Calibration grid size: "+Integer.toString(dimX)+"x"+Integer.toString(dimY)+"x2");
		
		dPxtoNmTable =  sml.ptable.getValueAsDouble(DOMConstants.Col_Xnm, 0)/sml.ptable.getValueAsDouble(DOMConstants.Col_X, 0);
		if(Math.abs(dPxtoNmTable-dPxtoNm)>0.1)
		{
			IJ.error("Pixel size of calibration " +Double.toString(dPxtoNm)+" nm is not equal to Results table pixel of "+Double.toString(dPxtoNmTable)+" nm");
			return;
		}
		//ok, seems pixel size is ok,
		//let's load results
		fref = sml.ptable.getColumnAsDoubles(DOMConstants.Col_FrameN);
		nPatN=fref.length;
		xyref = new double [2][nPatN];
		xyref[0] = sml.ptable.getColumnAsDoubles(DOMConstants.Col_X);
		xyref[1] = sml.ptable.getColumnAsDoubles(DOMConstants.Col_Y);
		
		//apply shift
    	for(i=0;i<nPatN;i++)
    	{
    		xyref[0][i] +=dXYShift[0];
    		xyref[1][i] +=dXYShift[1];
    	}
    	
    	//apply b-spline transform
    	 xyref = bspline_transform_slow_2d(O_trans, Spacing, xyref);
    	
    	IJ.showStatus("Updating table with corrected X and Y values...");
 		//update table
 		sml.ptable_lock.lock();
 		for(i=0;i<nPatN;i++)
 		{
			sml.ptable.setValue(DOMConstants.Col_Xnm, i, xyref[0][i]*dPxtoNmTable);
			sml.ptable.setValue(DOMConstants.Col_Ynm, i, xyref[1][i]*dPxtoNmTable);
			sml.ptable.setValue(DOMConstants.Col_X, i, xyref[0][i]);
			sml.ptable.setValue(DOMConstants.Col_Y, i, xyref[1][i]);
 		}
 		sml.ptable_lock.unlock();
		//sml_.ptable.updateResults();
		sml.showTable();
		IJ.showStatus("Updating table with corrected X and Y values...done");
    }
    /** 
     * function applying b-spline grid transform to the image/stack data
     * */
    void ColocCal_transformImage()
    {
    	
    	//int i_new, j_new;
    	//ImageProcessor ipfinal;
    
    	//let's get the image
		imp = IJ.getImage();
		
		if (imp==null)
		{
		    IJ.noImage();
		    return;
		}
		else if (imp.getType() != ImagePlus.GRAY8 && imp.getType() != ImagePlus.GRAY16 && imp.getType() != ImagePlus.GRAY32) 
		{
		    IJ.error("8, 16 or 32 bit greyscale image required");
		    return;
		}
		IJ.log(" --- DoM plugin version " + DOMConstants.DOMversion+ " --- ");
		IJ.log("Applying chromatic correction to image/stack...");
		IJ.log("Calibration date: "+ reportDate);
		IJ.log("Shift in X: "+Integer.toString((int)dXYShift[0])+" px and Y:"+Integer.toString((int)dXYShift[1])+" px");
		IJ.log("Calibration grid size: "+Integer.toString(dimX)+"x"+Integer.toString(dimY)+"x2");
		
		
		interpolate_image();
		
		return;
		//ipfinal = interpolate_bicubic();
		/*
		
		for(nSlice=0;nSlice<nStackSize;nSlice++)
		{
			imp.setSliceWithoutUpdate(nSlice+1);
			//ip = ;
			ipfinal = interpolate_image(imp.getProcessor().duplicate());
			if(imp.getType()==ImagePlus.GRAY8)
				finStack.addSlice(ipfinal.convertToByte(false));
			if(imp.getType()==ImagePlus.GRAY16)
				finStack.addSlice(ipfinal.convertToShort(false));
			else
				finStack.addSlice(ipfinal);
		}
		new ImagePlus(imp.getTitle()+" (chromatic corrected)",finStack).show();
		*/
		
    }
    
    /**
     * bicubic image interpolation 
     * */
    void interpolate_image()
    {
    	int trImHeight,trImWidth;
    	double [][] xyOrig;
    	double [][] xyTrans;
    	//double [][] xyBas;
    	int i,j;//,k;
    	double newx,newy;
    	double dVal;
    	int nSlice, nStackSize;
    	ImageStack finStack;
    	
    	//int i_new,j_new;
    	
    	ImageProcessor ip2;
    	FloatProcessor ipn;
    	
    	
    	// Make all x,y indices
		trImWidth = imp.getWidth();
		trImHeight = imp.getHeight();
		nStackSize = imp.getStackSize();
		finStack = new ImageStack(trImWidth,trImHeight);
	
		
		//ipn = (FloatProcessor) ip.convertToFloat();
		
		//account for translation		
		xyOrig=new double[2][trImWidth*trImWidth]; 
		for (i=0;i<trImWidth;i++)
			for (j=0;j<trImHeight;j++)
			{
				xyOrig[0][i+j*trImWidth] = i-dXYShift[0];
				xyOrig[1][i+j*trImWidth] = j-dXYShift[1];

			}
		
		IJ.showStatus("Calculating image transform...");
		// Calculate the transformation of all image coordinates by the b-spline grid
		xyTrans = bspline_transform_slow_2d(O_trans, Spacing, xyOrig);

		IJ.showStatus("Applying correction to images...");
		for(nSlice=0;nSlice<nStackSize;nSlice++)
		{
				IJ.showProgress(nSlice, nStackSize);
				imp.setSliceWithoutUpdate(nSlice+1);
				ipn = imp.getProcessor().duplicate().convertToFloatProcessor();
		
				ip2 = new FloatProcessor(trImWidth,trImHeight);
				
				//add bilinear? is it quicker?!!
				ipn.setInterpolationMethod(ImageProcessor.BICUBIC);
				for (i=0;i<trImWidth;i++)
					for (j=0;j<trImHeight;j++)
					{
						//newx=xyTrans[0][i+j*trImWidth];
						//newy=xyTrans[1][i+j*trImWidth];
						newx=2*xyOrig[0][i+j*trImWidth]-xyTrans[0][i+j*trImWidth];
						newy=2*xyOrig[1][i+j*trImWidth]-xyTrans[1][i+j*trImWidth];
						
						/*if(newx<0)
							newx=0;*/
						dVal = ipn.getInterpolatedPixel(newx, newy);
						//dVal = ip.getf((int)xyBas[0][i+j*trImWidth],(int)xyBas[1][i+j*trImWidth]);
						ip2.setf(i,j,(float)dVal);
					}
		
				if(imp.getType()==ImagePlus.GRAY8)
					finStack.addSlice(ip2.convertToByte(false));
				if(imp.getType()==ImagePlus.GRAY16)
					finStack.addSlice(ip2.convertToShort(false));
				else
					finStack.addSlice(ip2);
		}
		IJ.showProgress(1.0); // (nSlice, nStackSize);
		IJ.showStatus("Applying correction to images...done.");
		new ImagePlus(imp.getTitle()+" (chromatic corrected)",finStack).show();
		return;
				
		
    }
    /**
     * saves chromatic calibration to file
     * */   
    void ColocCal_Save()
    {
		SaveDialog sd = new SaveDialog("Save chromatic calibration", "chromatic_calibration", ".txt");
        String path = sd.getDirectory();
        int i,j,h;
        String sVals;
        if (path==null)
        	return;
        String filename = path+sd.getFileName();
        PrintWriter writer;
		try {
			writer = new PrintWriter(filename, "UTF-8");
			
			writer.println("CCdate\t"+reportDate);
			writer.println("CC_dPxtoNm\t"+Double.toString(dPxtoNm));
			writer.println("CC_dMheight\t"+Double.toString(dMheight));
			writer.println("CC_dMwidth\t"+Double.toString(dMwidth));
			writer.println("CC_Spacing0\t"+Integer.toString(Spacing[0]));
			writer.println("CC_Spacing1\t"+Integer.toString(Spacing[1]));
			writer.println("CC_dimX\t"+Integer.toString(dimX));
			writer.println("CC_dimY\t"+Integer.toString(dimY));
			writer.println("CC_dXYShift0\t"+Double.toString(dXYShift[0]));
			writer.println("CC_dXYShift1\t"+Double.toString(dXYShift[1]));
			writer.println("Coeff map");
			for(h=0;h<2;h++)
			{

				
				for(j=0;j<dimY;j++)
				{
					sVals = "";
					for(i=0;i<dimX;i++)
					{
						sVals = sVals +Double.toString(O_trans[h][i][j])+"\t"; 
					}
					writer.println(sVals);
				}
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
		IJ.log("Chromatic calibration saved.");
		return;
    }
    /**
     * loads chromatic calibration from file
     * */   
    void ColocCal_Load()
    {
    	
    	String sCalibration;
		String [] sCalSplittedRows;
		String [] sCalSplittedVals;
		String delimsn = "[\n]+";
		String delimst = "[\t]+";
		int i,j,h;
		String sFieldName;
		
		sCalibration = IJ.openAsString("");
		if(sCalibration==null || sCalibration.equals(""))
			return;
		sCalSplittedRows = sCalibration.split(delimsn);
		
		sCalSplittedVals =sCalSplittedRows[0].split(delimst);
		if(!(sCalSplittedVals[0].equals("CCdate")))
		{
			IJ.error("Not able to detect a valid chromatic calibration, try another file.");
			return;
		}
		
		Prefs.set("SiMoLOc.CC_calDat",sCalSplittedVals[1]);
		sCalSplittedVals =sCalSplittedRows[1].split(delimst);
		Prefs.set("SiMoLOc.CC_dPxtoNm",Double.parseDouble(sCalSplittedVals[1]));
		sCalSplittedVals =sCalSplittedRows[2].split(delimst);
		Prefs.set("SiMoLOc.CC_dMheight",Double.parseDouble(sCalSplittedVals[1]));
		sCalSplittedVals =sCalSplittedRows[3].split(delimst);
		Prefs.set("SiMoLOc.CC_dMwidth",Double.parseDouble(sCalSplittedVals[1]));
		sCalSplittedVals =sCalSplittedRows[4].split(delimst);
		Prefs.set("SiMoLOc.CC_Spacing0",Integer.parseInt(sCalSplittedVals[1]));	
		sCalSplittedVals =sCalSplittedRows[5].split(delimst);
		Prefs.set("SiMoLOc.CC_Spacing1",Integer.parseInt(sCalSplittedVals[1]));	
		
		sCalSplittedVals =sCalSplittedRows[6].split(delimst);
		dimX = Integer.parseInt(sCalSplittedVals[1]);
		Prefs.set("SiMoLOc.CC_dimX",dimX);	
		sCalSplittedVals =sCalSplittedRows[7].split(delimst);
		dimY = Integer.parseInt(sCalSplittedVals[1]);
		Prefs.set("SiMoLOc.CC_dimY",dimY);			
		
		sCalSplittedVals =sCalSplittedRows[8].split(delimst);
		Prefs.set("SiMoLOc.CC_dXYShift0",Double.parseDouble(sCalSplittedVals[1]));	
		sCalSplittedVals =sCalSplittedRows[9].split(delimst);
		Prefs.set("SiMoLOc.CC_dXYShift1",Double.parseDouble(sCalSplittedVals[1]));	
		
		for(h=0;h<2;h++)
		{

			for(j=0;j<dimY;j++)
			{
				sCalSplittedVals =sCalSplittedRows[(11+j)+h*dimY].split(delimst);
				for(i=0;i<dimX;i++)
				{	
					sFieldName = "SiMoLOc.O_trans_" + Integer.toString(h)+"_"+ Integer.toString(i)+"_"+ Integer.toString(j);
					Prefs.set(sFieldName,Double.parseDouble(sCalSplittedVals[i]));
				}
			}
			
		}
		IJ.log("Chromatic calibration loaded.");
		return;
    }
    /** Function outputs coordinates of reference
     * and warped channels to Results table*/
    void OutputResults()
    {
    	double dl;
    	sml.ptable.reset(); // erase particle table
    	
    	//transform coordinates
    	xywarp = bspline_transform_slow_2d(O_trans, Spacing, xymov);
    	
    	
    	sml.ptable_lock.lock();
 		for(int i=0;i<nColocPatN;i++)
 		{
 			sml.ptable.incrementCounter();
 			sml.ptable.addValue("X_(nm)_reference", (xystat[0][i])*dPxtoNm);
 			sml.ptable.addValue("Y_(nm)_reference", (xystat[1][i])*dPxtoNm);
 			sml.ptable.addValue("X_(nm)_warped_before", (xymov[0][i]-dXYShift[0])*dPxtoNm);
 			sml.ptable.addValue("Y_(nm)_warped_before", (xymov[1][i]-dXYShift[1])*dPxtoNm);
 			sml.ptable.addValue("X_(nm)_warped_translation", (xymov[0][i])*dPxtoNm);
 			sml.ptable.addValue("Y_(nm)_warped_translation", (xymov[1][i])*dPxtoNm);
 			sml.ptable.addValue("X_(nm)_warped_final", (xywarp[0][i])*dPxtoNm);
 			sml.ptable.addValue("Y_(nm)_warped_final", (xywarp[1][i])*dPxtoNm);
 			dl = Math.sqrt(Math.pow(xystat[0][i]-xymov[0][i]+dXYShift[0], 2)+Math.pow(xystat[1][i]-xymov[1][i]+dXYShift[1], 2))*dPxtoNm;
 			sml.ptable.addValue("dl(ref_-_before)", dl);
 			dl = Math.sqrt(Math.pow(xystat[0][i]-xymov[0][i], 2)+Math.pow(xystat[1][i]-xymov[1][i], 2))*dPxtoNm;
 			sml.ptable.addValue("dl(ref_-_translation)", dl);
 			dl = Math.sqrt(Math.pow(xystat[0][i]-xywarp[0][i], 2)+Math.pow(xystat[1][i]-xywarp[1][i], 2))*dPxtoNm;
 			sml.ptable.addValue("dl(ref_-_final)", dl);

 		}
 		sml.ptable_lock.unlock();
 		sml.showTable();
    	
    }
}
