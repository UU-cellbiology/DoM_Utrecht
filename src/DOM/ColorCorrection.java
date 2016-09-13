package DOM;

import java.awt.Color;
import java.util.ArrayList;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Undo;
import ij.gui.Arrow;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import ij.process.ColorProcessor;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.StackProcessor;


public class ColorCorrection implements PlugIn{
	
	ImagePlus imp;
	ImagePlus imp_refMax;
	ImagePlus imp_warpMax;
	
	ImagePlus mark_ref_imp;
	ImagePlus mark_warp_imp;

	
	ImageStack is_temp;
	Roi refROI, warpedROI;
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
	int dimX,dimY;
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
	
	//final NonBlockingGenericDialog nb = new NonBlockingGenericDialog("Fit Z-calibration"); 
	@Override
	public void run(String arg) {
	
		int i;
		
		IJ.register(ColorCorrection.class);
	
		//show parameters dialog
		if(!dlg.dglColorCalibration())
			return;
		dCCSNR= dlg.dCCSNR;
		dCCDist = dlg.dCCDist;
		
		NonBlockingGenericDialog nb = new NonBlockingGenericDialog("Choose image and ROI of REFERENCE channel");
		
		nb.addMessage("Choose calibration image and ROI region for REFERENCE channel");
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
		mark_ref_imp = imp;
		if(!(refROI==null))
		{
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

		dPxtoNm =  sml.ptable.getValueAsDouble(DOMConstants.Col_Xnm, 0)/sml.ptable.getValueAsDouble(DOMConstants.Col_X, 0);
		
		
		nb = new NonBlockingGenericDialog("Choose image and ROI of WARPED channel");
		nb.addMessage("Choose calibration image and ROI region for WARPED channel");
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
		mark_warp_imp = imp;
		if(!(warpedROI==null))
		{
			imp=imp.duplicate();
		}
		zProj = new ZProjector(imp);
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
		
		
		IJ.log("Aligning maximum projection pictures...");
		
		AdjustDimensions();
		dXYShift = calcShiftFFTCorrelation(imp_refMax.getProcessor(),imp_warpMax.getProcessor());
				
		IJ.showProgress(1.0);
		
		IJ.log("X shift: "+Double.toString(dXYShift[0])+ " px");
		IJ.log("Y shift: "+Double.toString(dXYShift[1]) + " px");
		IJ.log("Finding corresponding images of particles in both channels...");
		
		//update warped particles coordinates
		for(i=0;i<fwarp.length;i++)
		{
			
			xywarp[0][i]+=dXYShift[0];
			xywarp[1][i]+=dXYShift[1];
			//xymov[0][i]+=dXYShift[0];
			//xymov[1][i]+=dXYShift[1];
			
		}	
		
		//find images of particles in both channels
		CCfindParticles();
		
		//mark them on images (what if not possible?!)
		
		IJ.log("In total "+Integer.toString(nColocPatN)+" particles found.");
		
		
		IJ.log("Calculating B-spline transform...");
		O_trans = calculateBsplineTransform(5.0);
		IJ.log("..done.");
		
		if(dlg.bCCShowMap)
		{
			GenerateMap();
		}
		
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
	     IJ.log("Initial average distance: " + Double.toString(dDistMean*dPxtoNm) + " nm");

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
  	 // dimX=(int) ((dMwidth+(Spacing[0]*3))/Spacing[0])+1;
	   //   dimY=(int) ((dMheight+(Spacing[1]*3))/Spacing[1])+1;

		lsizeX=(int) ((dMwidth+(Spacing[0]*3))/Spacing[0])+1;
		lsizeY=(int) ((dMheight+(Spacing[1]*3))/Spacing[1])+1;
		O_new = new double[2][lsizeX][lsizeY];
		 // Set the final refined matrix
//		for(i=0;i<O_newB_szx;i++)
	//		for(j=0;j<O_newB_szy;j++)
				for(i=0;i<lsizeX;i++)
					for(j=0;j<lsizeY;j++)

				for(h=0;h<2;h++)
					O_new[h][i][j] = O_newB_linear[i + j*O_newB_szx +h*O_newB_szx*O_newB_szy];
					//O_newB[h][i][j] = O_newB_linear[i + j*O_newB_szx +h*O_newB_szx*O_newB_szy];
				
		//some extra verification test
		
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
				//TODO PLUS ONE IN ORIGINAL CODE  in both directions hmmmm ??? verify
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
		//int [][] indexx = new int [nColocPatN][16];
		for(i=0;i<nColocPatN;i++)
			for(j=0;j<16;j++)
			{
				//indexx[i][j] = gx[i]+ix[j];
				indexx[j+(i*16)] = gx[i]+ix[j];
				//indexy[j][i] = gy[i]+iy[j];
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
	
	//TODO add subpixel precision (?)
	/** Function calculates shift in X and Y between images using
	 * cross correlation calculated through FFT transform
	 * */
	public double [] calcShiftFFTCorrelation(ImageProcessor ip1, ImageProcessor ip2)
	{
		double [] xyshift = new double [2];
		
		
		int i,j;
		float dRefRe,dRefIm, dWarpRe, dWarpIm;

		//int imheight, imwidth;
		int imheightFFT, imwidthFFT;
	    FHT fht1,fht2;
	    
		//imheight=ip1.getHeight();
		//imwidth =ip1.getWidth();
		
	    fht1 = new FHT(pad(ip1));
	    fht2 = new FHT(pad(ip2));
	    //fht1.originalHeight=imheight;
	    //fht1.originalWidth=imwidth;
	    
	    //fht2.originalHeight=imheight;
	    //fht2.originalWidth=imwidth;
	    
	    fht1.transform();
	    fht2.transform();
	    ImageStack ct1 = fht1.getComplexTransform();
	    ImageStack ct2 = fht2.getComplexTransform();
	    ImageStack istemp;
		
	    FloatProcessor [] is_refFFT = new FloatProcessor[2];
		FloatProcessor [] is_warpFFT = new FloatProcessor[2];
		FloatProcessor [] is_FFTMult = new FloatProcessor[2];
		imheightFFT=ct1.getHeight();
		imwidthFFT =ct1.getWidth();
		for(i=0;i<2;i++)
		{
			is_refFFT[i]=(FloatProcessor) ct1.getProcessor(i+1);
		}
		for(i=0;i<2;i++)
		{
			is_warpFFT[i]=(FloatProcessor) ct2.getProcessor(i+1);
			is_FFTMult[i] = new FloatProcessor(imwidthFFT,imheightFFT);
		}
			
		//calculate convolution between two FFT
		for(i=0;i<imwidthFFT;i++)
			for(j=0;j<imheightFFT;j++)
			{
				dRefRe=is_refFFT[0].getf(i,j);
				dRefIm=is_refFFT[1].getf(i,j);
				
				dWarpRe=is_warpFFT[0].getf(i,j);
				//conjugate
				dWarpIm=(float) ((-1.0)*is_warpFFT[1].getf(i,j));
		
				
				is_FFTMult[0].setf(i,j,dRefRe*dWarpRe-dWarpIm*dRefIm);
				is_FFTMult[1].setf(i,j,dRefRe*dWarpIm+dWarpRe*dRefIm);
				
				
			}

		istemp = new ImageStack(imwidthFFT, imheightFFT);
		


		istemp.addSlice("Real", is_FFTMult[0]);
		istemp.addSlice("Imaginary", is_FFTMult[1]);
		
		
		//ct1 = doComplexInverseTransform(istemp,imwidth, imheight);
		ct1 = doComplexInverseTransform(istemp,0, 0);
		
		//ImagePlus imp1 = new ImagePlus("Complex of 1", ct1);
		//imp1.show();
		 
	    
		for(i=0;i<2;i++)
		{
			is_FFTMult[i]=(FloatProcessor) ct1.getProcessor(i+1);
		}
		
		int niMax=0,njMax=0;
		double nAbsMax=(-1.0)*Double.MAX_VALUE;
		double dVal;
		for(i=0;i<imwidthFFT;i++)
			for(j=0;j<imheightFFT;j++)
			{
				dVal=Math.sqrt(is_FFTMult[0].getf(i,j)*is_FFTMult[0].getf(i,j)+is_FFTMult[1].getf(i,j)*is_FFTMult[1].getf(i,j));
				
				if(dVal>0)
				{
					dVal++;
				}
				if(dVal>nAbsMax)
				{
					nAbsMax=dVal;
					niMax=i;
					njMax=j;
				}
			}
		
		if((niMax+1)>Math.floor(0.5*imwidthFFT))
		{
			niMax = niMax - imwidthFFT;
		}
		/*else
		{
			niMax = niMax;
		}*/
		if((njMax+1)>Math.floor(0.5*imheightFFT))
		{
			njMax = njMax - imheightFFT;
		}
		/*
		else
		{
			njMax = njMax;
		}*/
		
		xyshift[0]=(double)niMax;
		xyshift[1]=(double)njMax;
	    return xyshift;
		
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
	 * some parts from ij.plugin.FFT code
	 *  
	 * */
    ImageProcessor pad(ImageProcessor ip) 
    {
        int originalWidth = ip.getWidth();
        int originalHeight = ip.getHeight();
        int maxN = Math.max(originalWidth, originalHeight);
        int i = 2;
        while(i<maxN) i *= 2;
        if (i==maxN && originalWidth==originalHeight) 
        {
        
            return ip;
        }
        maxN = i;
   
        ImageStatistics stats = ImageStatistics.getStatistics(ip, ImageStatistics.MEAN, null);
        ImageProcessor ip2 = ip.createProcessor(maxN, maxN);
        ip2.setValue(stats.mean);
        ip2.fill();
        ip2.insert(ip, 0, 0);
  
        Undo.reset();
        //new ImagePlus("padded", ip2.duplicate()).show();
        return ip2;
    }
	/**	Swap quadrants 1 and 3 and 2 and 4 of the specified ImageProcessor 
	so the power spectrum origin is at the center of the image.
	<pre>
	    2 1
	    3 4
	</pre>
*/
    void swapQuadrants(ImageStack stack) {
        FHT fht = new FHT(new FloatProcessor(1, 1));
        for (int i=1; i<=stack.getSize(); i++)
            fht.swapQuadrants(stack.getProcessor(i));
    }
    
    /** Complex to Complex Inverse Fourier Transform
    *   @author Joachim Wesner
    */
    void c2c2DFFT(float[] rein, float[] imin, int maxN, float[] reout, float[] imout) 
    {
            FHT fht = new FHT(new FloatProcessor(maxN,maxN));
            float[] fhtpixels = (float[])fht.getPixels();
            // Real part of inverse transform
            for (int iy = 0; iy < maxN; iy++)
                  cplxFHT(iy, maxN, rein, imin, false, fhtpixels);
            fht.inverseTransform();
            // Save intermediate result, so we can do a "in-place" transform
            float[] hlp = new float[maxN*maxN];
            System.arraycopy(fhtpixels, 0, hlp, 0, maxN*maxN);
            // Imaginary part of inverse transform
            for (int iy = 0; iy < maxN; iy++)
                  cplxFHT(iy, maxN, rein, imin, true, fhtpixels);
            fht.inverseTransform();
            System.arraycopy(hlp, 0, reout, 0, maxN*maxN);
            System.arraycopy(fhtpixels, 0, imout, 0, maxN*maxN);
      }
    
    /** Build FHT input for equivalent inverse FFT
    *   @author Joachim Wesner
    */
    void cplxFHT(int row, int maxN, float[] re, float[] im, boolean reim, float[] fht) 
    {
            int base = row*maxN;
            int offs = ((maxN-row)%maxN) * maxN;
            if (!reim) {
                  for (int c=0; c<maxN; c++) {
                        int l =  offs + (maxN-c)%maxN;
                        fht[base+c] = ((re[base+c]+re[l]) - (im[base+c]-im[l]))*0.5f;
                  }
            } else {
                  for (int c=0; c<maxN; c++) {
                        int l = offs + (maxN-c)%maxN;
                        fht[base+c] = ((im[base+c]+im[l]) + (re[base+c]-re[l]))*0.5f;
                  }
            }
      }
    
    ImageStack doComplexInverseTransform(ImageStack stack, int width, int height ) 
    {
   
        int maxN = stack.getWidth();
        swapQuadrants(stack);
        float[] rein = (float[])stack.getPixels(1);
        float[] imin = (float[])stack.getPixels(2);
        float[] reout= new float[maxN*maxN];
        float[] imout = new float[maxN*maxN];
        c2c2DFFT(rein, imin, maxN, reout, imout);
        ImageStack stack2 = new ImageStack(maxN, maxN);
        swapQuadrants(stack);
        stack2.addSlice("Real", reout);
        stack2.addSlice("Imaginary", imout);
        stack2 = unpad(stack2,width,height);
        /*String name = WindowManager.getUniqueName(imp.getTitle().substring(10));
        ImagePlus imp2 = new ImagePlus(name, stack2);
        imp2.getProcessor().resetMinAndMax();
        imp2.show();*/
        return stack2;
    }
    
    ImageStack unpad(ImageStack stack, int width, int height) 
    {
      
 
        if (width==0 || height==0 || (width==stack.getWidth()&&height==stack.getHeight()))
            return stack;
        StackProcessor sp = new StackProcessor(stack, null);
        ImageStack stack2 = sp.crop(0, 0, width, height);
        return stack2;
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
     * Function generates distortion map after displacement
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
}
