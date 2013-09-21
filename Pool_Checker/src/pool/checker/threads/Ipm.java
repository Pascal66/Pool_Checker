package pool.checker.threads;

import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;

import org.opencv.calib3d.Calib3d;
import org.opencv.core.Core;
import org.opencv.core.Core.MinMaxLocResult;
import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.core.Point;
import org.opencv.core.Scalar;
import org.opencv.core.Size;

import android.content.Context;
import android.util.Log;

public class Ipm {
    // http://dpl.ceegs.ohio-state.edu/papers/lai_asprs_2009.pdf

    public class CameraInfos {
	private Context context;

	private static final String FILENAME = "CameraInfo";

	private double focalLengthX_arg = 3.5 * 81.35;

	// Intrinsic Matrice f3.5 ??
	// double[][] data = { { 2.8 * 100, 0, 360 }, { 0, 585.5, 240 },
	// { 0, 0, 1 } };
	// Mat cameraMatrix = new Mat(3, 3, CvType.CV_16S);
	// cameraMatrix.put(0, 0, data[0]);
	// cameraMatrix.put(1, 0, data[1]);
	// cameraMatrix.put(2, 0, data[2]);

	private double focalLengthY_arg = 3.5 * 167.3;
	private double opticalCenterX_arg = 360;
	private double opticalCenterY_arg = 240;
	private double cameraHeight_arg = 2;
	private double pitch_arg = 0;
	private double yaw_arg = 0;
	private double imageWidth_arg = 720;
	private double imageHeight_arg = 480;
	/** focal length in x and y */
	public Point focalLength = new Point();

	/**
	 * optical center coordinates in image frame (origin is (0,0) at top
	 * left)
	 */
	public Point opticalCenter = new Point();
	/** height of camera above ground */
	public double cameraHeight;
	/** pitch angle in radians (+ve downwards) */
	public double pitch;
	/** yaw angle in radians (+ve clockwise) */
	public double yaw;
	/** width of images */
	public double imageWidth;
	/** height of images */
	public double imageHeight;
	/** Near focus distance: */
	public Object near_dist;
	/** Far focus distance: */
	public double far_dist;
	/** Depth of field: */
	public double depth_of_field;
	/** Depth of focus: */
	public double depth_of_focus;
	/** Hyperfocal distance: */
	public double hyper_dist;
	/** FOV (horizontal) (degrees): */
	public double fov_h;
	/** FOV (vertical) (degrees): */
	public double fov_v;
	/** FOV (diagonal) (degrees): */
	public double fov_d;
	// private final String TAG = getClass().getSimpleName();
	public CameraInfos() {
	    Log.i(TAG, "Camerainfos init");
	    // init the strucure
	    // focalLength.x = focalLengthX_arg;
	    // focalLength.y = focalLengthY_arg;
	    // opticalCenter.x = opticalCenterX_arg;
	    // opticalCenter.y = opticalCenterY_arg;
	    // cameraHeight = cameraHeight_arg;
	    // pitch = pitch_arg * Math.PI / 180;
	    // yaw = yaw_arg * Math.PI / 180;
	    // imageWidth = imageWidth_arg;
	    // imageHeight = imageHeight_arg;

	    // Text = text_arg;
	}

	public CameraInfos getCameraInfo() {
	    CameraInfos cachedCameraInfo = new CameraInfos();
	    try {
		FileInputStream fileReader = context.openFileInput(FILENAME);
		ObjectInputStream CameraInfosReader = new ObjectInputStream(
			fileReader);
		cachedCameraInfo = (CameraInfos) CameraInfosReader.readObject();
	    } catch (Exception e) {
		Log.e("CameraInfos",
			"couldn't read cameraInfos: " + e.getMessage());
		// return null;
	    }

	    return cachedCameraInfo;
	}

	public void setCameraInfo() {

	    // TODO
	    /* try { FileOutputStream fileWriter = context.openFileOutput(FILENAME, Context.MODE_PRIVATE);
	     * ObjectOutputStream cameraInfosWriter = new ObjectOutputStream(fileWriter);
	     * cameraInfosWriter.writeObject(cameraInfos);
	     * fileWriter.close();
	     * } catch (FileNotFoundException e) {
	     * Log.e("CameraInfos",
	     * "file not found, couldn't cache cameraInfos: " + e.getMessage());
	     * } catch (IOException e) {
	     * Log.e("CameraInfos",
	     * "io exception, couldn't cache cameraInfos: " + e.getMessage()); } */
	}
    }
    /** Structure to hold the info about IPM transformation */
    public class IPMInfo {
	/** min and max x-value on ground in world coordinates */
	public double[] xLimits;

	/** min and max y-value on ground in world coordinates */
	public double[] yLimits;
	/**
	 * conversion between mm in world coordinate on the ground
	 * /*in x-direction and pixel in image
	 */
	public double xScale;
	/**
	 * conversion between mm in world coordinate on the ground
	 * /*in y-direction and pixel in image
	 */
	public double yScale;
	/** width */
	public int width;
	/* height */
	public int height;
	/**
	 * portion of image height to add to y-coordinate of
	 * /*vanishing point
	 */
	public double vpPortion;
	/** Left point in original image of region to make IPM for */
	public double ipmLeft;
	/** Right point in original image of region to make IPM for */
	public double ipmRight;
	/** Top point in original image of region to make IPM for */
	public double ipmTop;
	/** Bottom point in original image of region to make IPM for */
	public double ipmBottom;
	/** interpolation to use for IPM (0: bilinear, 1:nearest neighbor) */
	public int ipmInterpolation;
	public IPMInfo() {
	    xLimits = new double[2];
	    yLimits = new double[2];
	    xScale = 1;
	    yScale = 1;
	    width = 0;
	    height = 0;
	    vpPortion = 0;
	    ipmLeft = ipmRight = ipmTop = ipmBottom = 0;
	    ipmInterpolation = 0;
	}
    }
    // Note: internally, all computations are done in feet. Converted to meters
    // as needed for display.
    private static final double mm_2_feet = 1 / 304.8;

    private static final double feet_2_mm = 304.8;

    /***
     * \file InversePerspectiveMapping.cc \author Mohamed Aly
     * <malaa@caltech.edu> \date 11/29/2006
     */
    private final String TAG = getClass().getSimpleName();

    private static double VP_PORTION = 0.05;
    CameraInfos cameraParams = new CameraInfos();

    public IPMInfo ipmParams = new IPMInfo();

    /**
     * We are assuming the world coordinate frame center is at the camera, the
     * ground plane is at height -h, the X-axis is going right, the Y-axis is
     * going forward, the Z-axis is going up. The camera is looking forward with
     * optical axis in direction of Y-axis, with possible pitch angle (above or
     * below the Y-axis) and yaw angle (left or right). The camera coordinates
     * have the same center as the world, but the Xc-axis goes right, the
     * Yc-axis goes down, and the Zc-axis (optical cxis) goes forward. The
     * uv-plane of the image is such that u is horizontal going right, v is
     * vertical going down. The image coordinates uv are such that the pixels
     * are at half coordinates i.e. first pixel is (.5,.5) ...etc where the
     * top-left point is (0,0) i.e. the tip of the first pixel is (0,0)
     */
    // http://www.eee.nuigalway.ie/Research/car/documents/docualain_issc10.pdf
    public Ipm() {
    }

    public double areaTarget() {
	/** Angle from Normal to Surface Radians */
	double theta = /*Math.toRadians(66.2); //*/-cameraParams.pitch;
	/** Height above Surface (HO) metres */
	double HO = cameraParams.cameraHeight;
	/** HBE = Horizontal Half-Angle FOV Degres */
	double HBE = Math.toRadians(cameraParams.fov_h / 2);
	/** BHG1 = BHG2 = Vertical Half-Angle FOV Degres */
	double BHG = Math.toRadians(cameraParams.fov_v / 2);
	/** near metres */
	double AO = HO * Math.tan(theta - HBE);
	/** far metres */
	double DO = HO * Math.tan(theta + HBE);
	/** rayon = CO - AO metres */
	double b = (DO - AO) / 2;
	/** Point de visée metres */
	double CO = b + AO;

	double HC = Math.hypot(HO, CO);

	double a = HC * Math.tan(HBE);

	double G1O = HO * Math.tan(theta - BHG);
	double G2O = HO * Math.tan(theta + BHG);

	double CG1 = Math.abs(CO - G1O);
	double CG2 = Math.abs(CO - G2O);

	double area = (a / b)
		* (CG1 * Math.sqrt(b * b - CG1 * CG1) + Math.asin(Math.abs(CG1
			/ b)) * b * b)
		+ (a / b)
		* (CG2 * Math.sqrt(b * b - CG2 * CG2) + Math.asin(Math.abs(CG2
			/ b)) * b * b);
Log.i(TAG, "a: "+a+" b "+b);
	return area;
    }

    void computeDof(double distance) {
	/** Focus/Subject distance: */
	double dist = distance;
	dist = dist * 1000 * mm_2_feet;
	/** Circle of confusion (mm): */
	double coc_mm = 0.016;
	double coc = coc_mm * mm_2_feet;
	/** Lens focal length (mm): */
	double flen_mm = 30;
	double flen = flen_mm * mm_2_feet;
	/** Aperture: */
	double aperture = 4.3;

	/**
	 * Depth of field formula from:
	 * "Focal Encylopedia of Photography", Macmillan Company, New York, 1960
	 */
	double tmp = (aperture * coc * (dist - flen)) / (flen * flen);
	double dist_far = (tmp >= 1) ? 0. : dist / (1 - tmp);
	double dist_near = dist / (1 + tmp);
	double depth_of_field = (tmp < 1) ? (dist_far - dist_near) : 0.;
	double hyper_focal = (flen * flen) / (aperture * coc);

	double depth_of_focus = (flen * flen * dist)
		/ (hyper_focal * (dist - flen));

	// airy_disk_diam = 2.43932 x light_wave_length x fnum

	double airy_disk_diam = 2.43932 * (510 * .000001) * aperture;

	cameraParams.near_dist = feet_2_mm * (dist_near);
	cameraParams.far_dist = feet_2_mm * (dist_far);
	cameraParams.depth_of_field = feet_2_mm * (depth_of_field);
	cameraParams.depth_of_focus = feet_2_mm * (depth_of_focus) / 10;
	cameraParams.hyper_dist = feet_2_mm * (hyper_focal);
    }

    /**
     * This calculator computes the angular field of view for a lens of a specified focal length on a 35mm camera. For most modern consumer
     * level digital SLR cameras, a focal length multiplier of greater than 1 is appropriate because these cameras have a smaller sensor
     * than a 35mm negative. For these cameras a focal length multiplier of approximately 1.5-1.6 is appropriate. Note: This calculator
     * assumes a standard width/height image ratio of 3:2.
     */
    public void computeFov() {
	double flen = 60;
	double film_width = 36;
	double film_height = 24;

	double flen_mult = 1.5;

	// Account for focal length multiplier (actually, a film/sensor size multiplier)
	//film_width = film_width / flen_mult;
	//film_height = film_height / flen_mult;
	//double film_diag = (Math.sqrt((film_width * film_width) + (film_height * film_height)));

	double fov_h = (2 * Math.atan(film_width / (2 * flen)) * 180 / Math.PI);
	double fov_v = (2 * Math.atan(film_height / (2 * flen)) * 180 / Math.PI);
	//double fov_d = (2 * Math.atan(film_diag / (2 * flen)) * 180 / Math.PI);
	cameraParams.fov_h = (fov_h);
	cameraParams.fov_v = (fov_v);
	//cameraParams.fov_d = (fov_d);
//Log.e(TAG, "fov_h"+cameraParams.fov_h+" fovv "+cameraParams.fov_v);
    }

    private double GET_CV_MAT_ELEM(Mat mat, int row, int col) {
	// CV_8U and CV_8S -> byte[],
	// CV_16U and CV_16S -> short[],
	// CV_32S -> int[],
	// CV_32F -> float[],
	// CV_64F-> double[].
	// M(i,j) = ((float*)(mat->data.ptr + mat->step*i))[j]
	// double buff[] = new double[(int) (mat.total() * mat.channels())];
	double buff[] = mat.get(row, col/* , buff */);
	return buff[0];
    }

    /**
     * warpPerspective ?
     * This function returns the Inverse Perspective Mapping of the input image,
     * assuming a flat ground plane, and given the camera parameters.
     * 
     * @param inImage the input image
     * @param outImage the output image in IPM
     * @param ipmInfo the returned IPM info for the transformation
     * @param focalLength focal length (in x and y direction)
     * @param cameraInfo the camera parameters
     * @param outPoints indices of points outside the image
     * 
     * @return outImage the output image in IPM
     */
    public Mat mcvGetIPM(Mat inImage/* , Mat outImage, List<Point> outPoints */) {
	Mat outImage = new Mat();
	List<Point> outPoints = new ArrayList<Point>();

	// check input images types
	// CvMat inMat, outMat;
	// cvGetMat(inImage, &inMat);
	// cvGetMat(outImage, &outMat);
	// cout << CV_MAT_TYPE(inImage.type) << " " <<
	// CV_MAT_TYPE(FLOAT_MAT_TYPE) << " "
	// << CV_MAT_TYPE(INT_MAT_TYPE)<<"\n";

	// TODO if (!(CV_ARE_TYPES_EQ(inImage, outImage) &&
	// (CV_MAT_TYPE(inImage.type)==CV_MAT_TYPE(FLOAT_MAT_TYPE) ||
	// (CV_MAT_TYPE(inImage.type)==CV_MAT_TYPE(INT_MAT_TYPE)))))
	// {
	// cerr << "Unsupported image types in mcvGetIPM";
	// exit(1);
	// }

	// get size of input image
	double u, v;
	v = inImage.height();
	u = inImage.width();

	// get the vanishing point
	Point vp = new Point();
	vp = mcvGetVanishingPoint();
	vp.y = Math.max(0, vp.y);
	// vp.y = 30;

	// get extent of the image in the xfyf plane
	double eps = ipmParams.vpPortion * v;// VP_PORTION*v;
	ipmParams.ipmLeft = Math.max(0, ipmParams.ipmLeft);
	ipmParams.ipmRight = Math.min(u - 1, ipmParams.ipmRight);
	ipmParams.ipmTop = Math.max(vp.y + eps, ipmParams.ipmTop);
	ipmParams.ipmBottom = Math.min(v - 1, ipmParams.ipmBottom);
	double uvLimitsp[] = { vp.x, ipmParams.ipmRight, ipmParams.ipmLeft,
		vp.x, ipmParams.ipmTop, ipmParams.ipmTop, ipmParams.ipmTop,
		ipmParams.ipmBottom };
	// {vp.x, u, 0, vp.x,
	// vp.y+eps, vp.y+eps, vp.y+eps, v};
	Mat uvLimits = new Mat(2, 4, CvType.CV_32F/* FLOAT_MAT_TYPE *//* ,
								       * uvLimitsp */);
	uvLimits.put(0, 0, uvLimitsp);
	// get these points on the ground plane
	// Mat xyLimitsp = cvCreateMat(2, 4, CvType.CV_32F/*FLOAT_MAT_TYPE*/);
	Mat xyLimitsp = new Mat(2, 4, CvType.CV_32F/* FLOAT_MAT_TYPE */);

	Mat xyLimits = xyLimitsp;
	xyLimits = mcvTransformImage2Ground(uvLimits/* , xyLimits */);
	// SHOW_MAT(xyLimitsp, "xyLImits");

	// get extent on the ground plane
	Mat row1, row2;

	row1 = xyLimits.row(0); // cvGetRow(xyLimits, row1, 0);
	row2 = xyLimits.row(1); // cvGetRow(xyLimits, row2, 1);

	double xfMax, xfMin, yfMax, yfMin;
	// cvMinMaxLoc(row1, (double) xfMin, (double) xfMax, 0, 0, 0);
	// cvMinMaxLoc(row2, (double) yfMin, (double) yfMax, 0, 0, 0);
	MinMaxLocResult MinMaxLocX = Core.minMaxLoc(row1);
	MinMaxLocResult MinMaxLocY = Core.minMaxLoc(row2);
	xfMin = MinMaxLocX.minVal;
	xfMax = MinMaxLocX.maxVal;
	yfMin = MinMaxLocY.minVal;
	yfMax = MinMaxLocY.maxVal;

	int outRow = outImage.height();
	int outCol = outImage.width();

	double stepRow = (yfMax - yfMin) / outRow;
	double stepCol = (xfMax - xfMin) / outCol;

	// construct the grid to sample
	// Mat xyGrid = cvCreateMat(2, outRow*outCol,
	// CvType.CV_32F/*FLOAT_MAT_TYPE*/);
	Mat xyGrid = new Mat(2, outRow * outCol, CvType.CV_32F);

	int i, j;
	double x, y;
	// fill it with x-y values on the ground plane in world frame
	for (i = 0, y = yfMax - .5 * stepRow; i < outRow; i++, y -= stepRow)
	    for (j = 0, x = xfMin + .5 * stepCol; j < outCol; j++, x += stepCol) {
		SET_CV_MAT_ELEM(xyGrid, 0, i * outCol + j, x);
		SET_CV_MAT_ELEM(xyGrid, 1, i * outCol + j, y);
	    }
	// get their pixel values in image frame
	// Mat uvGrid = cvCreateMat(2, outRow*outCol,
	// CvType.CV_32F/*FLOAT_MAT_TYPE*/);
	Mat uvGrid = new Mat(2, outRow * outCol, CvType.CV_32F);

	uvGrid = mcvTransformGround2Image(xyGrid/* , uvGrid */);
	// now loop and find the nearest pixel value for each position
	// that's inside the image, otherwise put it zero
	double ui, vi;
	// get mean of the input image
	Scalar means = Core.mean(inImage); // CvScalar means = cvAvg(inImage);
	double[] mean = means.val;
	// generic loop to work for both float and int matrix types
	// #define MCV_GET_IPM(type)
	for (i = 0; i < outRow; i++)
	    for (j = 0; j < outCol; j++) {
		/* get pixel coordinates */
		ui = GET_CV_MAT_ELEM(uvGrid, /* FLOAT_MAT_ELEM_TYPE, */0, i
			* outCol + j);
		vi = GET_CV_MAT_ELEM(uvGrid, /* FLOAT_MAT_ELEM_TYPE, */1, i
			* outCol + j);
		/* check if out-of-bounds */
		/* if (ui<0 || ui>u-1 || vi<0 || vi>v-1) \ */
		if (ui < ipmParams.ipmLeft || ui > ipmParams.ipmRight
			|| vi < ipmParams.ipmTop || vi > ipmParams.ipmBottom) {
		    // CV_MAT_ELEM(outImage, type, i, j) = (type)mean;
		    outImage.put(i, j, mean);
		}
		/* not out of bounds, then get nearest neighbor */
		else {
		    /* Bilinear interpolation */
		    if (ipmParams.ipmInterpolation == 0) {
			int x1 = (int) (ui), x2 = (int) (ui + 1);
			int y1 = (int) (vi), y2 = (int) (vi + 1);
			x = ui - x1;
			y = vi - y1;
			double val = GET_CV_MAT_ELEM(inImage, y1, x1) * (1 - x)
				* (1 - y) + GET_CV_MAT_ELEM(inImage, y1, x2)
				* x * (1 - y)
				+ GET_CV_MAT_ELEM(inImage, y2, x1) * (1 - x)
				* y + GET_CV_MAT_ELEM(inImage, y2, x2) * x * y;
			SET_CV_MAT_ELEM(outImage, i, j, val);

		    }
		    /* nearest-neighbor interpolation */
		    else
			SET_CV_MAT_ELEM(
				outImage,
				i,
				j,
				GET_CV_MAT_ELEM(inImage, (int) (vi + .5),
					(int) (ui + .5)));
		}
		// TODO Out of memory error fpr outPoints
		/* if (outPoints != null
		 * && (ui < ipmParams.ipmLeft + 10
		 * || ui > ipmParams.ipmRight - 10
		 * || vi < ipmParams.ipmTop || vi > ipmParams.ipmBottom - 2))
		 * outPoints.add(new Point(j, i)); // .push_back(Point(j, i)); */}
	// TODO if
	// (/*CV_MAT_TYPE*/(inImage.type())==CvType.CV_32F/*FLOAT_MAT_TYPE*/)
	// {
	// MCV_GET_IPM(CvType.CV_32F/*FLOAT_MAT_TYPE*/)
	// }
	// else
	// {
	// MCV_GET_IPM(CvType.CV_16U/*INT_MAT_ELEM_TYPE*/)
	// }
	// return the ipm info
	ipmParams.xLimits[0] = GET_CV_MAT_ELEM(xyGrid, 0, 0);
	ipmParams.xLimits[1] = GET_CV_MAT_ELEM(xyGrid, 0, (outRow - 1) * outCol
		+ outCol - 1);
	ipmParams.yLimits[1] = GET_CV_MAT_ELEM(xyGrid, 1, 0);
	ipmParams.yLimits[0] = GET_CV_MAT_ELEM(xyGrid, 1, (outRow - 1) * outCol
		+ outCol - 1);
	ipmParams.xScale = 1 / stepCol;
	ipmParams.yScale = 1 / stepRow;
	ipmParams.width = outCol;
	ipmParams.height = outRow;

	// clean
	xyLimits.release(); // cvReleaseMat(xyLimitsp);
	xyGrid.release(); // cvReleaseMat(xyGrid);
	uvGrid.release(); // cvReleaseMat(uvGrid);

	return outImage;
    }

    /**
     * Gets the extent of the image on the ground plane given the camera
     * parameters
     * 
     * @param cameraInfo the input camera info
     * @param ipmInfo the IPM info containing the extent on ground plane:
     *        xLimits & yLimits only are changed
     */
    public void mcvGetIPMExtent() {
	// get size of input image
	double /* FLOAT */u, v;
	v = cameraParams.imageHeight;
	u = cameraParams.imageWidth;

	// get the vanishing point
	Point /* FLOAT_POINT2D */vp;
	vp = mcvGetVanishingPoint();
	vp.y = Math.max(0, vp.y);

	// get extent of the image in the xfyf plane
	double eps = VP_PORTION * v;
	double uvLimitsp[] = { vp.x, u, 0, vp.x, vp.y + eps, vp.y + eps,
		vp.y + eps, v };
	Mat uvLimits = new Mat(2, 4, CvType.CV_32F/* , uvLimitsp */);
	uvLimits.put(0, 0, uvLimitsp);
	// get these points on the ground plane
	// Mat xyLimitsp = cvCreateMat(2, 4, CvType.CV_32F);

	Mat xyLimitsp = new Mat(2, 4, CvType.CV_32F);
	Mat xyLimits = xyLimitsp;
	xyLimits = mcvTransformImage2Ground(uvLimits);
	// SHOW_MAT(xyLimitsp, "xyLImits");
	// Log.i(TAG, xyLimits.dump());
	// get extent on the ground plane
	Mat row1, row2;
	row1 = xyLimits.row(0); // cvGetRow(xyLimits, row1, 0);
	row2 = xyLimits.row(1); // cvGetRow(xyLimits, row2, 1);
	double xfMax, xfMin, yfMax, yfMin;
	// cvMinMaxLoc(row1, (double) xfMin, (double) xfMax, 0, 0, 0);
	// cvMinMaxLoc(row2, (double) yfMin, (double) yfMax, 0, 0, 0);

	MinMaxLocResult MinMaxLocX = Core.minMaxLoc(row1);
	MinMaxLocResult MinMaxLocY = Core.minMaxLoc(row2);
	xfMin = MinMaxLocX.minVal;
	xfMax = MinMaxLocX.maxVal;
	yfMin = MinMaxLocY.minVal;
	yfMax = MinMaxLocY.maxVal;
	// return
	ipmParams.xLimits[0] = xfMin;
	ipmParams.xLimits[1] = xfMax;
	ipmParams.yLimits[1] = yfMax;
	ipmParams.yLimits[0] = yfMin;
    }

    /**
     * Computes the vanishing point in the image plane uv. It is the point of
     * intersection of the image plane with the line in the XY-plane in the
     * world coordinates that makes an angle yaw clockwise (form Y-axis) with
     * Y-axis
     * 
     * @param cameraInfo the input camera parameter
     * @return the computed vanishing point in image frame
     */
    public Point /* FLOAT_POINT2D */mcvGetVanishingPoint() {
	
	// get the vp in world coordinates
	double[] vpp = {
		Math.sin(cameraParams.yaw) / Math.cos(cameraParams.pitch),
		Math.cos(cameraParams.yaw) / Math.cos(cameraParams.pitch), 0 };
	Mat vp = new Mat(3, 1, CvType.CV_32F/* , vpp */);
	vp.put(0, 0, vpp);

	// transform from world to camera coordinates
	//
	// rotation matrix for yaw
	double /* FLOAT_MAT_ELEM_TYPE */tyawp[] = { Math.cos(cameraParams.yaw),
		-Math.sin(cameraParams.yaw), 0, Math.sin(cameraParams.yaw),
		Math.cos(cameraParams.yaw), 0, 0, 0, 1 };
	Mat tyaw = new Mat(3, 3, CvType.CV_32F/* , tyawp */);
	tyaw.put(0, 0, tyawp);
	// rotation matrix for pitch
	double /* FLOAT_MAT_ELEM_TYPE */tpitchp[] = { 1, 0, 0, 0,
		-Math.sin(cameraParams.pitch), -Math.cos(cameraParams.pitch),
		0, Math.cos(cameraParams.pitch), -Math.sin(cameraParams.pitch) };
	Mat transform = new Mat(3, 3, CvType.CV_32F/* , tpitchp */);
	transform.put(0, 0, tpitchp);
	// combined transform
	Core.multiply(transform, tyaw, transform); // cvMatMul(transform, tyaw,
						   // transform);

	// Intrinsic Matrice
	// transformation from (xc, yc) in camra coordinates
	// to (u,v) in image frame matrix to shift optical center and focal
	// length
	double /* FLOAT_MAT_ELEM_TYPE */t1p[] = { cameraParams.focalLength.x,
		0, cameraParams.opticalCenter.x, 0, cameraParams.focalLength.y,
		cameraParams.opticalCenter.y, 0, 0, 1 };
	Mat t1 = new Mat(3, 3, CvType.CV_32F/* FLOAT_MAT_TYPE *//* , t1p */);
	t1.put(0, 0, t1p);

	// combine transform
	Core.multiply(t1, transform, transform);
	// cvMatMul(t1, transform, transform);
	// transform
	// Core.multiply(vp, transform, vp); //inversion des matricces a
	// multiplier? non // cvMatMul(transform, vp, vp);
	Mat src3 = new Mat();
	Core.gemm(transform, vp, 1, src3, 0, vp, 0);
	//
	// clean and return
	//
	Point ret = new Point();
	ret.x = vp.get(0, 0)[0]; // ret.x = cvGetReal1D(vp, 0);
	ret.y = vp.get(1, 0)[0]; // ret.y = cvGetReal1D(vp, 1);
	return ret;
    }

    /**
     * Initializes the cameraInfo structure with data read from the conf file
     * 
     * @param fileName the input camera conf file name
     * @param cameraInfo the returned camera parametrs struct TOOO
     * 
     * @return
     */
    public CameraInfos mcvInitCameraInfo(Context context) {
	// parsed camera data
	// read the data
	Log.i(TAG, "Camerainfos init");

	// assert (cameraInfoParser_configfile(fileName, camInfo, 0, 1, 1) ==
	// 0);
	// if (cameraParams == null) {

	// init the strucure
	cameraParams.focalLength.x = cameraParams.focalLengthX_arg;
	cameraParams.focalLength.y = cameraParams.focalLengthY_arg;
	cameraParams.opticalCenter.x = cameraParams.opticalCenterX_arg;
	cameraParams.opticalCenter.y = cameraParams.opticalCenterY_arg;
	cameraParams.cameraHeight = cameraParams.cameraHeight_arg;
	cameraParams.pitch = cameraParams.pitch_arg * Math.PI / 180;
	cameraParams.yaw = cameraParams.yaw_arg * Math.PI / 180;
	cameraParams.imageWidth = cameraParams.imageWidth_arg;
	cameraParams.imageHeight = cameraParams.imageHeight_arg;
	// }
	return cameraParams;
    }

    /**
     * Converts a point from IPM pixel coordinates into world coordinates
     * 
     * @param point in/out point
     * @param ipmInfo the ipm info from mcvGetIPM
     * 
     * @return
     */
    Point mcvPointImIPM2World(Point point) {
	// x-direction
	point.x /= ipmParams.xScale;
	point.x += ipmParams.xLimits[0];
	// y-direction
	point.y /= ipmParams.yScale;
	point.y = ipmParams.yLimits[1] - point.y;

	return point;
    }

    /**
     * Scales the cameraInfo according to the input image size
     * 
     * @param cameraInfo the input/return structure
     * @param size the input image size
     */
    void mcvScaleCameraInfo(Size size) {
	// compute the scale factor
	double scaleX = size.width / cameraParams.imageWidth;
	double scaleY = size.height / cameraParams.imageHeight;
	// scale
	cameraParams.imageWidth = size.width;
	cameraParams.imageHeight = size.height;
	cameraParams.focalLength.x *= scaleX;
	cameraParams.focalLength.y *= scaleY;
	cameraParams.opticalCenter.x *= scaleX;
	cameraParams.opticalCenter.y *= scaleY;
    }

    /**
     * Transforms points from the ground plane (z=-h) in the world frame into
     * points on the image in image frame (uv-coordinates)
     * 
     * @param inPoints 2xN array of input points on the ground in world
     *        coordinates
     * @param outPoints 2xN output points in on the image in image coordinates
     * @param cameraInfo the camera parameters
     * 
     * @return
     */
    Mat mcvTransformGround2Image(Mat inPoints/* , Mat outPoints */) {
	Mat outPoints = new Mat();
	// add two rows to the input points
	// Mat inPoints3 = cvCreateMat(inPoints.rows()+1, inPoints.cols(),
	// inPoints.type() /*cvGetElemType(inPoints)*/);
	Mat inPoints3 = new Mat(inPoints.rows() + 1, inPoints.cols(),
		inPoints.type());
	// copy inPoints to first two rows
	Mat inPoints2, inPointsr3;

	inPoints2 = inPoints3.rowRange(0, 2); // cvGetRows(inPoints3, inPoints2,
					      // 0, 2);
	inPointsr3 = inPoints3.row(2); // cvGetRow(inPoints3, inPointsr3, 2);

	Scalar sCH = new Scalar(-cameraParams.cameraHeight);
	inPointsr3.setTo(sCH); // cvSet(inPointsr3,
			       // cvRealScalar(-cameraInfo.cameraHeight));

	inPoints.copyTo(inPoints2); // cvCopy(inPoints, inPoints2);
	// create the transformation matrix
	double c1 = Math.cos(cameraParams.pitch);
	double s1 = Math.sin(cameraParams.pitch);
	double c2 = Math.cos(cameraParams.yaw);
	double s2 = Math.sin(cameraParams.yaw);
	double matp[] = {
		cameraParams.focalLength.x * c2 + c1 * s2
			* cameraParams.opticalCenter.x,
		-cameraParams.focalLength.x * s2 + c1 * c2
			* cameraParams.opticalCenter.x,
		-s1 * cameraParams.opticalCenter.x,

		s2
			* (-cameraParams.focalLength.y * s1 + c1
				* cameraParams.opticalCenter.y),
		c2
			* (-cameraParams.focalLength.y * s1 + c1
				* cameraParams.opticalCenter.y),
		-cameraParams.focalLength.y * c1 - s1
			* cameraParams.opticalCenter.y,

		c1 * s2, c1 * c2, -s1 };
	Mat mat = new Mat(3, 3, CvType.CV_32F/* FLOAT_MAT_TYPE *//* , matp */);
	mat.put(0, 0, matp);
	// multiply
	// Core.multiply(mat, inPoints3, inPoints3); // cvMatMul(mat, inPoints3,
	// inPoints3);
	Mat src3 = new Mat();
	Core.gemm(mat, inPoints3, 1, src3, 0, inPoints3, 0);

	// divide by last row of inPoints4
	for (int i = 0; i < inPoints.cols(); i++) {
	    double div = GET_CV_MAT_ELEM(inPointsr3, 0, i);
	    SET_CV_MAT_ELEM(inPoints3, 0, i, GET_CV_MAT_ELEM(inPoints3, 0, i)
		    / div);
	    SET_CV_MAT_ELEM(inPoints3, 1, i, GET_CV_MAT_ELEM(inPoints3, 1, i)
		    / div);
	}
	// put back the result into outPoints
	inPoints2.copyTo(outPoints); // cvCopy(inPoints2, outPoints);
	// clear
	inPoints3.release(); // cvReleaseMat(inPoints3);

	return outPoints;
    }

    /**
     * Transforms points from the image frame (uv-coordinates) into the real
     * world frame on the ground plane (z=-height)
     * 
     * @param inPoints input points in the image frame
     * @param cameraInfo the input camera parameters
     * @return
     */
    Mat mcvTransformImage2Ground(Mat inPoints/* ,Mat outPoints */) {
	Mat outPoints = new Mat();
	// add two rows to the input points
	// Mat inPoints4 = cvCreateMat(inPoints.rows()+2, inPoints.cols(),
	// inPoints.type()/*cvGetElemType(inPoints)*/);
	Mat inPoints4 = new Mat(inPoints.rows() + 2, inPoints.cols(),
		inPoints.type());
	// copy inPoints to first two rows
	Mat inPoints2, inPoints3, inPointsr4, inPointsr3;

	inPoints2 = inPoints4.rowRange(0, 2); // cvGetRows(inPoints4, inPoints2,
					      // 0, 2);
	inPoints3 = inPoints4.rowRange(0, 3); // cvGetRows(inPoints4, inPoints3,
					      // 0, 3);
	inPointsr3 = inPoints4.row(2); // cvGetRow(inPoints4, inPointsr3, 2);
	inPointsr4 = inPoints4.row(3); // cvGetRow(inPoints4, inPointsr4, 3);

	Scalar One = new Scalar(1);
	inPointsr3.setTo(One); // cvSet(inPointsr3, cvRealScalar(1));
	inPoints.copyTo(inPoints2); // cvCopy(inPoints, inPoints2);
	// create the transformation matrix
	double c1 = Math.cos(cameraParams.pitch);
	double s1 = Math.sin(cameraParams.pitch);
	double c2 = Math.cos(cameraParams.yaw);
	double s2 = Math.sin(cameraParams.yaw);
	double matp[] = {
		-cameraParams.cameraHeight * c2 / cameraParams.focalLength.x,
		cameraParams.cameraHeight * s1 * s2
			/ cameraParams.focalLength.y,
		(cameraParams.cameraHeight * c2 * cameraParams.opticalCenter.x / cameraParams.focalLength.x)
			- (cameraParams.cameraHeight * s1 * s2
				* cameraParams.opticalCenter.y / cameraParams.focalLength.y)
			- cameraParams.cameraHeight * c1 * s2,

		cameraParams.cameraHeight * s2 / cameraParams.focalLength.x,
		cameraParams.cameraHeight * s1 * c2
			/ cameraParams.focalLength.y,
		(-cameraParams.cameraHeight * s2 * cameraParams.opticalCenter.x / cameraParams.focalLength.x)
			- (cameraParams.cameraHeight * s1 * c2
				* cameraParams.opticalCenter.y / cameraParams.focalLength.y)
			- cameraParams.cameraHeight * c1 * c2,

		0,
		cameraParams.cameraHeight * c1 / cameraParams.focalLength.y,
		(-cameraParams.cameraHeight * c1 * cameraParams.opticalCenter.y / cameraParams.focalLength.y)
			+ cameraParams.cameraHeight * s1,

		0,
		-c1 / cameraParams.focalLength.y,
		(c1 * cameraParams.opticalCenter.y / cameraParams.focalLength.y)
			- s1, };
	Mat mat = new Mat(4, 3, CvType.CV_32F/* , matp */);
	mat.put(0, 0, matp);
	// multiply
	// Core.multiply(mat, inPoints3, inPoints4); // cvMatMul(mat, inPoints3,
	// // inPoints4);
	Mat src3 = new Mat();
	Core.gemm(mat, inPoints3, 1, src3, 0, inPoints4, 0);
	// divide by last row of inPoints4
	for (int i = 0; i < inPoints.cols(); i++) {
	    double div = GET_CV_MAT_ELEM(inPointsr4, 0, i);
	    SET_CV_MAT_ELEM(inPoints4, 0, i, GET_CV_MAT_ELEM(inPoints4, 0, i)
		    / div);
	    SET_CV_MAT_ELEM(inPoints4, 1, i, GET_CV_MAT_ELEM(inPoints4, 1, i)
		    / div);
	}
	// put back the result into outPoints
	inPoints2.copyTo(outPoints); // cvCopy(inPoints2, outPoints);
	// clear
	inPoints4.release(); // cvReleaseMat(inPoints4);

	return outPoints;
    }

    /**
     * Converts from IPM pixel coordinates into world coordinates
     * 
     * @param inMat input matrix 2xN
     * @param outMat output matrix 2xN
     * @param ipmInfo the ipm info from mcvGetIPM
     */
    Mat mcvTransformImIPM2Ground(Mat inMat/* , Mat outMat */) {
	Mat outMat = new Mat();

	Mat mat;
	mat = outMat;
	if (inMat != mat) {
	    inMat.copyTo(mat); // cvCopy(inMat, mat);
	}

	// work on the x-direction i.e. first row
	Mat row;
	row = mat.row(0); // cvGetRow(mat, row, 0);
	// cvConvertScale(row, row, 1. / ipmInfo.xScale, ipmInfo.xLimits[0]);
	row.convertTo(row, row.type(), 1. / ipmParams.xScale,
		ipmParams.xLimits[0]);
	// work on y-direction
	row = mat.row(1); // cvGetRow(mat, row, 1);
	// cvConvertScale(row, row, -1. / ipmInfo.yScale, ipmInfo.yLimits[1]);
	row.convertTo(row, row.type(), -1. / ipmParams.yScale,
		ipmParams.yLimits[1]);

	return row;

    }

    /**
     * Converts from IPM pixel coordinates into Image coordinates
     * 
     * @param inMat input matrix 2xN
     * @param outMat output matrix 2xN
     * @param ipmInfo the ipm info from mcvGetIPM
     * @param cameraInfo the camera info
     * 
     * @return
     */
    Mat mcvTransformImIPM2Im(Mat inMat/* , Mat outMat */) {
	Mat outMat = new Mat();
	// convert to world coordinates
	outMat = mcvTransformImIPM2Ground(inMat/* , outMat */);

	// convert to image coordinates
	outMat = mcvTransformGround2Image(outMat/* , outMat */);

	return outMat;
    }

    private Mat SET_CV_MAT_ELEM(Mat mat, int row, int col, double value) {
	// CV_8U and CV_8S -> byte[],
	// CV_16U and CV_16S -> short[],
	// CV_32S -> int[],
	// CV_32F -> float[],
	// CV_64F-> double[].
	// M(i,j) = ((float*)(mat->data.ptr + mat->step*i))[j]
	// double buff[] = new double[(int) (mat.total() * mat.channels())];
	mat.put(row, col, value);
	return mat;
    }
}
