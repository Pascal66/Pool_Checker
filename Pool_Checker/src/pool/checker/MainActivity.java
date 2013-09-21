/* Copyright (c) 2011 Regents of the University of California.
 * All rights reserved.
 * This software was developed at the University of California, Irvine.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in
 * the documentation and/or other materials provided with the
 * distribution.
 * 3. All advertising materials mentioning features or use of this
 * software must display the following acknowledgment:
 * "This product includes software developed at the University of
 * California, Irvine by Nicolas Oros, Ph.D.
 * (http://www.cogsci.uci.edu/~noros/)."
 * 4. The name of the University may not be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 * 5. Redistributions of any form whatsoever must retain the following
 * acknowledgment:
 * "This product includes software developed at the University of
 * California, Irvine by Nicolas Oros, Ph.D.
 * (http://www.cogsci.uci.edu/~noros/)."
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL THE UNIVERSITY OR THE PROGRAM CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * ***************************
 * *************************************************************************** */
package pool.checker;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import net.priveyes.ia.controller.FuzzyController;
import net.priveyes.ia.controller.PIDController;

import org.opencv.android.BaseLoaderCallback;
import org.opencv.android.CameraBridgeViewBase.CvCameraViewFrame;
import org.opencv.android.CameraBridgeViewBase.CvCameraViewListener2;
import org.opencv.android.JavaCameraView;
import org.opencv.android.LoaderCallbackInterface;
import org.opencv.android.OpenCVLoader;
import org.opencv.core.Core;
import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.core.MatOfByte;
import org.opencv.core.MatOfFloat;
import org.opencv.core.MatOfInt;
import org.opencv.core.MatOfPoint;
import org.opencv.core.MatOfPoint2f;
import org.opencv.core.Point;
import org.opencv.core.Rect;
import org.opencv.core.Scalar;
import org.opencv.core.Size;
import org.opencv.highgui.Highgui;
import org.opencv.imgproc.Imgproc;

import pool.checker.threads.Ipm;
import pool.checker.threads.Ipm.CameraInfos;
import pool.checker.threads.SensorFusion;
import pool.checker.threads.Solar_Thread;
import pool.checker.threads.ThreadManager;
import android.annotation.SuppressLint;
import android.app.Activity;
import android.content.Intent;
import android.location.Location;
import android.media.ExifInterface;
import android.net.Uri;
import android.os.Bundle;
import android.os.Environment;
import android.os.Handler;
import android.os.Message;
import android.util.DisplayMetrics;
import android.util.Log;
import android.view.Display;
import android.view.MotionEvent;
import android.view.View;
import android.view.Window;
import android.view.WindowManager;
import android.widget.CompoundButton;
import android.widget.CompoundButton.OnCheckedChangeListener;
import android.widget.LinearLayout.LayoutParams;
import android.widget.Toast;

import com.hourlyweather.forecast.ForecastFetcherTask.ResultCallback;
import com.hourlyweather.forecast.ForecastHour;
import com.hourlyweather.forecast.HourlyForecast;

/**
 * This is the main entry point of the program. It handles the connection to the
 * Zeemotes, as well as starting all the other threads.
 * 
 * @author Alex Flynn
 */
public class MainActivity extends Activity implements CvCameraViewListener2,
	SensorFusion.ResultCallback, Solar_Thread.ResultCallback,
	ResultCallback {

    private static final int DRAW_RESULT_BITMAP = 10;
    private static final int DRAW_SUN_POS = 20;
    private static final int DRAW_NORTH_POS = 30;
    private static final int WEATHER_READY = 40;

    // The threshold value for the lower and upper color limits
    public static final double THRESHOLD_LOW = 35;
    public static final double THRESHOLD_HIGH = 35;

    private final String TAG = getClass().getSimpleName();

    private Handler mUiHandler;
    private ThreadManager threadManager;

    public static final int FORECAST_HOUR_SPAN = 36;

    /** {@link Mat} to hold the OpenCV images. */
    private Mat mRgba, mGray, mHsv, mPool, mShadow, mSun, mIntermediateMat;

    private List<Byte> byteStatus;
    private List<Point> outPixels, cornersThis, cornersPrev;

    private Mat matOpFlowPrev;
    private Mat matOpFlowThis;
    private Mat mMOP2fptsSafe;

    private MatOfPoint2f mMOP2fptsPrev;
    private MatOfPoint2f mMOP2fptsThis;

    private MatOfPoint2f mMOP2fHarris;

    MatOfByte mMOBStatus;
    MatOfFloat mMOFerr;

    private Point pt, pt1, pt2;
    private Scalar colorRed, colorGreen, colorYellow, colorNew;
    private Scalar colorBlack, colorWhite;
    private int iGFFTMax = 40, iFileOrdinal = 0;

    private boolean STOP_THREAD, inUse = false;

    /** The IP address of the server to send feed. */
    // private String ip_address;
    /** The host activity. */
    // private MainActivity app;
    /** The Utils thread to control the camera. */
    // private UtilsThread utils;
    /** A list of MatOfPointcorners provided by OpenCV. */
    private MatOfPoint MatOfPointcorners;
    private FuzzyController fuzzyController;
    private PIDController Pid_Controller;

    private ThreadManager manager;

    private Mat lines;
    private Mat img_edge;
    private Size ksize;
    private Scalar min_green;
    private Scalar max_green;
    private Mat mask, contour;
    private MatOfPoint mMOP;
    private Scalar max_shadow;
    private Scalar min_shadow;
    private Scalar min_sun;
    private Scalar max_sun;
    public Scalar max_pool, scalar_max, min_pool_rgb;
    public Scalar min_pool, scalar_min, max_pool_rgb;

    private Point mSelectedPoint = null;

    private Scalar mLowerColorLimit;
    private Scalar mUpperColorLimit;

    private Scalar red_max;
    private Scalar red_min;
    private Scalar red_rot;

    double perimeter = 0;

    public static Point SunPos;
    public static Scalar Orientation;
    public Location location;

    private static JavaCameraView mOpenCvCameraView;
    private double Temp;

    private Mat mErodeKernel;
    private boolean bShootNow = false, bDisplayTitle = true;

    private boolean bNewColor = false;
    private String string, sShotText;

    private double dTextScaleFactor;
    String filename;

    private double[] vecHoughLines;
    Rect MinRoi;
    private Point center;
    private WindowManager mWindowManager;
    private Display mDisplay;
    private List<Mat> pChannel;
    private CameraInfos cameraParam;
    private double x1;
    private double x2;
    private double y1;
    private double y2;

    static List<Double> temperature;

    private Mat mMixedMat;
    private Mat mCompMat;

    Ipm ipm = new Ipm();

    ExifInterface mExif;

    // Handle initialization error
    // Dynamic init
    //static { if (!OpenCVLoader.initDebug()) { } }

    /** This class will receive a callback once the OpenCV library is loaded. */
    // private static final class OpenCVLoaderCallback extends BaseLoaderCallback {
    public BaseLoaderCallback OpenCVLoaderCallback = new BaseLoaderCallback(
	    this) {
	// private Context mContext;

	/* public OpenCVLoaderCallback(Context context) { super(context); mContext = context; } */
	@Override
	public void onManagerConnected(int status) {
	    switch (status) {
	    case LoaderCallbackInterface.SUCCESS: {
		Log.i(TAG, "OpenCV loaded successfully");

		mOpenCvCameraView.enableView();

		if (threadManager == null) {
		    threadManager = new ThreadManager(MainActivity.this);

		} else {
		    threadManager.restartAll();
		}
	    }
		break;
	    default:
		super.onManagerConnected(status);
		break;
	    }
	}
    };
    private ExifInterface exif;
    private Mat mPoolRgb;
    private ColorProbability colorProbability;
    private boolean contourOk;
    private CondensationAlgo condensation;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
	super.onCreate(savedInstanceState);

	/** Create a Handler that we can post messages to so we avoid having to use anonymous Runnables and runOnUiThread() instead */
	mUiHandler = new Handler(getMainLooper(), new UiCallback());
    }

    @Override
    protected void onResume() {
	super.onResume();
	Log.i(TAG, "onResume Trying to load OpenCV library");
	try {

	    if (!OpenCVLoader.initAsync(OpenCVLoader.OPENCV_VERSION_2_4_6,
		    this, OpenCVLoaderCallback))
		Toast.makeText(this, "cannot connect to opencv ",
			Toast.LENGTH_SHORT).show();
	} catch (Exception e) {
	    Log.e(TAG,
		    "can not load opencv, an exception has been thrown while.");
	    e.printStackTrace();
	} finally {
	    // getWindow().setFlags(WindowManager.LayoutParams.FLAG_FULLSCREEN,
	    // WindowManager.LayoutParams.FLAG_FULLSCREEN);
	    getWindow()
		    .addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);
	    requestWindowFeature(Window.FEATURE_NO_TITLE);
	    setContentView(R.layout.main);

	    mOpenCvCameraView = (JavaCameraView) findViewById(R.id.java_surface_view);

	    mOpenCvCameraView.setVisibility(View.VISIBLE);
	    // mOpenCvCameraView.setCvCameraViewListener(this);

	    mOpenCvCameraView.setLayoutParams(new LayoutParams(
		    android.view.ViewGroup.LayoutParams.MATCH_PARENT,
		    android.view.ViewGroup.LayoutParams.MATCH_PARENT));

	    mOpenCvCameraView.setCvCameraViewListener(this);
	}
    }

    @Override
    public void onDestroy() {
	super.onDestroy();
	if (mOpenCvCameraView != null)
	    mOpenCvCameraView.disableView();
    }

    @Override
    protected void onPause() {
	super.onPause();
	if (mOpenCvCameraView != null)
	    mOpenCvCameraView.disableView();
    }

    @Override
    protected void onStop() {
	super.onStop();
	if (threadManager != null) {
	    threadManager.stopAll();
	}
	this.finish();
    }

    private class toggleListener implements OnCheckedChangeListener {
	@Override
	public void onCheckedChanged(CompoundButton buttonView,
		boolean isChecked) {
	    if (isChecked) {
		if (threadManager == null) {
		    threadManager = new ThreadManager(MainActivity.this);

		} else {
		    threadManager.restartAll();
		}
		Toast.makeText(MainActivity.this, "Start streaming",
			Toast.LENGTH_SHORT).show();
	    } else {
		threadManager.stopAll();
		Toast.makeText(MainActivity.this, "Stop streaming",
			Toast.LENGTH_SHORT).show();
	    }
	}
    }

    // public void onResultMatrixReady(Bitmap resultBitmap) {
    // mUiHandler.obtainMessage(DRAW_RESULT_BITMAP,
    // resultBitmap).sendToTarget();
    // }
    @Override
    public void onOrientationReady(Scalar orientation) {
	mUiHandler.obtainMessage(DRAW_NORTH_POS, orientation).sendToTarget();
    }

    @Override
    public void onSolarReady(Point mSunPos) {
	mUiHandler.obtainMessage(DRAW_SUN_POS, mSunPos).sendToTarget();
    }

    @Override
    public void onWeatherReady(HourlyForecast forecast) {
	mUiHandler.obtainMessage(WEATHER_READY, forecast).sendToTarget();
    }

    /** This Handler callback is used to */
    private class UiCallback implements Handler.Callback {

	@Override
	public boolean handleMessage(Message message) {
	    if (message.what == DRAW_NORTH_POS) {
		setOrientation((Scalar) message.obj);
	    }
	    if (message.what == DRAW_SUN_POS) {
		setSun((Point) message.obj);
		// Log.w(TAG, "Received msg " + message.obj.toString());
	    }
	    if (message.what == WEATHER_READY) {
		setWeather((HourlyForecast) message.obj);
	    }
	    // Log.w(TAG, "Received msg " + message.obj.toString());

	    return true;
	}
    }

    @Override
    public boolean onTouchEvent(final MotionEvent event) {
	bShootNow = false;
	if (event.getEventTime() - event.getDownTime() > 1000)
	    bShootNow = true;
	// else
	// pickColorFromTap(event);
	return false; // don't need more than one touch event
    }

    /** Calculate the point in the preview frame from the tap point on the screen */
    private void pickColorFromTap(MotionEvent event) {

	mLowerColorLimit = null;
	mUpperColorLimit = null;
	mSelectedPoint = new Point(event.getX(), event.getY());
	bNewColor = true;
    }

    @Override
    public void onCameraViewStarted(int width, int height) {
	// TODO Auto-generated method stub
	Log.i(TAG, "init");

	// byteColourTrackCentreHue = new byte[3];
	// green = 60 // mid yellow 27
	// byteColourTrackCentreHue[0] = 27;
	// byteColourTrackCentreHue[1] = 100;
	// byteColourTrackCentreHue[2] = (byte)255;
	byteStatus = new ArrayList<Byte>();

	// channels = new ArrayList<Integer>();
	// channels.add(0);
	colorRed = new Scalar(255, 0, 0, 255);
	colorGreen = new Scalar(0, 255, 0, 255);
	colorYellow = new Scalar(255, 255, 0, 255);
	colorNew = new Scalar(200, 0, 0, 255);
	colorBlack = new Scalar(0);
	colorWhite = new Scalar(255, 255, 255);
	// contours = new ArrayList<MatOfPoint>();
	// corners = new ArrayList<Point>();
	cornersThis = new ArrayList<Point>();
	cornersPrev = new ArrayList<Point>();

	// faces = new MatOfRect();

	// histSize = new MatOfInt(25);

	// iHueMap = new ArrayList<Integer>();
	// iHueMap.add(0);
	// iHueMap.add(0);
	// lines = new Mat();

	// mApproxContour = new MatOfPoint2f();
	// mContours = new Mat();
	// mHist = new Mat();
	// mGray = new Mat();
	// mHSVMat = new Mat();
	// mIntermediateMat = new Mat();
	// mMatRed = new Mat();
	// mMatGreen = new Mat();
	// mMatBlue = new Mat();
	// mMatRedInv = new Mat();
	// mMatGreenInv = new Mat();
	// mMatBlueInv = new Mat();
	// MOIone = new MatOfInt(0);

	// MOFrange = new MatOfFloat(0f, 256f);
	// mMOP2f1 = new MatOfPoint2f();
	// mMOP2f2 = new MatOfPoint2f();
	mMOP2fptsPrev = new MatOfPoint2f();
	mMOP2fptsThis = new MatOfPoint2f();
	mMOP2fptsSafe = new MatOfPoint2f();
	mMOFerr = new MatOfFloat();
	mMOBStatus = new MatOfByte();
	MatOfPointcorners = new MatOfPoint();
	// mRgba = new Mat();
	// mROIMat = new Mat();
	// mFaceDest = new Mat();
	// mFaceResized = new Mat();
	// matFaceHistogramPrevious = new Mat();
	// matFaceHistogramThis= new Mat();
	matOpFlowThis = new Mat();
	matOpFlowPrev = new Mat();

	pt = new Point(0, 0);
	// pt1 = new Point (0, 0);
	pt2 = new Point(0, 0);

	fuzzyController = new FuzzyController();

	ksize = new Size(5, 5);
	min_green = new Scalar(30, 50, 50, 255); // Pelouse
	max_green = new Scalar(70, 256, 250, 255); //

	/* One characteristic of shadow is your intensity (channel V) is low. It
	 * lie in a Range between lower and high Value, you will need discover
	 * your. Other characteristic is the saturation is lower. Frequently the
	 * object (caster) have high saturation. */
	min_shadow = new Scalar(113, 52, 0, 255); // Shadow pool
	max_shadow = new Scalar(146, 157, 104, 255); //

	min_sun = new Scalar(55, 5, 130, 255); // Sun pool
	max_sun = new Scalar(75, 115, 250, 255); //

	min_pool = new Scalar(113, 52, 104, 255); // Piscine 113- 52-104
	max_pool = new Scalar(146, 157, 254, 255); // 146-157-254

	min_pool_rgb = new Scalar(0, 33, 55, 255);
	max_pool_rgb = new Scalar(150, 205, 231, 255); // R55 G110 B130 Shadow

	red_min = new Scalar(114, 135, 135, 255);
	red_max = new Scalar(142, 255, 255, 255);
	int rotation = 128 - 255; // 255 = red
	red_rot = new Scalar(rotation, 0, 0, 255);

	DisplayMetrics dm = this.getResources().getDisplayMetrics();
	int densityDpi = dm.densityDpi;
	dTextScaleFactor = (densityDpi / 240.0) * 0.9;
	string = "";

	mRgba = new Mat(height, width, CvType.CV_8UC1);
	mGray = new Mat();
	mIntermediateMat = new Mat(height, width, CvType.CV_8UC1);

	mMixedMat = new Mat(height, width, CvType.CV_8UC1);
	mCompMat = new Mat(height, width, CvType.CV_8UC1);

	mHsv = new Mat();
	mPool = new Mat();
	mPoolRgb = new Mat();
	mShadow = new Mat();

	mErodeKernel = Imgproc
		.getStructuringElement(Imgproc.MORPH_CROSS, ksize);

	center = new Point(mRgba.cols() / 2, mRgba.rows() / 2);

	try {
	    cameraParam = ipm.mcvInitCameraInfo(this);
	    ipm.computeFov();

	    Log.i(TAG, "OK");
	} catch (Exception e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	}
	outPixels = new ArrayList<Point>();

	pt1 = new Point();
	pt2 = new Point();

	temperature = new ArrayList<Double>();

	colorProbability = new ColorProbability();
	
//	condensation = new CondensationAlgo();
    }

    @Override
    public void onCameraViewStopped() {
	// TODO Auto-generated method stub

    }

    @Override
    public Mat onCameraFrame(CvCameraViewFrame inputFrame) {
	// Voir Gaussian->SobelX-Threshold->Houghline probabilistic->Blob
	/* inUse = true; mRgba.put(0, 0, datadata); inUse = false; */
	mRgba = inputFrame.rgba();
	//mGray = inputFrame.gray();
	//if(contourOk) mRgba = colorProbability.simplestColorBalance(mRgba, 0.01);
	// ipm.mcvGetIPM(mGray, mGray, outPixels);

	// Mat mG = new Mat();
	// Mat mR = new Mat();
	// Mat mB = new Mat();

	// New Dilate Pascal
	/* // Core.extractChannel(mRgba, mR, 0);
	 * Core.extractChannel(mRgba, mG, 1);
	 * // Core.extractChannel(mRgba, mB, 2);
	 * 
	 * mG = NewDilate(mG);
	 * Core.bitwise_not(mRgba, mRgba, mG);
	 * 
	 * // Core.insertChannel(mG, mRgba, 1); */

	/* try { Log.e(TAG, "Solar" + SunPos.toString()); Log.e(TAG, "Orient" +
	 * Orientation.toString()); } catch (Exception e) { // TODO
	 * Auto-generated catch block // e.printStackTrace(); } */
	// Imgproc.cvtColor(mRgba, mGray, Imgproc.COLOR_RGBA2GRAY, 4);
	// in bootstrap Imgproc.cvtColor(mRgba, mHsv, Imgproc.COLOR_RGB2HSV_FULL, 4);

	// Imgproc.GaussianBlur(mHsv, mHsv, ksize, 0); // mPool src CV_8U
	// Here i'm only using the external contours and by eroding we make the draw a teeny bit faster and the result a lot
	// smoother on the rough edges where the colour fades out of range by losing a lot of the little spiky corners.
	// Imgproc.erode(mHsv, mHsv, mErodeKernel);
	// in bootstrap Imgproc.morphologyEx(mHsv, mHsv, Imgproc.MORPH_RECT, mErodeKernel);

	// Size kernel_size = new Size(3, 3);
	// Mat kernel = Imgproc.getStructuringElement(Imgproc.CV_SHAPE_RECT, kernel_size/* , point_kernel */);

	// Imgproc.dilate(mG, mMixedMat,new Mat()/* kernel */);
	// Core.compare(mG, mMixedMat, mCompMat, Core.CMP_EQ);

	// Core.multiply(mG, mCompMat, mG/*, 1 / 255.0*/);
	// mG.copyTo(mG, mCompMat);
	// Core.insertChannel(mR, mRgba, 0);
	// Core.insertChannel(mCompMat, mRgba, 1);
	// Core.insertChannel(mB, mRgba, 2);

	// final ImgPlus<T> imgPlus = cellValue.getImgPlus();
	// PickImagePeaks_2 peakPicker = new PickImagePeaks_2(mRgba);
	// peakPicker.process();
	// Img<BitType> Mat res = peakPicker.getResult();
	// ArrayList<long[]> peaks = peakPicker.getPeakList();
	// RandomAccess<BitType> access = res.randomAccess();
	// for (long[] l : peaks) {
	// access.setPosition(l);
	// access.get().set(true);
	// }
	// return m_imgCellFactory.createCell(res, imgPlus);

	// Scalar org_min_pool = min_pool.clone();
	// Scalar org_max_pool = max_pool.clone();
	// If we have selected a new point, get the color range and decide the
	// color range
	/* if (bNewColor && mSelectedPoint != null) { double[] selectedColor =
	 * mHsv.get((int) mSelectedPoint.x, (int) mSelectedPoint.y);
	 * bNewColor = false;
	 * // We check the colors in a 5x5 pixels square (Region Of Interest) and get the average from that
	 * // if (mSelectedPoint.x < 2) { mSelectedPoint.x = 2;
	 * }
	 * // else if (mSelectedPoint.x >= (720 - 2)) {
	 * mSelectedPoint.x = 720 - 2; }
	 * // if (mSelectedPoint.y < 2) { mSelectedPoint.y = 2; }
	 * // else if (mSelectedPoint.y >= (480 - 2)) { mSelectedPoint.y = 480 - 2;
	 * } //New ROI Pascal // ROI (Region Of Interest) is used to find the average value around the point we
	 * clicked.
	 * // This will reduce the risk of getting "freak" values if the pixel where we clicked has an unexpected value
	 * Rect roiRect =
	 * new Rect((int) (mSelectedPoint.x - 3), (int) (mSelectedPoint.y - 3), 7, 7); // Get the Matrix representing the ROI
	 * Mat roi =
	 * mHsv.submat(roiRect); // Calculate the mean value of the the ROI
	 * matrix // Scalar sumColor = Core.mean(roi); MatOfDouble mean = new
	 * MatOfDouble(); MatOfDouble stdev = new MatOfDouble();
	 * Core.meanStdDev(roi, mean, stdev); Scalar sMean = new
	 * Scalar(mean.toArray()[0], mean.toArray()[1], mean.toArray()[2]);
	 * Scalar sDev = new Scalar(stdev.toArray()[0], stdev.toArray()[1],
	 * stdev.toArray()[2]); Log.i(TAG, sMean + " " + sDev); // double[]
	 * sumColorValues = sumColor.val; // sumColorValues[0] =
	 * sumColorValues[0]<THRESHOLD_LOW * 3?THRESHOLD_LOW * //
	 * 3:sumColorValues[0]; // sumColorValues[1] = //
	 * sumColorValues[1]<THRESHOLD_LOW?THRESHOLD_LOW:sumColorValues[1]; //
	 * sumColorValues[2] = //
	 * sumColorValues[2]<THRESHOLD_LOW?THRESHOLD_LOW:sumColorValues[2]; //
	 * Decide on the color range based on the mean value from the ROI if
	 * (selectedColor != null) { // mLowerColorLimit = new
	 * Scalar(sumColorValues[0] - THRESHOLD_LOW * 3, // sumColorValues[1] -
	 * THRESHOLD_LOW, sumColorValues[2] // - THRESHOLD_LOW); //
	 * mUpperColorLimit = new Scalar(sumColorValues[0] + THRESHOLD_HIGH * 3,
	 * // sumColorValues[1] + THRESHOLD_HIGH, sumColorValues[2] // +
	 * THRESHOLD_HIGH); mLowerColorLimit = new Scalar(sMean.val[0] -
	 * sDev.val[0], sMean.val[1] - sDev.val[1], sMean.val[2] - sDev.val[2]);
	 * mUpperColorLimit = new Scalar(sMean.val[0] + sDev.val[0],
	 * sMean.val[1] + sDev.val[1], sMean.val[2] + sDev.val[2]); //
	 * Log.i(TAG, "New Color" + mLowerColorLimit.toString() + " " // +
	 * mUpperColorLimit.toString()); min_pool = mLowerColorLimit; max_pool =
	 * mUpperColorLimit; } else { min_pool = org_min_pool; max_pool =
	 * org_max_pool; } } */
	// in bootstrap Core.inRange(mHsv, min_pool, max_pool, mPool);
	// Core.inRange(mHsv, min_shadow, max_shadow, mShadow);
	// Core.inRange(mHsv, min_sun, max_sun, mSun);

	// http://www.morethantechnical.com/2013/03/05/skin-detection-with-probability-maps-and-elliptical-boundaries-opencv-wcode/#more-1247
	// take the pixels that are inside the ranges in both colorspaces
	// in bootstrap Core.inRange(mRgba, min_pool_rgb, max_pool_rgb, mPoolRgb);
	// combine the masks
	// in bootstrap Core.bitwise_and(mPool, mPoolRgb, mPool);
	// if (!colorProbability.isInitialized() || !contourOk) {
	 if (!contourOk) mPool = colorProbability.boostrap(mRgba);
	
	//condensation.initialise();
	// Log.i("BootStrap mPool", mPool.toString());
	// }else {
	// Log.i("bootstrap mPool ok", mPool.dump());
	// mPool = colorProbability.predict(mRgba, mPool);
	// Log.i("predict mPool", mPool.dump());
	// }
	// Log.i(TAG, mPool.toString());
//	for (int i = 0; i < 3; i++) { // predict-train 3 times for convergence
	    // Imgproc.medianBlur(mPool, mPool, 3);
//	    colorProbability.train(mRgba, mPool);
//	    Core.bitwise_and(mPool, colorProbability.predict(mRgba, mPool), mPool);
//	}

	// After we get a bootstrap, we can start training the model to get a better segmentation.
	// Essentially this means calculating 2 histograms, one for skin pixels and one for non-skin pixels.
	// colorProbability.train(mRgba, mPool);

	// mPool = colorProbability.predict(mPool);
	// Log.i("Proba", mPool.toString());
	List<MatOfPoint> contours = new ArrayList<MatOfPoint>();
	Mat hierarchy = new Mat();
	// Find contour
	Imgproc.findContours(mPool, contours, hierarchy, Imgproc.RETR_TREE,
		Imgproc.CHAIN_APPROX_TC89_KCOS);

	double maxArea = -1;
	int maxAreaIdx = -1;
	for (int idx = 0; idx < contours.size(); idx++) {
	    contour = contours.get(idx);
	    double contourarea = Imgproc.contourArea(contour);
	    if (contourarea > maxArea) {
		maxArea = contourarea;
		maxAreaIdx = idx;
	    }
	}
	if (maxAreaIdx > 0 && contours.size() > 5) {

	    // Get all the point of this contour
	    MatOfPoint2f mMOP2f = new MatOfPoint2f();
	    mMOP = contours.get(maxAreaIdx);
	    mMOP.convertTo(mMOP2f, CvType.CV_32FC2);

	    // Get the bound of eclipse
	    if (mMOP2f.total() > 5) {
		contourOk = true;
		/* List<MatOfPoint> hullPoints = new ArrayList<MatOfPoint>();
		 * List<Point> hullPointList = new ArrayList<Point>();
		 * MatOfPoint hullPointMat = new MatOfPoint();
		 * MatOfInt hull = new MatOfInt(); Imgproc.convexHull(mMOP, hull);
		 * for (int j = 0; j < hull.toList().size(); j++) {
		 * hullPointList.add(mMOP.toList().get(hull.toList().get(j))); }
		 * hullPointMat.fromList(hullPointList);
		 * hullPoints.add(hullPointMat);
		 * Imgproc.drawContours(mRgba, hullPoints, maxAreaIdx, colorYellow, 1); */
		MatOfInt hull = new MatOfInt();
		Imgproc.convexHull(mMOP, hull); // find the convex hull of
						// contour
		// Hull points' indices
		int[] intlist = hull.toArray();
		List<Point> plist = new ArrayList<Point>();
		plist.clear();
		for (int i = 0; i < intlist.length; i++) {
		    plist.add(mMOP.toList().get(hull.toList().get(i)));
		}
		// MatOfPoint2f contour2f = new MatOfPoint2f();
		// MatOfPoint2f mat_hull2f = new MatOfPoint2f();
		// MatOfPoint mat_hull_approximated = new MatOfPoint();
		MatOfPoint mat_hull = new MatOfPoint();

		List<MatOfPoint> hull_List = new ArrayList<MatOfPoint>();

		mat_hull.fromList(plist);
		hull_List.clear();
		hull_List.add(mat_hull/* _approximated */);
		// Imgproc.drawContours(mRgba, hull_List, -1, colorYellow);

		// mat_hull.convertTo(contour2f, CvType.CV_32FC2);
		// Imgproc.approxPolyDP(contour2f, mat_hull2f,
		// 0.1*Imgproc.arcLength(contour2f, true), true);
		// mat_hull2f.convertTo(mat_hull_approximated, CvType.CV_32S);
		// hull_List.clear();
		// hull_List.add(mat_hull_approximated);
		// Imgproc.drawContours(mRgba, hull_List, -1, colorYellow);

		Rect boundRect = Imgproc.boundingRect(mMOP);
		// Mon masque
		Mat magic = new Mat(mRgba.size(), CvType.CV_8UC1, colorBlack);
		Imgproc.drawContours(magic, hull_List, -1, colorWhite, 1);
		// Inutile pour le moment
		// Moments mu = Imgproc.moments(mMOP, false);
		// Point mc = new Point(mu.get_m10() / mu.get_m00(), mu.get_m01()/ mu.get_m00());
		// Core.circle(mRgba, mc, 4, colorGreen, -1, 8, 0);

		Size hullSize = boundRect.size();

		Mat hough_lines = new Mat();
		double minLineLength = 0.25 * Math.min(hullSize.width,
			hullSize.height);

		Imgproc.HoughLinesP(magic, hough_lines,
			Imgproc.CV_HOUGH_PROBABILISTIC, Math.PI / 180, 70,
			minLineLength, minLineLength * 0.10); // was 70
		if (!hough_lines.empty()) {
		    for (int x = 0; x < Math.min(hough_lines.cols(), 6); x++) {
			vecHoughLines = hough_lines.get(0, x);

			if (vecHoughLines.length == 0)
			    break;

			x1 = vecHoughLines[0];
			y1 = vecHoughLines[1];
			x2 = vecHoughLines[2];
			y2 = vecHoughLines[3];

			pt1.x = x1;
			pt1.y = y1;
			pt2.x = x2;
			pt2.y = y2;

			Core.line(mRgba, pt1, pt2, colorRed, 2);
		    }
		    vecHoughLines = null;
		}
		// double angle = atan2(y2 - y1, x2 - x1) * 180.0 / CV_PI;

		/* RotatedRect rotatedRect = Imgproc.fitEllipse(mMOP2f); //
		 * rotated rectangle Point[] rect_points = new Point[4];
		 * rotatedRect.points(rect_points); for (int j = 0; j < 4; j++)
		 * { Core.line(mRgba, rect_points[j], rect_points[(j + 1) % 4],
		 * colorGreen, 1); } */
		Core.rectangle(mRgba, boundRect.tl(), boundRect.br(), colorRed);
		// perimeter = Imgproc.arcLength(mMOP2f, true);

		// On reprend le masque pour shadow
		/* dx = GaussianFilter[ImageData[bin], 5, {0, 1}];
		 * dy = GaussianFilter[ImageData[bin], 5, {1, 0}];
		 * 
		 * edgePixels = Position[ImageData[edges], 1];
		 * normalDirection = Transpose[{Extract[dy, edgePixels], Extract[dx, edgePixels]}]; */
		// Imgproc.GaussianBlur(magic, dx, hullSize, sigmaX, sigmaY);

		// Point vp = new Point();
		// ipm.mcvGetIPMExtent();
		// vp = ipm.mcvGetVanishingPoint();
		// Log.i(TAG, vp.x + " "+ vp.y);
		// vp.x = vp.x + 360;
		// vp.y =vp.y + 480;
		// Core.line(mRgba, boundRect.tl() /*Point(0,480)*/, vp,
		// colorWhite, 2);
		// Core.line(mRgba, boundRect.br(),/*Point(720,480),*/ vp,
		// colorWhite, 2);
		// DrawCross(mRgba, colorWhite, vp);
		double ratio = (720 * 480)
			/ (boundRect.width * boundRect.height);
		ShowTitle(
			string = String.format("SURF: %2.1f %2.1f cm", ratio,
				ipm.areaTarget() * 100), 4, colorRed);
		// Rect MinRoi = Imgproc.boundingRect(mMOP);
		// Get the Matrix representing the ROI
		Mat roi = mRgba/* mHsv */.submat(boundRect);
		mPool = colorProbability.newRoi(roi);
		/* On oublie un moment ce système d'expansion
		 * MatOfDouble mean = new MatOfDouble();
		 * MatOfDouble stdev = new MatOfDouble();
		 * Core.meanStdDev(roi, mean, stdev);
		 * Scalar sMean = new Scalar(mean.toArray()[0], mean.toArray()[1],
		 * mean.toArray()[2]);
		 * Scalar sDev = new Scalar(stdev.toArray()[0],
		 * stdev.toArray()[1], stdev.toArray()[2]);
		 * // Log.i(TAG, sMean+" "+sDev);
		 * mLowerColorLimit = new Scalar(sMean.val[0] - sDev.val[0],
		 * sMean.val[1] - sDev.val[1], sMean.val[2] - sDev.val[2]);
		 * mUpperColorLimit = new Scalar(sMean.val[0] + sDev.val[0],
		 * sMean.val[1] + sDev.val[1], sMean.val[2] + sDev.val[2]);
		 * // Log.i(TAG, "New Color" + mLowerColorLimit.toString() + " "
		 * // + mUpperColorLimit.toString());
		 * 
		 * min_pool = mLowerColorLimit;
		 * max_pool = mUpperColorLimit;
		 * 
		 * // min_pool = new Scalar(113, 52, 104, 255); // Piscine
		 * max_pool = new Scalar(146, 157, 254, 255);
		 * if (min_pool.val[0] < org_min_pool.val[0] - 1. * THRESHOLD_LOW)
		 * min_pool.val[0] = org_min_pool.val[0];
		 * if (max_pool.val[0] > org_max_pool.val[0] + 1. * THRESHOLD_HIGH)
		 * min_pool.val[0] = org_max_pool.val[0];
		 * if (min_pool.val[2] < org_min_pool.val[2] - 1. * THRESHOLD_LOW)
		 * min_pool.val[2] = org_min_pool.val[2];
		 * if (max_pool.val[2] > org_max_pool.val[2] + 1. * THRESHOLD_HIGH)
		 * min_pool.val[2] = org_max_pool.val[2]; */
	    } else
		contourOk = false;
	}

	// Core.ellipse(mGray, rotatedRect, colorRed);
	// Core.rectangle( mGray, boundRect.tl(), boundRect.br(), colorGreen, 2,
	// 8, 0 );

	// Imgproc.drawContours(mRgba, contours, maxAreaIdx, colorRed, 1);

	// distance(mm) = (focal length (mm) * real height of the object (mm) * camera frame height in device (pixels) )/
	// ( image height (pixels) * sensor height (mm))
	// where Sensorheight for perticular device can be obtain from
	// http://cameraimagesensor.com/
	// for phone and tablet camera it is 3.42mm

	// The intrinsic matrix contains the focal lengths fx and fy. These are
	// the
	// distances, measured in pixels, from the camera's principal point to
	// the
	// image plane. If you also have the image size, you can get at the FOV
	// with
	// some trigonometry:

	// theta_x/2 = tan_inv( (width/2) / fx )
	// theta_y/2 = tan_inv( (height/2) / fy )

	// Shadow Part
	/* One characteristic of shadow is your intensity (channel V) is low. It
	 * lie in a Range between lower and high Value, you will need discover
	 * your. Other characteristic is the saturation is lower. Frequently the
	 * object (caster) have high saturation. */
	// tan(sun) = Hauteur mur / longueur ombre
	/* hierarchy = new Mat(); // Find contour Imgproc.findContours(mShadow,
	 * contours, hierarchy, Imgproc.RETR_TREE,
	 * Imgproc.CHAIN_APPROX_SIMPLE);
	 * maxArea = -1; maxAreaIdx = -1; for (int idx = 0; idx <
	 * contours.size(); idx++) { Mat contour = contours.get(idx); double
	 * contourarea = Imgproc.contourArea(contour); if (contourarea >
	 * maxArea) { maxArea = contourarea; maxAreaIdx = idx; } } */
	/* if (maxAreaIdx > 0 && contours.size() > 5) { // Get all the point of
	 * this contour MatOfPoint2f mMOP2f = new MatOfPoint2f(); mMOP =
	 * contours.get(maxAreaIdx); mMOP.convertTo(mMOP2f, CvType.CV_32FC2); if
	 * (mMOP2f.total() > 5) { RotatedRect rotatedRect =
	 * Imgproc.fitEllipse(mMOP2f); // rotated rectangle Point[] rect_points
	 * = new Point[4]; rotatedRect.points(rect_points); for (int j = 0; j <
	 * 4; j++) Core.line(mRgba, rect_points[j], rect_points[(j + 1) % 4],
	 * colorGreen, 1); } } */
	// Imgproc.drawContours(mRgba, contours, maxAreaIdx, colorRed, 1);

	/* threshold = Imgproc .threshold(mPool, mPool, 0, 500,
	 * Imgproc.THRESH_OTSU); Imgproc.Canny(mPool, img_edge, threshold *
	 * 0.4,
	 * threshold); Mat kernel =
	 * Imgproc.getStructuringElement(Imgproc.MORPH_ELLIPSE, ksize);
	 * Imgproc.morphologyEx(img_edge, img_edge, Imgproc.MORPH_CLOSE,
	 * kernel); kernel.release(); int minLineSize = 30; int lineGap = 4; */
	// Imgproc.HoughLinesP(img_edge, lines, rho = 1, theta = Math.PI / 180,
	// 100, minLineSize, lineGap);

	// Filter the lines we are not interested based on their angle
	// Extract the angles from the remaining lines
	// Calculate the vanish points relative to horizon measured from sensor
	// Identify the most probably farthest point from all the samples

	// Appartenance d'un point à un contour :
	// dist = cv2.pointPolygonTest(cnt,(50,50),True);
	// Min enclosing circle CvType.CV_32S
	/* Imgproc.minEnclosingCircle(mMOP2f, center, radius);
	 * Core.circle(mGray, center, Integer.parseInt(radius.toString()) ,
	 * colorRed, 2); */

	/* def do_green_line_detection(img, framenum=0): */

	/* The time-to-contact the FOE - focus of expansion - is the
	 * intersection point of all the pixel flows. In the more general case,
	 * the FoE is the epipole which can be found once you have estimated the
	 * Fundamental matrix. OpenCV allows you to find the fundamental matrix
	 * and the epipoles. One way to compute the TTC wold be like this: -
	 * compute the distance from the focus of expansion to the center of the
	 * image: D_foe = abs(x_cim - x_foe) + abs(y_cim - y_foe); - compute the
	 * TTC by substituting the value 100 with D_foe: ttc =
	 * D_foe/sum(Vmag(:)); */

	/* Imgproc.cvtColor(mRgba, mGray, Imgproc.COLOR_YUV420sp2RGB, 4); //
	 * ***FEATURE FINDING*** if (MatOfPointcorners == null) {
	 * MatOfPointcorners = new MatOfPoint(); } // MatOfPointcorners.clear();
	 * Imgproc.cvtColor(mGray, mGray, Imgproc.COLOR_RGB2GRAY);
	 * Imgproc.goodFeaturesToTrack(mGray, MatOfPointcorners, 20, 0.1, 5);
	 * for (Point p : MatOfPointcorners.toArray()) { Core.circle(mRgba, p,
	 * 8, circleCol, 1); } */
	// Convert and send the image.
	// Imgproc.cvtColor(mRgba, mGray, Imgproc.COLOR_YUV420sp2RGB,
	// 4);

	// notifyResultCallback(mGray);

	// Utils.matToBitmap(mGray, mBitmap, true);
	double north = 0, pitch = 0, roll = 0;
	Point perim = new Point();
	double rayon = Math.hypot(480, 720) / 2;
	double sun = 0;
	double sunx = 0;
	double suny = 0;
	double Temp = 0;

	try { // from : Orientation = new Scalar(mAzimuth, mPitch, mRoll);
	    north = Orientation.val[0];
	    pitch = Orientation.val[1];
	    roll = Orientation.val[2];

	} catch (Exception e) {
	    north = 0;
	    pitch = 0;
	    roll = 0;
	}

	try {
	    sunx = SunPos.x;
	    suny = SunPos.y;
	} catch (Exception e) {
	    sunx = 0;
	    suny = 0;
	}
	try {
	    Temp = getTemp();
	    // Log.i(TAG, "Temp Avg: " + Temp);
	    ShowTitle(String.format("TEMP: %2.1f", Temp), 8, colorRed);
	} catch (Exception e) {
	    Temp = 0.;
	} finally {
	    // Point perim = new Point();
	    // double rayon = Math.hypot(480, 720) / 2;
	    north = 270 + (Math.abs(north - 270));
	    sun = north + sunx;
	    // Fond out the co-ordinates in the circle perimeter for north and
	    // draw the line from center
	    perim.x = center.x + rayon * Math.cos(Math.toRadians(north));
	    perim.y = center.y + rayon * Math.sin(Math.toRadians(north));
	    Core.line(mRgba, center, perim, colorYellow);
	    ShowTitle(
		    string = String.format("HEAD: %2.1f %2.1f %2.1f",
			    north % 360, pitch, roll), 1, colorRed);

	    perim.x = center.x + rayon * Math.cos(Math.toRadians(sun));
	    perim.y = center.y + rayon * Math.sin(Math.toRadians(sun));
	    Core.line(mRgba, center, perim, colorRed);
	    ShowTitle(string = String.format("SUN: %2.1f %2.1f", sunx, suny),
		    2, colorRed);

	    cameraParam.pitch = -roll * Math.PI / 180; // Inversé à cause de
						       // landscape ?!
	    cameraParam.yaw = pitch * Math.PI / 180;
	}

	// Correction du roulis ...
	// Mat r = Imgproc.getRotationMatrix2D(center, pitch, 1.0);
	// Imgproc.warpAffine(mRgba, mRgba, r, mRgba.size());

	long lMilliShotTime = 0;
	if (bShootNow) {
	    // get the time of the attempt to save a screenshot
	    lMilliShotTime = System.currentTimeMillis();
	    bShootNow = false;

	    // try it, and set the screen text accordingly.
	    // this text is shown at the end of each frame until
	    // 1.5 seconds has elapsed
	    if (SaveImage(mRgba)) {
		sShotText = "SCREENSHOT SAVED";
	    } else {
		sShotText = "SCREENSHOT FAILED";
	    }

	}

	if (System.currentTimeMillis() - lMilliShotTime < 1500)
	    ShowTitle(sShotText, 3, colorRed);

	return mRgba;
    }

    /* SAA, SEC */
    public static void setSun(Point pSun) {
	SunPos = pSun;
    }

    public static Point getSun() {
	return SunPos;
    }

    public Location getLocation() {
	return location;
    }

    /** mAzimuth, mPitch, mRoll */
    public static void setOrientation(Scalar sOr) {
	Orientation = sOr;
	// mOpenCvCameraView.setCameraDistance(distance * dm.density);
	// mOpenCvCameraView.setRotationX((float) sOr.val[1]);
	// mOpenCvCameraView.setRotationY((float) sOr.val[2]);
    }

    /** Is called each temp received (36) */
    public void setWeather(HourlyForecast forecast) {

	ForecastHour forecastHour;

	try {
	    for (int i = 1; i < forecast.getSize(); i++) {
		forecastHour = forecast.getForecastHours()[i];
		this.Temp = 0;
		try {
		    this.Temp = forecastHour.getTemp();
		} catch (Exception e) {
		    // TODO Auto-generated catch block
		    // e.printStackTrace();
		}
		Log.e(TAG, this.Temp + "");
		if (this.Temp != 0)
		    temperature.add(this.Temp);
		// Temp += Temp;
		// Temp /= i;
		// Log.w(TAG, " " +forecastHour.getTemp());
	    }
	} catch (Exception e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	}
	Log.i(TAG, "Got Forecast" + temperature.size());
    }

    public double getTemp() {
	double Average = 0;
	for (double temp : temperature)
	    Average += temp;
	this.Temp = Average / temperature.size();
	return Temp;
    }

    @SuppressLint("SimpleDateFormat")
    public boolean SaveImage(Mat mat) {

	Imgproc.cvtColor(mat, mIntermediateMat, Imgproc.COLOR_RGBA2BGR, 3);

	File path = Environment
		.getExternalStoragePublicDirectory(Environment.DIRECTORY_PICTURES);

	filename = "PoolChecker_";
	SimpleDateFormat fmt = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
	Date date = new Date(System.currentTimeMillis());
	String dateString = fmt.format(date);
	filename += dateString + "-" + iFileOrdinal;
	filename += ".jpg";

	File file = new File(path, filename);

	Boolean bool = null;
	filename = file.toString();
	MatOfInt jpgParam = new MatOfInt();
	jpgParam.put(0, 1, Highgui.CV_IMWRITE_JPEG_QUALITY);
	jpgParam.put(0, 2, 95);
	bool = Highgui.imwrite(filename, mIntermediateMat);

	Log.i(TAG, filename);
	// http://stackoverflow.com/questions/2170214/image-saved-to-sdcard-doesnt-appear-in-androids-gallery-app
	sendBroadcast(new Intent(
		Intent.ACTION_MEDIA_MOUNTED,
		Uri.parse("file://"
			+ Environment
				.getExternalStoragePublicDirectory(Environment.DIRECTORY_PICTURES))));

	// Upload to server task
	new Thread(new Runnable() {
	    @Override
	    public void run() {
		runOnUiThread(new Runnable() {
		    @Override
		    public void run() {
		    }
		});
		// int response = uploadFile(filename);
		// System.out.println("RES : " + response);
	    }
	}).start();

	// Log.e(TAG, mOpenCvCameraView.getMatrix().toString());
	// if (bool == false)
	// Log.d("Baz", "Fail writing image to external storage");

	return bool;
    }

    public void ShowTitle(String s, int iLineNum, Scalar color) {
	Core.putText(mRgba, s, new Point(10,
		(int) (dTextScaleFactor * 60 * iLineNum)),
		Core.FONT_HERSHEY_SIMPLEX, dTextScaleFactor, color, 2);
    }

    public void setLocation(Location location) {
	// TODO Auto-generated method stub
	this.location = location;
    }

    public int uploadFile(String sourceFileUri) {
	String upLoadServerUri = "http://priveyes.craym.eu/pool_checker/upload_media_test.php";
	String fileName = sourceFileUri;

	HttpURLConnection conn = null;
	DataOutputStream dos = null;
	String lineEnd = "\r\n";
	String twoHyphens = "--";
	String boundary = "*****";
	int bytesRead, bytesAvailable, bufferSize;
	byte[] buffer;
	int maxBufferSize = 1 * 1024 * 1024;
	File sourceFile = new File(sourceFileUri);
	if (!sourceFile.isFile()) {
	    Log.e("uploadFile", "Source File Does not exist");
	    return 0;
	}
	int serverResponseCode = 0;
	try { // open a URL connection to the <span id="IL_AD8"
	      // class="IL_AD">Servlet</span>
	    FileInputStream fileInputStream = new FileInputStream(sourceFile);
	    URL url = new URL(upLoadServerUri);
	    conn = (HttpURLConnection) url.openConnection(); // Open a HTTP
							     // connection to
							     // the URL
	    conn.setDoInput(true); // Allow Inputs
	    conn.setDoOutput(true); // Allow Outputs
	    conn.setUseCaches(false); // Don't use a Cached Copy
	    conn.setRequestMethod("POST");
	    conn.setRequestProperty("Connection", "Keep-Alive");
	    conn.setRequestProperty("ENCTYPE", "multipart/form-data");
	    conn.setRequestProperty("Content-Type",
		    "multipart/form-data;boundary=" + boundary);
	    conn.setRequestProperty("uploaded_file", fileName);
	    dos = new DataOutputStream(conn.getOutputStream());

	    dos.writeBytes(twoHyphens + boundary + lineEnd);
	    dos.writeBytes("Content-Disposition: form-data; name=\"uploaded_file\";filename=\""
		    + fileName + "\"" + lineEnd);
	    dos.writeBytes(lineEnd);

	    bytesAvailable = fileInputStream.available(); // create a buffer of
							  // maximum
							  // size

	    bufferSize = Math.min(bytesAvailable, maxBufferSize);
	    buffer = new byte[bufferSize];

	    // read file and write it into form...
	    bytesRead = fileInputStream.read(buffer, 0, bufferSize);

	    while (bytesRead > 0) {
		dos.write(buffer, 0, bufferSize);
		bytesAvailable = fileInputStream.available();
		bufferSize = Math.min(bytesAvailable, maxBufferSize);
		bytesRead = fileInputStream.read(buffer, 0, bufferSize);
	    }

	    // send multipart form data necesssary after file data...
	    dos.writeBytes(lineEnd);
	    dos.writeBytes(twoHyphens + boundary + twoHyphens + lineEnd);

	    // Responses from the server (code and message)
	    serverResponseCode = conn.getResponseCode();
	    String serverResponseMessage = conn.getResponseMessage();

	    Log.i("uploadFile", "HTTP Response is : " + serverResponseMessage
		    + ": " + serverResponseCode);
	    if (serverResponseCode == 200) {
		runOnUiThread(new Runnable() {
		    @Override
		    public void run() {
			ShowTitle("File Upload Completed.", 3, colorRed);
			// tv.setText("File Upload Completed.");
			// Toast.makeText(UploadImageDemo.this,
			// "File Upload Complete.",
			// Toast.LENGTH_SHORT).show();
		    }
		});
	    }

	    // close the streams //
	    fileInputStream.close();
	    dos.flush();
	    dos.close();

	} catch (MalformedURLException ex) {
	    // dialog.dismiss();
	    ex.printStackTrace();
	    // Toast.makeText(UploadImageDemo.this, "MalformedURLException",
	    // Toast.LENGTH_SHORT).show();
	    Log.e("Upload file to server", "error: " + ex.getMessage(), ex);
	} catch (Exception e) {
	    // dialog.dismiss();
	    e.printStackTrace();
	    // Toast.makeText(UploadImageDemo.this, "Exception : " +
	    // e.getMessage(),
	    // Toast.LENGTH_SHORT).show();
	    Log.e("Upload file to server Exception",
		    "Exception : " + e.getMessage(), e);
	}
	// dialog.dismiss();
	return serverResponseCode;
    }

    public void DrawCross(Mat mat, Scalar color, Point pt) {
	int iCentreCrossWidth = 24;
	pt1 = new Point();
	pt2 = new Point();

	pt1.x = pt.x - (iCentreCrossWidth >> 1);
	pt1.y = pt.y;
	pt2.x = pt.x + (iCentreCrossWidth >> 1);
	pt2.y = pt.y;

	Core.line(mat, pt1, pt2, color, 1);

	pt1.x = pt.x;
	pt1.y = pt.y + (iCentreCrossWidth >> 1);
	pt2.x = pt.x;
	pt2.y = pt.y - (iCentreCrossWidth >> 1);

	Core.line(mat, pt1, pt2, color, 1);

    }
}
