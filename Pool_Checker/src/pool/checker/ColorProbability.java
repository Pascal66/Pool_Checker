package pool.checker;

import java.util.ArrayList;
import java.util.List;

import org.opencv.core.Core;
import org.opencv.core.Core.MinMaxLocResult;
import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.core.MatOfFloat;
import org.opencv.core.MatOfInt;
import org.opencv.core.Scalar;
import org.opencv.core.Size;
import org.opencv.imgproc.Imgproc;

import android.util.Log;

// package pool.checker;
/* SkinProbablilityMaps.h
 * CurveMatching
 * 
 * Created by roy_shilkrot on 2/8/13.
 * Copyright (c) 2013 MIT
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 * ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH
 * THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 * 
 * based on section 2.2.2 of "A survey of skin-color modeling and detection methods", Kukumanu et al 2007 */

// #include "AbstractSkinDetector.h"

public class ColorProbability {
    // private:
    public Mat pool_hist;
    public Mat non_skin_hist;

    float theta_thresh;
    MatOfInt hist_bins;
    Scalar low_range;
    Scalar high_range;
    Scalar range_dist;

    private Scalar max_pool, min_pool_rgb;
    private Scalar min_pool, max_pool_rgb;

    public MatOfFloat ranges;
    public MatOfInt channels, histSize;
    public boolean initialized = false;
    
    Mat hsv;
    
    Size ksize;
    Mat mErodeKernel;

    /** constructor ... */
    public ColorProbability() {
	this.pool_hist = new Mat();
	this.non_skin_hist = new Mat();

	setTheta(8);
	hist_bins = new MatOfInt(50, 50); // new Scalar(50, 50);
	low_range = new Scalar(0.2 * 255, 0.3 * 255);
	high_range = new Scalar(0.4 * 255, 0.5 * 255);

	range_dist = new Scalar(high_range.val[0] - low_range.val[0],
		high_range.val[1] - low_range.val[1]);
	// range_dist.val[0] = high_range.val[0] - low_range.val[0];
	// range_dist.val[1] = high_range.val[1] - low_range.val[1];

	channels = new MatOfInt();
	histSize = new MatOfInt();
	ranges = new MatOfFloat(0f, 256f);

	min_pool = new Scalar(113, 52, 104, 255); // Piscine 113- 52-104
	max_pool = new Scalar(146, 157, 254, 255); // 146-157-254

	min_pool_rgb = new Scalar(0, 33, 55, 255);
	max_pool_rgb = new Scalar(150, 205, 231, 255); // R55 G110 B130 Shadow
	
	ksize = new Size(3, 3);
	
	mErodeKernel = Imgproc
		.getStructuringElement(Imgproc.MORPH_CROSS, ksize);
	
	hsv = new Mat();
    }

    /* simple threshold on HSV and normalized-RGB taken from:
     * "A survey of skin-color modeling and detection methods", P. Kakumanu, S. Makrogiannis, N. Bourbakis 2006 */
    /**
     * For bostratpping (namely, getting an initial guess for skin pixels) I used the method from Gomez & Morales 2002:
     * 
     * @param rgb
     * @return mask from initial hsv and rgb range for the pool
     */
    public Mat boostrap(Mat rgb) {
	
	Imgproc.cvtColor(rgb, hsv, Imgproc.COLOR_RGB2HSV_FULL/* COLOR_BGR2HSV */);

	// Mat nrgb = getNormalizedRGB(rgb);

	Mat mask_hsv = new Mat();
	Mat mask_nrgb = new Mat();

	Core.inRange(hsv, min_pool, max_pool, mask_hsv);
	Core.inRange(rgb, min_pool_rgb, max_pool_rgb, mask_nrgb);
	// org H=[0,50], S= [0.20,0.68] and V= [0.35,1.0]
	// org Core.inRange(hsv, new Scalar(0, 0.2 * 255.0, 0.35 * 255.0), new Scalar( 50.0 / 2.0, 0.68 * 255.0, 1.0 * 255.0), mask_hsv);
	// org r = [0.36,0.465], g = [0.28,0.363]
	// org Core.inRange(nrgb, new Scalar(0, 0.28, 0.36), new Scalar(1.0, 0.363, 0.465), mask_nrgb);

	// rule from "Automatic Feature Construction and a Simple Rule Induction Algorithm for Skin Detection", Gomez & Morales 2002
	// //r/g > 1.185, rb/(r+g+b)^2 > 0.107, rg/(r+g+b)^2 > 0.112
	// inRange(nrgb, Scalar(1.185,0.107,0.112), Scalar::all(1000), mask_nrgb);

	Mat outputmask = new Mat();
	// outputmask = mask_hsv & mask_nrgb;
	Core.bitwise_and(mask_hsv, mask_nrgb, outputmask);
	mask_hsv.release();
	mask_nrgb.release();
	
	Imgproc.morphologyEx(outputmask, outputmask, Imgproc.MORPH_RECT, mErodeKernel);
	return outputmask;
    }

    public Mat newRoi(Mat img_roi) {
	ArrayList<Mat> list = new ArrayList<Mat>();
	list.add(img_roi);
	MatOfFloat ranges = new MatOfFloat(0f, 256f);
	MatOfInt histSize = new MatOfInt(25);
	channels = new MatOfInt(0); // new MatOfInt(1, 2);
	Imgproc.calcHist(list, channels, new Mat(0), pool_hist, histSize, ranges);
	Core.normalize(pool_hist, pool_hist);
	Imgproc.calcBackProject(list, channels, pool_hist, img_roi, ranges, 255);
	return img_roi;
    }

    public Mat predict_2(Mat img_rgb, Mat mask) {
	ArrayList<Mat> list = new ArrayList<Mat>();
	list.add(img_rgb);
	MatOfFloat ranges = new MatOfFloat(0f, 256f);
	MatOfInt histSize = new MatOfInt(25);
	channels = new MatOfInt(0); // new MatOfInt(1, 2);
	Imgproc.calcBackProject(list, channels, pool_hist, mask, ranges, 255);

	return channels;

    }

    /**
     * After we get a bootstrap, which you can see is pretty bad, we can start training the model to get a better segmentation. Essentially
     * this means calculating 2 histograms, one for skin pixels and one for non-skin pixels.
     * 
     * @param img_rgb
     * @param mask
     */
    public void train(Mat img_rgb, Mat mask) {
	// On fait le camshift une seule fois pour localiser l'objet et calculer son histogramme

	// Mat nrgb = getNormalizedRGB(img_rgb);
	/*Mat notmask = mask.clone();
	 * Core.bitwise_not(mask, notmask);
	 * 
	 * // pool_hist.setTo(0); non_skin_hist.setTo(0);
	 * // les hist sont à 250*250*32FC1
	 * this.skin_hist = calc_rg_hist(img_rgb, mask, hist_bins, low_range,
	 * high_range);
	 * this.non_skin_hist = calc_rg_hist(img_rgb, notmask, hist_bins,
	 * low_range, high_range);
	 * // Log.i("Prob ", pool_hist.toString());
	 * // create a probabilty density function
	 * Scalar skin_pixels = new Scalar(Core.countNonZero(mask)); // Min 0 max 345600
	 * Scalar non_skin_pixels = new Scalar(Core.countNonZero(notmask));
	 * Core.divide(this.skin_hist, skin_pixels, this.skin_hist);
	 * Core.divide(this.non_skin_hist, non_skin_pixels, this.non_skin_hist);
	 * /* for (int ubin = 0; ubin < hist_bins.toArray()[0]; ubin++) {
	 * for (int vbin = 0; vbin < hist_bins.toArray()[1]; vbin++) {
	 * this.skin_hist.get(ubin, vbin, value);
	 * if (value[0] > 0 || skin_pixels != 0) {
	 * this.skin_hist.put(ubin, vbin, value[0] / skin_pixels);
	 * }
	 * // if (pool_hist.at<float>(ubin,vbin) > 0) {pool_hist.at<float>(ubin,vbin) /= skin_pixels; }
	 * this.non_skin_hist.get(ubin, vbin, value);
	 * if (value[0] > 0 || non_skin_pixels != 0) {
	 * this.non_skin_hist.put(ubin, vbin, value[0]
	 * / non_skin_pixels);
	 * }
	 * // if (non_skin_hist.at<float>(ubin,vbin) > 0) {non_skin_hist.at<float>(ubin,vbin) /= non_skin_pixels;}
	 * }
	 * }
	 * // Log.i("Prob ", non_skin_hist.dump()); */
	// notmask.release();
	this.initialized = true;
    }

    /**
     * The SPM algorithm simply states you make a binary prediction (yes or no) for each pixel based on a Theta value (a variable given to
     * the
     * algorithm):
     * p(c|skin)p(c|non-skin)>Theta
     * The implementation is trivial if you already have the normalized histograms (and thus PDFs):
     * it good to re-train the model a number of time. What I mean is that I iterate a number of times with train-predict operations, where
     * every iteration I use the result of the prediction to train the model again:
     * 
     * @param img_rgb
     * @param mask
     * @return predicted mask
     */
    public Mat predict(Mat img_rgb, Mat mask) {
	if (Core.sumElems(mask).val[0] == 0)
	    return mask;
	/* http://answers.opencv.org/question/14978/how-to-editreading-pixels-from-matvec3b-using-java/ */
	// assume, that Mat with three channels
	// byte[] sourceBuffer = new byte[(int) (img_rgb.total() * img_rgb.channels())];
	Mat result_mask = mask.clone(); // new byte[(int) (mask.total() * mask.channels())];

	ArrayList<Mat> list = new ArrayList<Mat>();
	list.add(img_rgb);

	// img_rgb.get(0, 0, sourceBuffer);
	// mask.get(0, 0, result_mask);
	// float[] skin_hist_val = new float[1];
	// float[] non_skin_hist_val = new float[1];
	// int cols = img_rgb.cols();
	// int rows = img_rgb.rows();
	// for (int i = 0; i < rows; i++) {
	// for (int x = 0; x < cols * 3; x += 3) {
	// // for (int channelId = 0; channelId < 3; channelId++) {
	// byte chan1 /* nrgb(i)[1] */= sourceBuffer[i * cols + x + 1/* channelId */];
	// byte chan2 = sourceBuffer[i * cols + x + 2/* channelId */];
	// if (chan1 < low_range.val[0] || chan1 > high_range.val[0]
	// || chan2 < low_range.val[1]
	// || chan2 > high_range.val[1]) {
	// result_mask[i] = 0;
	// continue;
	// }
	// int gbin = (int) Math.round((chan1 - low_range.val[0])
	// / range_dist.val[0] * hist_bins.toArray()[0]);
	// int rbin = (int) Math.round((chan1 - low_range.val[1])
	// / range_dist.val[1] * hist_bins.toArray()[1]);

	Mat result_hist = new Mat();
	Core.divide(pool_hist, non_skin_hist, result_hist);
	Imgproc.threshold(result_hist, result_hist, theta_thresh, 255,
		Imgproc.THRESH_BINARY);
	Core.normalize(result_hist, result_hist);

	Imgproc.calcBackProject(list, channels, result_hist, result_mask,
		ranges, 1);
	// Imgproc.equalizeHist(result_hist, result_hist);
	// pool_hist.get(gbin, rbin, skin_hist_val);
	// if (skin_hist_val[0] > 0) {
	// non_skin_hist.get(gbin, rbin, non_skin_hist_val);

	// if (non_skin_hist_val[0] > 0) {
	// if ((skin_hist_val[0] / non_skin_hist_val[0]) > theta_thresh)
	// result_mask[i] = (byte) 255;
	// else
	// result_mask[i] = 0;
	// } else {
	// result_mask[i] = 0;
	// }
	// } else {
	// result_mask[i] = 0;
	// }

	// }
	// }
	// }
	// mask.put(0, 0, result_mask);
	// Mat_<Vec3f> nrgb = getNormalizedRGB(img_rgb).reshape(3, img_rgb.rows*img_rgb.cols);
	// Mat_<uchar> result_mask(nrgb.size());
	// for (int i=0; i<nrgb.rows; i++) {
	// if (nrgb.at<Vec3f>(i)[1] < low_range[0] || nrgb(i)[1] > high_range[0] ||
	// nrgb(i)[2] < low_range[1] || nrgb(i)[2] > high_range[1])
	// {
	// result_mask(i) = 0;
	// continue;
	// }
	// int gbin = cvRound((nrgb(i)[1] - low_range[0])/range_dist[0] * hist_bins[0]);
	// int rbin = cvRound((nrgb(i)[2] - low_range[1])/range_dist[1] * hist_bins[1]);
	// float skin_hist_val = pool_hist.at<float>(gbin,rbin);
	// if (skin_hist_val > 0) {
	// float non_skin_hist_val = non_skin_hist.at<float>(gbin,rbin);
	// if (non_skin_hist_val > 0) {
	// if((skin_hist_val / non_skin_hist_val) > theta_thresh)
	// result_mask(i) = 255;
	// else
	// result_mask(i) = 0;
	// } else {
	// result_mask(i) = 0;
	// }
	// } else {
	// result_mask(i) = 0;
	// }
	// }
	// output_mask = result_mask.reshape(1, img_rgb.rows);
	// }
	return result_mask;
    }

    private void setTheta(float t) {
	this.theta_thresh = t;
    }

    private Mat calc_rg_hist(Mat img, Mat mask, MatOfInt histbins, Scalar low,
	    Scalar high) {
	// hist_bins = new MatOfInt(250, 250); // new Scalar(250, 250);
	// low = new Scalar(0, 0);
	// high = new Scalar(1, 1);

	channels = new MatOfInt(1, 2); // new Scalar(1, 2);

	return calc_2D_hist(img, mask, channels, hist_bins, low, high);
    }

    /* AsbtractSkinDetector.h
     * CurveMatching
     * 
     * Created by roy_shilkrot on 2/9/13.
     * Copyright (c) 2013 MIT */

    private Mat calc_2D_hist(Mat img, Mat mask, MatOfInt wchannels,
	    MatOfInt hist_bins, Scalar low, Scalar high) {
	Mat hist = new Mat();

	// org int histSize[] = { bins.val[0], bins.val[1] };

	// org uranges[] = { low.val[0], high.val[0] };
	// org vranges[] = { low.val[1], high.val[1] };
	// org ranges[] = { uranges, vranges };

	// org int channels[] = {wchannels.val[0], wchannels.val[1]};

	/* calcHist( img, 1, channels, mask,
	 * hist, 2, histSize, ranges,
	 * true, // the histogram is uniform
	 * false ); */

	ArrayList<Mat> list = new ArrayList<Mat>();
	list.add(img);
	// MatOfInt channels = new MatOfInt(1, 2);
	// MatOfInt bins = new MatOfInt(50, 50);
	this.ranges = new MatOfFloat(0.0f, 255.0f, 0.0f, 255.0f);

	Imgproc.calcHist(list, wchannels, mask, hist, hist_bins, ranges, true);

	return hist;
    }

    private Mat getNormalizedRGB(Mat rgb) {
	assert (rgb.type() == CvType.CV_8UC3);
	Mat rgb32f = new Mat();
	rgb.convertTo(rgb32f, CvType.CV_32FC3);
	List<Mat> split_rgb = new ArrayList<Mat>(3);
	// org vector<Mat> split_rgb;

	Core.split(rgb32f, split_rgb);
	Mat sum_rgb = new Mat();
	Mat divide_rgb = new Mat();
	divide_rgb.setTo(new Scalar(0));

	Core.add(split_rgb.get(0), split_rgb.get(1), sum_rgb);
	Core.add(sum_rgb, split_rgb.get(2), sum_rgb);
	// org Mat sum_rgb = split_rgb[0] + split_rgb[1] + split_rgb[2];

	split_rgb.set(0, divide_rgb);
	// org split_rgb[0].setTo(0); //org split_rgb[0] / sum_rgb;

	Core.divide(split_rgb.get(1), sum_rgb, divide_rgb);
	split_rgb.set(1, divide_rgb);
	Core.divide(split_rgb.get(2), sum_rgb, divide_rgb);
	split_rgb.set(2, divide_rgb);
	// org split_rgb[1] = split_rgb[1] / sum_rgb;
	// org split_rgb[2] = split_rgb[2] / sum_rgb;
	// sum_rgb = sum_rgb.mul(sum_rgb);
	// split_rgb[0] = split_rgb[2] / split_rgb[1];
	// split_rgb[1] = split_rgb[2].mul(split_rgb[0]) / sum_rgb;
	// split_rgb[2] = split_rgb[2].mul(split_rgb[1]) / sum_rgb;

	Core.merge(split_rgb, rgb32f);
	return rgb32f;
    }

    public boolean isInitialized() {
	// TODO Auto-generated method stub
	return initialized;
    }

    // PP

    // https://code.google.com/p/js-handtracking/source/browse/trunk/src/cv.js
    public Mat NewDilate(Mat imageSrc) {
	// Mat imageDst = new Mat();
	Mat imageDst = imageSrc.clone();
	// imageSrc.convertTo(imageSrc, CvType.CV_64FC1);
	int size = (int) (imageSrc.total() * imageSrc.channels());
	byte[] src = new byte[size]; // use double[] instead of byte[]
	byte[] dst = new byte[size];
	imageSrc.get(0, 0, src);

	int width = imageSrc.width();
	int height = imageSrc.height();
	int[] offsets = { -width - 1, -width, -width + 1, -1, 1, width - 1,
		width, width + 1 };
	int klen = offsets.length;
	int pos = 0;
	int value;

	for (int i = 0; i < width; ++i) {
	    dst[pos++] = 0;
	}

	for (int i = 2; i < height; ++i) {
	    dst[pos++] = 0;

	    for (int j = 2; j < width; ++j) {
		value = src[pos];

		for (int k = 0; k < klen; ++k) {
		    // PickImagePeaks function :
		    /* Create an intermediate image. This image will contain a sort of signum operation of the difference
		     * along a given dimension of the input image. "Sort of" because 1 corresponds to greater than or
		     * equal to zero, while 0 corresponds to less than 0, rather than the traditions signum.
		     * I've written this method in this way in order that we don't have to care what order the
		     * cursor traverses the Img.
		     * imImagePushCursor.get().set(tc.compareTo(t0) >= 0); */// int set = 1;
									     // if (value > src[pos +
									     // offsets[k]])/*(peakVal.compareTo(inPeak.peakVal) == 1)*/ {
									     // set = -1;/* return -1;*/
									     // } else if (value == src[pos +
									     // offsets[k]])/*(peakVal.compareTo(inPeak.peakVal) == 0)*/ {
									     // set = 0; /*return 0;*/
									     // } else { set = 1; /* return 1;*/
									     // }
		    value = 0;
		    // if (set >= 0) value = 255;
		    if (value <= src[pos + offsets[k]])
			value = 255;
		    // Original dilate function : max
		    // value = Math.max(value, src[pos + offsets[k]]);
		}

		dst[pos++] = (byte) value;
	    }

	    dst[pos++] = 0;
	}

	for (int i = 0; i < width; ++i) {
	    dst[pos++] = 0;
	}

	imageDst.put(0, 0, dst);

	return imageDst;
    }

    // Implement simple white balance
    // http://scien.stanford.edu/pages/labsite/2010/psych221/projects/2010/JasonSu/simplestcb.html
    // ou miexu expliqué : http://www.ipol.im/pub/art/2011/llmps-scb/
    /**
     * simplestColorBalance(filename,outFile,satLevel,plot)
     * Performs color balancing via histogram normalization.
     * satLevel controls the percentage of pixels to clip to white and black.
     * Set plot = 0 or 1 to turn diagnostic plots on or off.
     * 
     * @return White balanced Rgb
     */

    public Mat simplestColorBalance(Mat imRGB, double satLevel) {
	imRGB.convertTo(imRGB, CvType.CV_32SC3);

	// full width histogram method
	// satLevel = .01; // percentage of the image to saturate to black or white, tweakable param
	// double q = [satLevel/2 1-satLevel/2];

	// org vector<Mat> split_rgb;

	// imRGB_orig = cbreshape(im_orig)*255;
	// imRGB = zeros(size(imRGB_orig));
	// N = size(imRGB_orig,2);
	// color = {'r','g','b'};
	// for ch = 1:3
	Mat temp = new Mat();
	for (int i = 0; i < 3; i++) {
	    Core.extractChannel(imRGB, temp, i);
	    temp = CalculateQuantiles(temp, satLevel);

	    Core.insertChannel(temp, imRGB, i);
	    temp.release();
	}

	// tiles = quantile(imRGB_orig(ch,:),q);
	// % [sum(imRGB_orig(ch,:)<tiles(1))/N,sum(imRGB_orig(ch,:)>tiles(2))/N] //check percentages are correct
	// imRGB(ch,:) = cbsaturate(imRGB_orig(ch,:),tiles); //saturate at the appropriate pts. in distribution

	// bottom = min(imRGB(ch,:)); top = max(imRGB(ch,:));
	// imRGB(ch,:) = (imRGB(ch,:)-bottom)*255/(top-bottom);
	//
	// print(gcf,'-dpng',[outFile '-fig' num2str(ch)])
	// end

	imRGB.convertTo(imRGB, CvType.CV_8UC3);

	return imRGB;

	// imwrite(cbunshape(imRGB,size(im_orig))/255,outFile,'png');
	// figure
	// imshow(cbunshape(imRGB,size(im_orig))/255)
	// title('Simplest Color Balance Corrected')
    }

    Mat CalculateQuantiles(Mat hist, double satLevel) { // Error handling. Works only for 1D histograms. KISS.
	// if (hist.dims > 1) return -1;
	// Scalar totalSum = Core.sumElems(hist);
	double totalSum = hist.total();
	int[] sourceBuffer = new int[(int) totalSum];
	hist.get(0, 0, sourceBuffer);
	int[] cumulativeHist = new int[512];
	// Mat cumulativeHist;
	float currentSum = 0.0f;
	// float currentSum = 0.0f;

	int vMin = 0, vMax = 255 - 1; // v1 (int) hist.total();
	// float vMin, vMax;
	// //Build the Cumulative Histogram
	for (int i = 0; i <= totalSum - 1; i++) {
	    int sourceVal = (int) sourceBuffer[i];

	    cumulativeHist[sourceVal] = sourceVal + 1;
	    // v1 currentSum += sourceBuffer[i];
	    // v1 cumulativeHist[i] = currentSum;
	}
	for (int i = 1; i <= 255; i++) {
	    cumulativeHist[i] = cumulativeHist[i] + cumulativeHist[i - 1];
	}
	// Log.i("Prob Quatile", "Sum: "+currentSum);

	// for(int i = 0; i < hist.dim[0].size; i++) {
	// currentSum += (float)hist.dim[0].at<float>(i);
	// cumulativeHist.at<float>(i) = currentSum;
	// }
	// //Calculate Lower Quantile

	// vMin = 0;
	while (cumulativeHist[vMin + 1] <= (totalSum * (satLevel / 2) / 100)) {
	    vMin = vMin + 1;
	}
	// while (cumulativeHist.at<float>(vMin + 1) <= (totalSum * (satLevel/2) / 100)){
	// vMin = vMin + 1;
	// }
	// vec.push_back(vMin);
	// //Calculate Upper Quantile

	// vMax = hist.dim[0].size-1;
	while (cumulativeHist[vMax - 1] > (totalSum * (1 - satLevel / 2) / 100)) {
	    vMax = vMax - 1;
	}
	// }
	// while (cumulativeHist.at<float>(vMax-1) > (totalSum * (1-satLevel/2) / 100)){
	// vMax = vMax - 1;
	// }
	// if(vMax < hist.dim[0].size -1){ vMax++; }
	// vec.push_back(vMax);
	// return 1;
	//
	// PP Saturate
	for (int i = 0; i < totalSum - 1; i++) {
	    if (sourceBuffer[i] <= vMin)
		sourceBuffer[i] = vMin;
	    if (sourceBuffer[i] >= vMax)
		sourceBuffer[i] = vMax;
	    // PP Rescale the pixels
	    sourceBuffer[i] = (sourceBuffer[i] - vMin) * 255 / (vMax - vMin);
	}
	hist.put(0, 0, sourceBuffer);
	return hist;
    }

    double FastArcTan(double x) {
	return (Math.PI / 4) * x - x * (Math.abs(x) - 1)
		* (0.2447 + 0.0663 * Math.abs(x));
    }
    // Color Conversion C1C2C3
    // http://www.jofcis.com/publishedpapers/2013_9_10_3783_3790.pdf
    /* On verra ça une autre fois
     * Core.extractChannel(mRgba, mR, 0);
     * Core.extractChannel(mRgba, mG, 1);
     * Core.extractChannel(mRgba, mB, 2);
     * 
     * Mat imageC1 = mR.clone();
     * Mat imageC2 = mR.clone();
     * Mat imageC3 = mR.clone();
     * 
     * Mat maxGB = new Mat();
     * Mat maxRB = new Mat();
     * Mat maxGR = new Mat();
     * Core.max(mG, mB, maxGB);
     * Core.max(mR, mB, maxRB);
     * Core.max(mG, mR, maxGR);
     * Mat divC1 = new Mat();
     * Mat divC2 = new Mat();
     * Mat divC3 = new Mat();
     * Core.divide(mR, maxGB, divC1);
     * Core.divide(mG, maxRB, divC2);
     * Core.divide(mB, maxGR, divC3);
     * 
     * // imageSrc.convertTo(imageSrc, CvType.CV_64FC1);
     * int size = (int) (mG.total() * mG.channels());
     * byte[] bC1 = new byte[size]; // use double[] instead of byte[]
     * byte[] bC2 = new byte[size];
     * byte[] bC3 = new byte[size];
     * 
     * divC1.get(0, 0, bC1);
     * divC2.get(0, 0, bC2);
     * divC3.get(0, 0, bC3);
     * 
     * for (int i = size; i-- > 0;) {
     * bC1[i] = (byte) FastArcTan(bC1[i]);
     * bC2[i] = (byte) FastArcTan(bC2[i]);
     * bC3[i] = (byte) FastArcTan(bC3[i]);
     * }
     * imageC1.put(0, 0, bC1);
     * imageC2.put(0, 0, bC2);
     * imageC3.put(0, 0, bC3);
     * 
     * // Core.insertChannel(imageC1, mRgba, 0);
     * // Core.insertChannel(imageC2, mRgba, 1);
     * // Core.insertChannel(imageC3, mRgba, 2);
     * 
     * Core.bitwise_not(mRgba, mRgba, imageC3); */

}
